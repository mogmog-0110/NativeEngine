// NativeEngine visual debugger -- a live, interactive tool in the spirit of the
// PhysX Visual Debugger: it STEPS the engine in real time (not a recording) and
// wraps it in an ImGui UI -- simulation controls, a scene hierarchy, a per-body
// inspector, live statistics, and 3D debug overlays (velocity arrows, contact
// points + normals, AABBs, selection highlight, periodic ghost images). Modern
// OpenGL 3.3 core (instanced Blinn-Phong, MSAA), right-handed / Y-up so the view
// matches NativeEngine.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

#include "gl_core.hpp"
#include "glmath.hpp"

#include "world.hpp"
#include "body.hpp"
#include "detect.hpp"

#include "imgui.h"
#include "backends/imgui_impl_glut.h"
#include "backends/imgui_impl_opengl3.h"

using namespace nv;
using namespace ne;

// ------------------------------------------------------------------ sim state
static World gWorld;
static int   gScene = 1;
static int   gStepsPerFrame = 4;    // dt = 1/240 * 4 ~ 1/60 s per frame -> real time
static bool  gPaused = false;
static int   gSelected = -1;
static double gSimTime = 0.0;
static std::mt19937_64 gRng(12345);

// visualization toggles (the debugger's "eyes")
static struct Viz {
    bool arrows = false, contacts = true, aabb = false, dimSleep = true;
    bool ghosts = false, grid = true, showCellAlways = false, joints = true;
    float velScale = 0.25f;
} gViz;

// An interactive WASD capsule character (scene 15). Not a rigid body: moved by
// collide-and-slide against the world, rendered separately.
static struct Character {
    bool active = false;
    V3 pos; V3 vel;
    double r = 0.4, hh = 0.6;
    bool grounded = false;
} gChar;
static bool gKey[256] = {false};

// camera / window
static float gAz = 0.7f, gEl = 0.45f, gDist = 40.0f;
static Vec3  gCenter{0, 0, 0};
static int   gLastX = 0, gLastY = 0; static bool gDrag = false;
static int   gWinW = 1440, gWinH = 900;
static bool  gShowCell = false;

static GLuint gProg = 0;
static GLint  gUVP = -1, gUCam = -1, gULight = -1, gUUnlit = -1;
static GLuint gLineProg = 0; static GLint gLUVP = -1;
static GLuint gLineVao = 0, gLineVbo = 0;
static std::vector<float> gLines;   // x,y,z,r,g,b per vertex, GL_LINES

struct Mesh {
    GLuint vao = 0, vbo = 0, ebo = 0, inst = 0;
    GLsizei indexCount = 0; GLenum prim = GL_TRIANGLES;
    std::vector<float> instData; int instCount = 0;
};
static Mesh gSphere, gCyl, gBox, gWire;

// ------------------------------------------------------------------ shaders
static const char* kVert =
    "#version 330 core\n"
    "layout(location=0) in vec3 aPos; layout(location=1) in vec3 aNormal;\n"
    "layout(location=3) in vec4 m0; layout(location=4) in vec4 m1;\n"
    "layout(location=5) in vec4 m2; layout(location=6) in vec4 m3;\n"
    "layout(location=7) in vec3 aColor;\n"
    "uniform mat4 uVP; out vec3 vN; out vec3 vW; out vec3 vC;\n"
    "void main(){ mat4 M=mat4(m0,m1,m2,m3); vec4 w=M*vec4(aPos,1.0);\n"
    "  vW=w.xyz; vN=mat3(M)*aNormal; vC=aColor; gl_Position=uVP*w; }\n";
static const char* kFrag =
    "#version 330 core\n"
    "in vec3 vN; in vec3 vW; in vec3 vC; out vec4 o;\n"
    "uniform vec3 uCam; uniform vec3 uLight; uniform int uUnlit;\n"
    "void main(){ if(uUnlit==1){ o=vec4(vC,1.0); return; }\n"
    "  vec3 N=normalize(vN); vec3 L=normalize(-uLight); vec3 V=normalize(uCam-vW);\n"
    "  vec3 H=normalize(L+V); float diff=max(dot(N,L),0.0);\n"
    "  float spec=pow(max(dot(N,H),0.0),48.0); float hemi=0.5+0.5*N.y;\n"
    "  vec3 amb=mix(vec3(0.16,0.17,0.20),vec3(0.34,0.36,0.40),hemi);\n"
    "  vec3 c=vC*(amb+diff*0.9)+vec3(1.0)*spec*0.35;\n"
    "  c=c/(c+vec3(1.0)); c=pow(c,vec3(1.0/2.2)); o=vec4(c,1.0); }\n";
static const char* kLineVert =
    "#version 330 core\n"
    "layout(location=0) in vec3 aPos; layout(location=1) in vec3 aCol;\n"
    "uniform mat4 uVP; out vec3 vCol;\n"
    "void main(){ vCol=aCol; gl_Position=uVP*vec4(aPos,1.0); }\n";
static const char* kLineFrag =
    "#version 330 core\n"
    "in vec3 vCol; out vec4 o; void main(){ o=vec4(vCol,1.0); }\n";

static GLuint compile(GLenum t, const char* src) {
    GLuint s = glCreateShader(t); glShaderSource(s, 1, &src, nullptr); glCompileShader(s);
    GLint ok = 0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) { char log[1024]; glGetShaderInfoLog(s, 1024, nullptr, log); std::fprintf(stderr, "shader: %s\n", log); }
    return s;
}
static GLuint link(const char* v, const char* f) {
    GLuint p = glCreateProgram(); GLuint vs = compile(GL_VERTEX_SHADER, v), fs = compile(GL_FRAGMENT_SHADER, f);
    glAttachShader(p, vs); glAttachShader(p, fs); glLinkProgram(p);
    GLint ok = 0; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) { char log[1024]; glGetProgramInfoLog(p, 1024, nullptr, log); std::fprintf(stderr, "link: %s\n", log); }
    glDeleteShader(vs); glDeleteShader(fs); return p;
}

// ------------------------------------------------------------------ meshes
static void upload(Mesh& m, const std::vector<float>& verts, const std::vector<unsigned>& idx, GLenum prim) {
    m.prim = prim; m.indexCount = (GLsizei)idx.size();
    glGenVertexArrays(1, &m.vao); glBindVertexArray(m.vao);
    glGenBuffers(1, &m.vbo); glBindBuffer(GL_ARRAY_BUFFER, m.vbo);
    glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_STATIC_DRAW);
    glGenBuffers(1, &m.ebo); glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(unsigned), idx.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0); glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float))); glEnableVertexAttribArray(1);
    glGenBuffers(1, &m.inst); glBindBuffer(GL_ARRAY_BUFFER, m.inst);
    const GLsizei stride = 19 * sizeof(float);
    for (int i = 0; i < 4; ++i) {
        glVertexAttribPointer(3 + i, 4, GL_FLOAT, GL_FALSE, stride, (void*)(size_t)(i * 4 * sizeof(float)));
        glEnableVertexAttribArray(3 + i); glVertexAttribDivisor(3 + i, 1);
    }
    glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, stride, (void*)(size_t)(16 * sizeof(float)));
    glEnableVertexAttribArray(7); glVertexAttribDivisor(7, 1);
    glBindVertexArray(0);
}
static void makeSphere(Mesh& m, int stacks = 20, int slices = 28) {
    std::vector<float> v; std::vector<unsigned> idx; const float PI = 3.14159265f;
    for (int i = 0; i <= stacks; ++i) { float phi = i / (float)stacks * PI;
        for (int j = 0; j <= slices; ++j) { float th = j / (float)slices * 2 * PI;
            float x = std::sin(phi) * std::cos(th), y = std::cos(phi), z = std::sin(phi) * std::sin(th);
            v.insert(v.end(), {x, y, z, x, y, z}); } }
    int row = slices + 1;
    for (int i = 0; i < stacks; ++i) for (int j = 0; j < slices; ++j) {
        unsigned a = i * row + j, b = a + row; idx.insert(idx.end(), {a, b, a + 1u, a + 1u, b, b + 1u}); }
    upload(m, v, idx, GL_TRIANGLES);
}
static void makeCylinder(Mesh& m, int slices = 28) {
    std::vector<float> v; std::vector<unsigned> idx; const float PI = 3.14159265f;
    for (int j = 0; j <= slices; ++j) { float th = j / (float)slices * 2 * PI, x = std::cos(th), z = std::sin(th);
        v.insert(v.end(), {x, -1, z, x, 0, z}); v.insert(v.end(), {x, 1, z, x, 0, z}); }
    for (int j = 0; j < slices; ++j) { unsigned a = j * 2u; idx.insert(idx.end(), {a, a + 1u, a + 2u, a + 2u, a + 1u, a + 3u}); }
    for (int sgn = -1; sgn <= 1; sgn += 2) { unsigned base = (unsigned)(v.size() / 6); float y = (float)sgn;
        v.insert(v.end(), {0, y, 0, 0, y, 0});
        for (int j = 0; j <= slices; ++j) { float th = j / (float)slices * 2 * PI; v.insert(v.end(), {std::cos(th), y, std::sin(th), 0, y, 0}); }
        for (int j = 0; j < slices; ++j) { if (sgn < 0) idx.insert(idx.end(), {base, base + 1u + j, base + 2u + j}); else idx.insert(idx.end(), {base, base + 2u + j, base + 1u + j}); } }
    upload(m, v, idx, GL_TRIANGLES);
}
static void makeBox(Mesh& m) {
    const float f[6][6] = {{1,0,0, 1,0,0},{-1,0,0, -1,0,0},{0,1,0, 0,1,0},{0,-1,0, 0,-1,0},{0,0,1, 0,0,1},{0,0,-1, 0,0,-1}};
    std::vector<float> v; std::vector<unsigned> idx;
    for (int face = 0; face < 6; ++face) {
        Vec3 n{f[face][0], f[face][1], f[face][2]};
        Vec3 u = std::fabs(n.y) > 0.5f ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        Vec3 t = n.cross(u); u = t.cross(n);
        Vec3 c = n;
        Vec3 corners[4] = {c + u + t, c - u + t, c - u - t, c + u - t};
        unsigned base = (unsigned)(v.size() / 6);
        for (auto& cc : corners) v.insert(v.end(), {cc.x, cc.y, cc.z, n.x, n.y, n.z});
        idx.insert(idx.end(), {base, base + 1u, base + 2u, base, base + 2u, base + 3u});
    }
    upload(m, v, idx, GL_TRIANGLES);
}
static void makeWireBox(Mesh& m) {
    std::vector<float> v;
    for (int i = 0; i < 8; ++i) { float x = (i & 1) ? 1.f : -1.f, y = (i & 2) ? 1.f : -1.f, z = (i & 4) ? 1.f : -1.f; v.insert(v.end(), {x, y, z, 0, 1, 0}); }
    std::vector<unsigned> idx = {0,1, 1,3, 3,2, 2,0, 4,5, 5,7, 7,6, 6,4, 0,4, 1,5, 2,6, 3,7};
    upload(m, v, idx, GL_LINES);
}

static void pushInstance(Mesh& m, const Mat4& model, const Vec3& col) {
    for (int i = 0; i < 16; ++i) m.instData.push_back(model.m[i]);
    m.instData.push_back(col.x); m.instData.push_back(col.y); m.instData.push_back(col.z); ++m.instCount;
}
static void drawMesh(Mesh& m, int unlit) {
    if (m.instCount == 0) return;
    glBindVertexArray(m.vao); glBindBuffer(GL_ARRAY_BUFFER, m.inst);
    glBufferData(GL_ARRAY_BUFFER, m.instData.size() * sizeof(float), m.instData.data(), GL_DYNAMIC_DRAW);
    glUniform1i(gUUnlit, unlit);
    glDrawElementsInstanced(m.prim, m.indexCount, GL_UNSIGNED_INT, 0, m.instCount);
    glBindVertexArray(0);
}
static Vec3 palette(int id, Shape s, bool sleeping, bool isStat) {
    if (isStat) return {0.30f, 0.31f, 0.34f};
    float t = (float)((id * 2654435761u) % 997) / 997.0f;
    Vec3 c;
    if (s == Shape::Sphere) c = {0.90f, 0.48f + 0.22f * t, 0.26f};
    else if (s == Shape::Box) c = {0.35f + 0.2f * t, 0.62f, 0.42f + 0.2f * t};
    else c = {0.28f + 0.15f * t, 0.55f + 0.2f * t, 0.80f};
    if (sleeping && gViz.dimSleep) c = c * 0.55f;
    return c;
}

// ------------------------------------------------------------------ debug lines
static void line(const Vec3& a, const Vec3& b, const Vec3& c) {
    gLines.insert(gLines.end(), {a.x, a.y, a.z, c.x, c.y, c.z, b.x, b.y, b.z, c.x, c.y, c.z});
}
static void cross3(const Vec3& p, float r, const Vec3& c) {
    line({p.x - r, p.y, p.z}, {p.x + r, p.y, p.z}, c);
    line({p.x, p.y - r, p.z}, {p.x, p.y + r, p.z}, c);
    line({p.x, p.y, p.z - r}, {p.x, p.y, p.z + r}, c);
}
static void aabb(const Vec3& ctr, float h, const Vec3& c) {
    Vec3 v[8];
    for (int i = 0; i < 8; ++i) v[i] = {ctr.x + ((i & 1) ? h : -h), ctr.y + ((i & 2) ? h : -h), ctr.z + ((i & 4) ? h : -h)};
    int e[12][2] = {{0,1},{1,3},{3,2},{2,0},{4,5},{5,7},{7,6},{6,4},{0,4},{1,5},{2,6},{3,7}};
    for (auto& ed : e) line(v[ed[0]], v[ed[1]], c);
}
static Vec3 v3(const V3& d) { return {(float)d.x, (float)d.y, (float)d.z}; }

// ------------------------------------------------------------------ scenes
static int gNextId = 0;
static const int kNumScenes = 15;
static double urand(double a, double b) { std::uniform_real_distribution<double> d(a, b); return d(gRng); }

static void addGround(double half = 30.0) {
    Body f = makeBox(gNextId++, {0, -1.0, 0}, Q{1, 0, 0, 0}, {half, 1.0, half}, 1.0);
    f.invMass = 0; f.invInertiaBody = {0, 0, 0}; f.dynamic = false;
    gWorld.bodies.push_back(f);
}
static std::size_t gKinBody = (std::size_t)-1;   // kinematic platform index (scene 12)

// Joint builders operating directly on gWorld (world anchors, world axis).
static void addBallJ(std::size_t a, std::size_t b, const V3& w) {
    Joint j; j.type = JointType::Ball; j.a = a; j.b = b;
    j.localA = gWorld.bodies[a].q.inverseRotate(w - gWorld.bodies[a].x);
    j.localB = gWorld.bodies[b].q.inverseRotate(w - gWorld.bodies[b].x);
    gWorld.joints.push_back(j);
}
static void addHingeJ(std::size_t a, std::size_t b, const V3& w, const V3& axis) {
    Joint j; j.type = JointType::Hinge; j.a = a; j.b = b;
    j.localA = gWorld.bodies[a].q.inverseRotate(w - gWorld.bodies[a].x);
    j.localB = gWorld.bodies[b].q.inverseRotate(w - gWorld.bodies[b].x);
    V3 ax = axis.normalized();
    j.axisA = gWorld.bodies[a].q.inverseRotate(ax);
    j.axisB = gWorld.bodies[b].q.inverseRotate(ax);
    j.refRel = gWorld.bodies[a].q * gWorld.bodies[b].q.conjugate();
    gWorld.joints.push_back(j);
}
static Q qAxisAngle(const V3& axis, double ang) {
    V3 a = axis.normalized(); double s = std::sin(ang * 0.5);
    return Q{std::cos(ang * 0.5), a.x * s, a.y * s, a.z * s};
}
static Body makeStaticBox(const V3& pos, const Q& rot, const V3& he) {
    Body b = makeBox(gNextId++, pos, rot, he, 1.0);
    b.invMass = 0; b.invInertiaBody = {0, 0, 0}; b.dynamic = false;
    return b;
}

static void openGravityWorld() {
    gWorld.openBoundary = true; gWorld.gravity = {0, -10, 0};
    gWorld.friction = 0.6; gWorld.restitution = 0.0; gWorld.sleepEnabled = true;
    gWorld.box.half = 40;
}

static void buildScene(int id) {
    gWorld = World{};
    gWorld.dt = 1.0 / 240.0;
    gNextId = 0; gScene = id; gShowCell = false; gSelected = -1; gSimTime = 0.0;
    gChar.active = false; gKinBody = (std::size_t)-1;
    gRng.seed(12345);

    switch (id) {
        case 1: {  // box pyramid on the ground (classic HelloWorld)
            openGravityWorld(); addGround();
            for (int level = 0; level < 6; ++level) {
                int n = 6 - level;
                for (int i = 0; i < n; ++i)
                    gWorld.bodies.push_back(makeBox(gNextId++,
                        {(i - (n - 1) / 2.0) * 2.05, 1.0 + level * 2.0, 0},
                        Q{1, 0, 0, 0}, {1, 1, 1}, 1.0));
            }
            break;
        }
        case 2: {  // brick wall knocked down by a projectile
            openGravityWorld(); gWorld.sleepEnabled = false; addGround();
            const int cols = 8, rows = 6;
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c) {
                    double off = (r % 2) ? 1.0 : 0.0;
                    gWorld.bodies.push_back(makeBox(gNextId++,
                        {(c - cols / 2.0) * 2.0 + off, 1.0 + r * 2.0, 0},
                        Q{1, 0, 0, 0}, {1.0, 1.0, 2.0}, 1.0));
                }
            Body ball = makeSphere(gNextId++, {0, 6, 22}, 1.6, 6.0);
            ball.v = {0, 0, -35};
            gWorld.bodies.push_back(ball);
            break;
        }
        case 3: {  // dominoes
            openGravityWorld(); addGround();
            const int N = 16;
            for (int i = 0; i < N; ++i)
                gWorld.bodies.push_back(makeBox(gNextId++,
                    {(i - N / 2.0) * 1.5, 2.0, 0}, Q{1, 0, 0, 0}, {0.25, 2.0, 1.2}, 1.0));
            gWorld.bodies[1].w = {0, 0, -3.0};
            gWorld.bodies[1].v = {2.0, 0, 0};
            break;
        }
        case 4: {  // sphere pile
            openGravityWorld(); gWorld.restitution = 0.1; gWorld.friction = 0.4; addGround();
            for (int i = 0; i < 140; ++i)
                gWorld.bodies.push_back(makeSphere(gNextId++, {urand(-6, 6), urand(3, 20), urand(-6, 6)}, 1.0, 1.0));
            break;
        }
        case 5: {  // cylinder pile
            openGravityWorld(); gWorld.restitution = 0.05; addGround();
            for (int i = 0; i < 60; ++i) {
                Q q{urand(-1, 1), urand(-1, 1), urand(-1, 1), urand(-1, 1)}; q = q.normalized();
                gWorld.bodies.push_back(makeCylinder(gNextId++, {urand(-6, 6), urand(3, 20), urand(-6, 6)}, q, 1.0, 4.0, 1.0));
            }
            break;
        }
        case 6: {  // hanging rope/chain (rigid joints)
            gWorld.openBoundary = true; gWorld.gravity = {0, -10, 0}; gWorld.box.half = 24;
            Body anchor = makeSphere(gNextId++, {0, 14, 0}, 0.3, 1.0);
            anchor.invMass = 0; anchor.invInertiaBody = {0, 0, 0}; anchor.dynamic = false;
            gWorld.bodies.push_back(anchor);
            const int N = 12; const double L = 1.4;
            for (int i = 1; i <= N; ++i) gWorld.bodies.push_back(makeSphere(gNextId++, {i * L, 14, 0}, 0.5, 1.0));
            for (int i = 0; i < N; ++i) { DistanceJoint j; j.a = i; j.b = i + 1; j.rest = L; gWorld.distanceJoints.push_back(j); }
            break;
        }
        case 7: {  // elastic sphere gas in a reflective BOX
            gWorld.box.half = 10; gWorld.restitution = 1.0; gShowCell = true;
            for (int i = 0; i < 130; ++i) {
                Body b = makeSphere(gNextId++, {urand(-8, 8), urand(-8, 8), urand(-8, 8)}, 0.8, 1.0);
                b.v = {urand(-5, 5), urand(-5, 5), urand(-5, 5)}; gWorld.bodies.push_back(b);
            }
            break;
        }
        case 8: {  // PERIODIC gas
            gWorld.box.half = 9; gWorld.box.periodic = true; gWorld.restitution = 1.0; gShowCell = true;
            for (int i = 0; i < 160; ++i) {
                Body b = makeSphere(gNextId++, {urand(-8.5, 8.5), urand(-8.5, 8.5), urand(-8.5, 8.5)}, 0.8, 1.0);
                b.v = {urand(-4, 4), urand(-4, 4), urand(-4, 4)}; gWorld.bodies.push_back(b);
            }
            // Zero the net momentum so the whole cloud does not drift + wrap as one
            // (standard MD practice: no centre-of-mass motion). PBC conserves the
            // total momentum, so an un-zeroed cloud streams across the box forever,
            // which reads as "everything drifts off" even though it is correct.
            V3 psum; double M = 0;
            for (const Body& b : gWorld.bodies) { double m = 1.0 / b.invMass; psum += b.v * m; M += m; }
            V3 vcm = psum * (1.0 / M);
            for (Body& b : gWorld.bodies) b.v -= vcm;
            break;
        }
        case 9: {  // two spheres colliding ACROSS the periodic boundary
            gWorld.box.half = 6; gWorld.box.periodic = true; gWorld.restitution = 1.0; gShowCell = true;
            Body a = makeSphere(gNextId++, {5.2, 0, 0}, 1.0, 1.0); a.v = {2, 0, 0};
            Body b = makeSphere(gNextId++, {-5.2, 0.3, 0}, 1.0, 1.0); b.v = {-2, 0, 0};
            gWorld.bodies.push_back(a); gWorld.bodies.push_back(b);
            break;
        }
        case 10: {  // ragdoll: capsules + a head linked by ball joints, dropped
            openGravityWorld(); addGround(); gWorld.sleepEnabled = false;
            std::size_t g = gWorld.bodies.size() - 1; (void)g;
            V3 c{0, 8, 0};
            std::size_t torso = gWorld.bodies.size();
            gWorld.bodies.push_back(makeCapsule(gNextId++, c, Q{1, 0, 0, 0}, 0.35, 0.9, 1.0));
            std::size_t head = gWorld.bodies.size();
            gWorld.bodies.push_back(makeSphere(gNextId++, c + V3{0, 1.6, 0}, 0.45, 1.0));
            addBallJ(torso, head, c + V3{0, 1.25, 0});
            // arms (along X) and legs (down), each a capsule on a ball joint
            for (int s = -1; s <= 1; s += 2) {
                std::size_t arm = gWorld.bodies.size();
                gWorld.bodies.push_back(makeCapsule(gNextId++, c + V3{s * 1.1, 0.9, 0},
                                        qAxisAngle({0, 0, 1}, PI / 2), 0.22, 0.7, 1.0));
                addBallJ(torso, arm, c + V3{s * 0.45, 1.1, 0});
                std::size_t leg = gWorld.bodies.size();
                gWorld.bodies.push_back(makeCapsule(gNextId++, c + V3{s * 0.4, -2.0, 0}, Q{1, 0, 0, 0}, 0.26, 0.9, 1.0));
                addBallJ(torso, leg, c + V3{s * 0.4, -1.0, 0});
            }
            break;
        }
        case 11: {  // heightfield terrain with objects tumbling down it
            openGravityWorld(); gWorld.restitution = 0.1;
            const int N = 24; const double sp = 2.0;
            std::vector<V3> verts; std::vector<int> idx;
            std::mt19937_64 rng(7);
            std::uniform_real_distribution<double> hh(0.0, 2.2);
            for (int r = 0; r < N; ++r)
                for (int c = 0; c < N; ++c) {
                    double d = std::sqrt((r - N / 2.0) * (r - N / 2.0) + (c - N / 2.0) * (c - N / 2.0));
                    double h = std::sin(r * 0.5) * std::cos(c * 0.5) * 1.2 + (d > 8 ? (d - 8) * 0.6 : 0) + hh(rng) * 0.2;
                    verts.push_back(V3{(c - N / 2.0) * sp, h, (r - N / 2.0) * sp});
                }
            for (int r = 0; r + 1 < N; ++r)
                for (int c = 0; c + 1 < N; ++c) {
                    int i00 = r * N + c, i01 = i00 + 1, i10 = i00 + N, i11 = i10 + 1;
                    idx.insert(idx.end(), {i00, i10, i11, i00, i11, i01});
                }
            auto md = std::make_shared<MeshData>(); md->verts = verts;
            for (std::size_t i = 0; i + 2 < idx.size(); i += 3) md->tris.push_back({idx[i], idx[i + 1], idx[i + 2]});
            md->build();
            gWorld.bodies.push_back(makeMesh(gNextId++, {0, 0, 0}, Q{1, 0, 0, 0}, md));
            for (int i = 0; i < 40; ++i) {
                int k = i % 3;
                if (k == 0) gWorld.bodies.push_back(makeSphere(gNextId++, {urand(-6, 6), urand(8, 16), urand(-6, 6)}, 0.6, 1.0));
                else if (k == 1) gWorld.bodies.push_back(makeBox(gNextId++, {urand(-6, 6), urand(8, 16), urand(-6, 6)}, Q{1, 0, 0, 0}, {0.6, 0.6, 0.6}, 1.0));
                else gWorld.bodies.push_back(makeCapsule(gNextId++, {urand(-6, 6), urand(8, 16), urand(-6, 6)}, Q{1, 0, 0, 0}, 0.4, 0.5, 1.0));
            }
            gWorld.box.half = 60;
            break;
        }
        case 12: {  // kinematic moving platform carrying a stack of boxes
            openGravityWorld(); addGround(); gWorld.friction = 0.9;
            std::size_t plat = gWorld.bodies.size();
            Body p = makeBox(gNextId++, {0, 1.5, 0}, Q{1, 0, 0, 0}, {3, 0.3, 3}, 1.0);
            p.invMass = 0; p.invInertiaBody = {0, 0, 0}; p.dynamic = false; p.kinematic = true;
            gWorld.bodies.push_back(p);
            gKinBody = plat;
            for (int i = 0; i < 4; ++i)
                gWorld.bodies.push_back(makeBox(gNextId++, {0, 2.3 + i * 1.05, 0}, Q{1, 0, 0, 0}, {0.5, 0.5, 0.5}, 1.0));
            break;
        }
        case 13: {  // suspension bridge of planks linked by hinges
            openGravityWorld(); gWorld.sleepEnabled = false;
            const int N = 12; const double pw = 1.0;
            std::size_t left = gWorld.bodies.size();
            Body la = makeBox(gNextId++, {-(N * pw) / 2 - 1, 6, 0}, Q{1, 0, 0, 0}, {1, 1, 2}, 1.0);
            la.invMass = 0; la.invInertiaBody = {0, 0, 0}; la.dynamic = false;
            gWorld.bodies.push_back(la);
            std::size_t prev = left;
            for (int i = 0; i < N; ++i) {
                std::size_t plank = gWorld.bodies.size();
                double x = -(N * pw) / 2 + pw * 0.5 + i * pw;
                gWorld.bodies.push_back(makeBox(gNextId++, {x, 6, 0}, Q{1, 0, 0, 0}, {pw * 0.5, 0.12, 2}, 1.0));
                addHingeJ(prev, plank, {x - pw * 0.5, 6, 0}, {0, 0, 1});
                prev = plank;
            }
            Body ra = makeBox(gNextId++, {(N * pw) / 2 + 1, 6, 0}, Q{1, 0, 0, 0}, {1, 1, 2}, 1.0);
            ra.invMass = 0; ra.invInertiaBody = {0, 0, 0}; ra.dynamic = false;
            std::size_t right = gWorld.bodies.size(); gWorld.bodies.push_back(ra);
            addHingeJ(prev, right, {(N * pw) / 2, 6, 0}, {0, 0, 1});
            for (int i = 0; i < 6; ++i)
                gWorld.bodies.push_back(makeSphere(gNextId++, {urand(-3, 3), 10.0 + i, 0}, 0.5, 2.0));
            gWorld.box.half = 30;
            break;
        }
        case 14: {  // compound bodies (dumbbells + tables) tumbling
            openGravityWorld(); addGround();
            for (int i = 0; i < 6; ++i) {
                std::vector<ChildShape> kids = {ChildShape::sphere({-1.1, 0, 0}, 0.6),
                                                ChildShape::sphere({1.1, 0, 0}, 0.6),
                                                ChildShape::cylinder({0, 0, 0}, qAxisAngle({0, 0, 1}, PI / 2), 0.2, 1.0)};
                gWorld.bodies.push_back(makeCompound(gNextId++, {urand(-5, 5), 4.0 + i * 2.5, urand(-5, 5)},
                                        qAxisAngle({urand(-1, 1), urand(-1, 1), urand(-1, 1)}, urand(0, 3)), kids, 1.0));
            }
            for (int i = 0; i < 4; ++i) {   // little tables
                std::vector<ChildShape> t = {ChildShape::box({0, 0.6, 0}, Q{1, 0, 0, 0}, {1.0, 0.15, 1.0})};
                for (int c = 0; c < 4; ++c) {
                    double sx = (c & 1) ? 0.8 : -0.8, sz = (c & 2) ? 0.8 : -0.8;
                    t.push_back(ChildShape::box({sx, 0, sz}, Q{1, 0, 0, 0}, {0.12, 0.6, 0.12}));
                }
                gWorld.bodies.push_back(makeCompound(gNextId++, {urand(-5, 5), 10.0 + i * 2.0, urand(-5, 5)}, Q{1, 0, 0, 0}, t, 1.0));
            }
            break;
        }
        case 15: {  // WASD capsule character on a box level (ramps + steps)
            openGravityWorld();
            gWorld.bodies.push_back(makeStaticBox({0, -0.5, 0}, Q{1, 0, 0, 0}, {25, 0.5, 25}));
            // a ramp
            gWorld.bodies.push_back(makeStaticBox({8, 1.2, 0}, qAxisAngle({0, 0, 1}, -0.5), {5, 0.3, 4}));
            // steps
            for (int i = 0; i < 5; ++i)
                gWorld.bodies.push_back(makeStaticBox({-8.0, 0.4 + i * 0.4, i * 1.0 - 2}, Q{1, 0, 0, 0}, {2, 0.2 + i * 0.2, 0.5}));
            // some walls / boxes to bump into
            gWorld.bodies.push_back(makeStaticBox({0, 1.5, 12}, Q{1, 0, 0, 0}, {12, 2, 0.5}));
            gWorld.bodies.push_back(makeStaticBox({0, 1.5, -12}, Q{1, 0, 0, 0}, {12, 2, 0.5}));
            for (int i = 0; i < 6; ++i)
                gWorld.bodies.push_back(makeBox(gNextId++, {urand(-4, 4), 1.0 + i, urand(-4, 4)}, Q{1, 0, 0, 0}, {0.5, 0.5, 0.5}, 1.0));
            gChar.active = true; gChar.pos = {0, 2.0, 0}; gChar.vel = {0, 0, 0};
            gWorld.box.half = 30;
            break;
        }
    }
    gDist = gWorld.box.half * (gShowCell ? 2.8f : 1.4f);
    if (gChar.active) gDist = 14.0f;
    gCenter = {0, gShowCell ? 0.0f : 5.0f, 0};
    // In a periodic scene, default the wrap-around ghosts ON so a body leaving one
    // face is visibly re-entering the opposite one -- otherwise the wrap looks like
    // a teleport and cross-boundary collisions look like bouncing off nothing.
    gViz.ghosts = gWorld.box.periodic;
}
static const char* kSceneNames[kNumScenes] = {
    "1 box pyramid", "2 brick wall + projectile", "3 dominoes", "4 sphere pile",
    "5 cylinder pile", "6 rope", "7 elastic gas (box)", "8 PERIODIC gas", "9 PERIODIC pair",
    "10 ragdoll (joints)", "11 terrain (heightfield)", "12 kinematic platform",
    "13 hinge bridge", "14 compound bodies", "15 CHARACTER (WASD)"};
static void dropSphere() {
    Body b = makeSphere(gNextId++, {urand(-4, 4), 20.0, urand(-4, 4)}, 1.2, 1.0);
    gWorld.bodies.push_back(b);
}

static bool isStatic(const Body& b) { return !b.dynamic && b.invMass == 0.0; }

// ------------------------------------------------------------------ picking
static void cameraRay(int mx, int my, Vec3& origin, Vec3& dir) {
    float aspect = gWinW / (float)(gWinH ? gWinH : 1);
    Vec3 eye{gCenter.x + gDist * std::cos(gEl) * std::sin(gAz),
             gCenter.y + gDist * std::sin(gEl),
             gCenter.z + gDist * std::cos(gEl) * std::cos(gAz)};
    Vec3 fwd = (gCenter - eye).norm();
    Vec3 right = fwd.cross({0, 1, 0}).norm();
    Vec3 up = right.cross(fwd);
    float tanHalf = std::tan(45.0f * 3.14159f / 180.0f * 0.5f);
    float ndcx = (2.0f * (mx + 0.5f) / gWinW - 1.0f) * tanHalf * aspect;
    float ndcy = (1.0f - 2.0f * (my + 0.5f) / gWinH) * tanHalf;
    origin = eye;
    dir = (fwd + right * ndcx + up * ndcy).norm();
}
static void pick(int mx, int my) {
    Vec3 o, d; cameraRay(mx, my, o, d);
    int best = -1; float bestT = 1e30f;
    for (size_t i = 0; i < gWorld.bodies.size(); ++i) {
        const Body& b = gWorld.bodies[i];
        Vec3 ctr = v3(b.x); float r = (float)b.boundingRadius();
        Vec3 m = o - ctr; float bb = m.dot(d), cc = m.dot(m) - r * r;
        if (cc > 0 && bb > 0) continue;
        float disc = bb * bb - cc; if (disc < 0) continue;
        float t = -bb - std::sqrt(disc); if (t < 0) t = 0;
        if (t < bestT) { bestT = t; best = (int)i; }
    }
    gSelected = best;   // -1 (empty click) deselects
}

// ------------------------------------------------------------------ UI
static float gKE[180] = {0}; static int gKEHead = 0;
static void pushKE(float e) { gKE[gKEHead] = e; gKEHead = (gKEHead + 1) % 180; }

static void buildUI() {
    ImGuiIO& io = ImGui::GetIO();
    // must cover every Shape enum value (out-of-range indexing crashes snprintf)
    static const char* shapeName[] = {"sphere", "cylinder", "box", "capsule",
                                      "convex", "plane", "mesh", "compound"};

    // ---- Simulation -------------------------------------------------------
    ImGui::SetNextWindowPos({10, 10}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize({330, 300}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Simulation");
    int sceneIdx = gScene - 1;
    if (ImGui::Combo("scene", &sceneIdx, kSceneNames, kNumScenes)) buildScene(sceneIdx + 1);
    if (ImGui::Button(gPaused ? "Play" : "Pause")) gPaused = !gPaused;
    ImGui::SameLine(); if (ImGui::Button("Step")) { gWorld.captureContacts = gViz.contacts; gWorld.step(); gSimTime += gWorld.dt; }
    ImGui::SameLine(); if (ImGui::Button("Reset")) buildScene(gScene);
    ImGui::SameLine(); if (ImGui::Button("Drop")) dropSphere();
    ImGui::SliderInt("substeps/frame", &gStepsPerFrame, 1, 16);
    ImGui::Text("dt = 1/%.0f s    sim t = %.2f s", 1.0 / gWorld.dt, gSimTime);
    ImGui::Separator();
    bool grav = (gWorld.gravity.y != 0.0);
    if (ImGui::Checkbox("gravity", &grav)) gWorld.gravity.y = grav ? -10.0 : 0.0;
    float rest = (float)gWorld.restitution, fric = (float)gWorld.friction;
    if (ImGui::SliderFloat("restitution", &rest, 0.0f, 1.0f)) gWorld.restitution = rest;
    if (ImGui::SliderFloat("friction", &fric, 0.0f, 1.5f)) gWorld.friction = fric;
    bool sleep = gWorld.sleepEnabled;
    if (ImGui::Checkbox("sleeping", &sleep)) gWorld.sleepEnabled = sleep;
    ImGui::TextDisabled("space pause | n next | r reset | d drop | h camera");
    ImGui::End();

    // ---- Visualization ----------------------------------------------------
    ImGui::SetNextWindowPos({10, 320}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize({330, 210}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Visualization");
    ImGui::Checkbox("velocity arrows", &gViz.arrows);
    ImGui::SameLine(); ImGui::SetNextItemWidth(90); ImGui::SliderFloat("scale", &gViz.velScale, 0.02f, 1.0f);
    ImGui::Checkbox("contact points + normals", &gViz.contacts);
    ImGui::Checkbox("bounding boxes", &gViz.aabb);
    ImGui::Checkbox("dim sleeping bodies", &gViz.dimSleep);
    ImGui::Checkbox("ground grid", &gViz.grid);
    ImGui::SameLine(); ImGui::Checkbox("joints", &gViz.joints);
    if (gChar.active) ImGui::TextDisabled("CHARACTER: W/A/S/D move, E jump (camera-relative)");
    if (gWorld.box.periodic) ImGui::Checkbox("periodic ghost images", &gViz.ghosts);
    else ImGui::Checkbox("always show domain box", &gViz.showCellAlways);
    ImGui::TextDisabled("click a body to select   yellow=contact  cyan=normal");
    ImGui::End();

    // ---- Statistics -------------------------------------------------------
    int awake = 0, sleeping = 0, statics = 0;
    for (const Body& b : gWorld.bodies) {
        if (isStatic(b)) ++statics; else if (b.sleeping) ++sleeping; else ++awake;
    }
    double ke = gWorld.totalKinetic();
    V3 P = gWorld.totalLinearMomentum(), L = gWorld.totalAngularMomentum();
    ImGui::SetNextWindowPos({(float)gWinW - 340, 10}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize({330, 250}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Statistics");
    ImGui::Text("%.1f FPS   (%.2f ms/frame)", io.Framerate, 1000.0f / io.Framerate);
    ImGui::Separator();
    ImGui::Text("bodies    : %zu", gWorld.bodies.size());
    ImGui::Text("awake %d   sleeping %d   static %d", awake, sleeping, statics);
    ImGui::Text("contacts  : %zu", gWorld.debugContacts.size());
    ImGui::Separator();
    ImGui::Text("kinetic E : %.3f", ke);
    ImGui::Text("|momentum|: %.3f", std::sqrt(P.dot(P)));
    ImGui::Text("|ang.mom.|: %.3f", std::sqrt(L.dot(L)));
    ImGui::PlotLines("KE", gKE, 180, gKEHead, nullptr, 0.0f, FLT_MAX, {310, 60});
    ImGui::End();

    // ---- Scene hierarchy --------------------------------------------------
    ImGui::SetNextWindowPos({(float)gWinW - 340, 270}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize({330, 300}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Scene Hierarchy");
    ImGui::Text("%zu bodies", gWorld.bodies.size());
    ImGui::BeginChild("list");
    for (size_t i = 0; i < gWorld.bodies.size(); ++i) {
        const Body& b = gWorld.bodies[i];
        char lbl[96];
        std::snprintf(lbl, sizeof lbl, "#%d  %s%s", b.id, shapeName[(int)b.shape],
                      isStatic(b) ? "  [static]" : (b.sleeping ? "  [sleep]" : ""));
        if (ImGui::Selectable(lbl, gSelected == (int)i)) gSelected = (int)i;
    }
    ImGui::EndChild();
    ImGui::End();

    // ---- Inspector --------------------------------------------------------
    ImGui::SetNextWindowPos({(float)gWinW - 340, 580}, ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize({330, 300}, ImGuiCond_FirstUseEver);
    ImGui::Begin("Inspector");
    if (gSelected >= 0 && gSelected < (int)gWorld.bodies.size()) {
        Body& b = gWorld.bodies[gSelected];
        ImGui::Text("body #%d  (%s)", b.id, shapeName[(int)b.shape]);
        ImGui::Text("%s%s", isStatic(b) ? "static" : "dynamic", b.sleeping ? ", sleeping" : "");
        if (!isStatic(b)) ImGui::Text("mass = %.3f", 1.0 / b.invMass);
        if (b.shape == Shape::Sphere) ImGui::Text("radius = %.3f", b.radius);
        else if (b.shape == Shape::Cylinder) ImGui::Text("r = %.3f  halfH = %.3f", b.radius, b.halfHeight);
        else ImGui::Text("half = (%.2f, %.2f, %.2f)", b.halfExtents.x, b.halfExtents.y, b.halfExtents.z);
        ImGui::Separator();
        float pos[3] = {(float)b.x.x, (float)b.x.y, (float)b.x.z};
        if (ImGui::DragFloat3("position", pos, 0.05f)) { b.x = {pos[0], pos[1], pos[2]}; World::wake(b); }
        float vel[3] = {(float)b.v.x, (float)b.v.y, (float)b.v.z};
        if (ImGui::DragFloat3("velocity", vel, 0.05f)) { b.v = {vel[0], vel[1], vel[2]}; World::wake(b); }
        float av[3] = {(float)b.w.x, (float)b.w.y, (float)b.w.z};
        if (ImGui::DragFloat3("ang. vel.", av, 0.05f)) { b.w = {av[0], av[1], av[2]}; World::wake(b); }
        ImGui::Text("quat wxyz  %.3f %.3f %.3f %.3f", b.q.w, b.q.x, b.q.y, b.q.z);
        if (b.sleeping) { ImGui::SameLine(); if (ImGui::SmallButton("wake")) World::wake(b); }
        if (ImGui::Button("give upward impulse")) { b.v.y += 8.0; World::wake(b); }
    } else {
        ImGui::TextDisabled("no body selected");
        ImGui::TextDisabled("click a body in the viewport or the hierarchy");
    }
    ImGui::End();
}

// ------------------------------------------------------------------ 3D overlays
static void buildDebugLines() {
    gLines.clear();
    if (gViz.grid && !gShowCell) {
        float e = 30.0f; Vec3 gc{0.24f, 0.26f, 0.30f};
        for (int i = -30; i <= 30; i += 3) { line({(float)i, 0, -e}, {(float)i, 0, e}, gc); line({-e, 0, (float)i}, {e, 0, (float)i}, gc); }
    }
    if (gViz.arrows) {
        for (const Body& b : gWorld.bodies) {
            if (isStatic(b)) continue;
            Vec3 p = v3(b.x), tip = p + v3(b.v) * gViz.velScale;
            float sp = (float)b.v.norm();
            Vec3 col = sp > 4 ? Vec3{1.0f, 0.4f, 0.2f} : Vec3{0.4f, 0.9f, 0.5f};
            line(p, tip, col);
        }
    }
    if (gViz.aabb) {
        for (const Body& b : gWorld.bodies)
            if (!isStatic(b)) aabb(v3(b.x), (float)b.boundingRadius(), {0.3f, 0.5f, 0.35f});
    }
    if (gViz.contacts) {
        for (const ContactViz& c : gWorld.debugContacts) {
            Vec3 p = v3(c.point);
            cross3(p, 0.18f, {1.0f, 0.95f, 0.2f});
            line(p, p + v3(c.normal) * 0.7f, {0.2f, 0.85f, 1.0f});
        }
    }
    if (gSelected >= 0 && gSelected < (int)gWorld.bodies.size())
        aabb(v3(gWorld.bodies[gSelected].x), (float)gWorld.bodies[gSelected].boundingRadius() * 1.05f, {1, 1, 1});

    // Level geometry drawn as line work: triangle meshes (wireframe) + convex AABBs.
    for (const Body& b : gWorld.bodies) {
        if (b.shape == Shape::Mesh && b.mesh) {
            const MeshData& md = *b.mesh; Vec3 mc{0.30f, 0.45f, 0.42f};
            for (size_t t = 0; t < md.tris.size(); ++t) {
                Vec3 w[3];
                for (int k = 0; k < 3; ++k) w[k] = v3(b.x + b.q.rotate(md.v((int)t, k)));
                line(w[0], w[1], mc); line(w[1], w[2], mc); line(w[2], w[0], mc);
            }
        } else if (b.shape == Shape::Convex) {
            aabb(v3(b.x), (float)b.boundingRadius(), {0.5f, 0.6f, 0.4f});
        }
    }
    // Joints: a line between the two anchor points.
    if (gViz.joints) {
        Vec3 jc{0.9f, 0.6f, 0.2f};
        for (const DistanceJoint& j : gWorld.distanceJoints)
            line(v3(gWorld.bodies[j.a].x + gWorld.bodies[j.a].q.rotate(j.localA)),
                 v3(gWorld.bodies[j.b].x + gWorld.bodies[j.b].q.rotate(j.localB)), jc);
        for (const Joint& j : gWorld.joints) {
            if (j.broken) continue;
            line(v3(gWorld.bodies[j.a].x + gWorld.bodies[j.a].q.rotate(j.localA)),
                 v3(gWorld.bodies[j.b].x + gWorld.bodies[j.b].q.rotate(j.localB)), jc);
        }
    }
}

static void sceneTick() {
    if (gKinBody < gWorld.bodies.size())   // scene 12: oscillate the platform
        gWorld.bodies[gKinBody].v = {4.0 * std::cos(gSimTime * 0.8), 0, 0};
}

// WASD capsule character: collide-and-slide against the world (see the engine's
// CharacterController; done inline here so it can share the viewer's gWorld).
static void updateChar(double dt) {
    if (!gChar.active) return;
    V3 fwd{-std::sin(gAz), 0, -std::cos(gAz)};
    V3 right{std::cos(gAz), 0, -std::sin(gAz)};
    V3 wish;
    if (gKey['w']) wish += fwd;   if (gKey['s']) wish -= fwd;
    if (gKey['d']) wish += right; if (gKey['a']) wish -= right;
    double wl = wish.norm(); if (wl > 1e-6) wish = wish * (5.0 / wl);
    gChar.vel.x = wish.x; gChar.vel.z = wish.z;
    gChar.vel.y += gWorld.gravity.y * dt;
    if (gChar.grounded && gKey['e']) { gChar.vel.y = 7.0; gChar.grounded = false; }
    gChar.pos += gChar.vel * dt;
    Body cap = makeCapsule(-1, gChar.pos, Q{1, 0, 0, 0}, gChar.r, gChar.hh, 1.0);
    gChar.grounded = false;
    for (int it = 0; it < 6; ++it) {
        double deepest = 0.0; V3 dn;
        for (const Body& b : gWorld.bodies) {
            if (b.sensor) continue;
            Contact ct = detectContact(cap, b, gWorld.box);
            if (ct.hit && ct.overlap > deepest) { deepest = ct.overlap; dn = ct.normal; }
        }
        if (deepest <= 0.01) break;
        gChar.pos += dn * (deepest - 0.01);
        if (dn.y > 0.6) { gChar.grounded = true; if (gChar.vel.y < 0) gChar.vel.y = 0; }
        cap.x = gChar.pos;
    }
    gCenter = {(float)gChar.pos.x, (float)gChar.pos.y + 1.0f, (float)gChar.pos.z};
}

static void display() {
    // step the engine live (unless paused)
    gWorld.captureContacts = gViz.contacts;
    if (!gPaused) {
        sceneTick();
        for (int i = 0; i < gStepsPerFrame; ++i) { gWorld.step(); gSimTime += gWorld.dt; }
        updateChar(gStepsPerFrame * gWorld.dt);
    }
    pushKE((float)gWorld.totalKinetic());

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGLUT_NewFrame();
    ImGui::NewFrame();
    buildUI();
    ImGui::Render();

    glEnable(GL_DEPTH_TEST); glDisable(GL_BLEND); glDisable(GL_SCISSOR_TEST);
    glClearColor(0.09f, 0.10f, 0.12f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    float aspect = gWinW / (float)(gWinH ? gWinH : 1);
    Vec3 eye{gCenter.x + gDist * std::cos(gEl) * std::sin(gAz),
             gCenter.y + gDist * std::sin(gEl),
             gCenter.z + gDist * std::cos(gEl) * std::cos(gAz)};
    Mat4 VP = Mat4::perspective(45.0f * 3.14159f / 180.0f, aspect, 0.1f, 5000.0f) *
              Mat4::lookAt(eye, gCenter, {0, 1, 0});

    glUseProgram(gProg);
    glUniformMatrix4fv(gUVP, 1, GL_FALSE, VP.m);
    glUniform3fv(gUCam, 1, &eye.x);
    Vec3 lightDir = Vec3{-0.4f, -1.0f, -0.3f}.norm();
    glUniform3fv(gULight, 1, &lightDir.x);

    if (gShowCell || gViz.showCellAlways) {
        float R = (float)gWorld.box.half;
        gWire.instData.clear(); gWire.instCount = 0;
        pushInstance(gWire, Mat4::scale({R, R, R}), {0.40f, 0.45f, 0.55f});
        drawMesh(gWire, 1);
    }

    gSphere.instData.clear(); gSphere.instCount = 0;
    gCyl.instData.clear(); gCyl.instCount = 0;
    gBox.instData.clear(); gBox.instCount = 0;

    // Push one primitive at a world transform. Capsule = cylinder + two cap spheres.
    auto pushShape = [&](Shape s, const Vec3& p, const Mat4& Rm, float r, float hh,
                         const Vec3& he, const Vec3& col) {
        Mat4 T = Mat4::translate(p);
        if (s == Shape::Sphere) pushInstance(gSphere, T * Rm * Mat4::scale({r, r, r}), col);
        else if (s == Shape::Box) pushInstance(gBox, T * Rm * Mat4::scale(he), col);
        else if (s == Shape::Cylinder) pushInstance(gCyl, T * Rm * Mat4::scale({r, hh, r}), col);
        else if (s == Shape::Capsule) {
            pushInstance(gCyl, T * Rm * Mat4::scale({r, hh, r}), col);
            Vec3 up{Rm.m[4], Rm.m[5], Rm.m[6]};   // body +Y in world (Rm column 1)
            pushInstance(gSphere, Mat4::translate(p + up * hh) * Mat4::scale({r, r, r}), col);
            pushInstance(gSphere, Mat4::translate(p - up * hh) * Mat4::scale({r, r, r}), col);
        }
        // Convex / Plane / Mesh are drawn as line geometry elsewhere.
    };
    auto emit = [&](const Body& b, const Vec3& shift, const Vec3& col) {
        Vec3 p{(float)b.x.x + shift.x, (float)b.x.y + shift.y, (float)b.x.z + shift.z};
        Mat4 Rm = Mat4::fromQuat((float)b.q.w, (float)b.q.x, (float)b.q.y, (float)b.q.z);
        if (b.shape == Shape::Compound) {
            for (const ChildShape& ch : b.children) {
                Vec3 lp{(float)ch.localPos.x, (float)ch.localPos.y, (float)ch.localPos.z};
                Vec3 wp = p + Vec3{Rm.m[0]*lp.x + Rm.m[4]*lp.y + Rm.m[8]*lp.z,
                                   Rm.m[1]*lp.x + Rm.m[5]*lp.y + Rm.m[9]*lp.z,
                                   Rm.m[2]*lp.x + Rm.m[6]*lp.y + Rm.m[10]*lp.z};
                Mat4 cRm = Rm * Mat4::fromQuat((float)ch.localRot.w, (float)ch.localRot.x, (float)ch.localRot.y, (float)ch.localRot.z);
                pushShape(ch.shape, wp, cRm, (float)ch.radius, (float)ch.halfHeight,
                          {(float)ch.halfExtents.x, (float)ch.halfExtents.y, (float)ch.halfExtents.z}, col);
            }
        } else {
            pushShape(b.shape, p, Rm, (float)b.radius, (float)b.halfHeight,
                      {(float)b.halfExtents.x, (float)b.halfExtents.y, (float)b.halfExtents.z}, col);
        }
    };
    float edge = (float)gWorld.box.edge();
    for (size_t i = 0; i < gWorld.bodies.size(); ++i) {
        const Body& b = gWorld.bodies[i];
        Vec3 col = palette(b.id, b.shape, b.sleeping, isStatic(b));
        if ((int)i == gSelected) col = col * 0.5f + Vec3{0.5f, 0.5f, 0.2f};
        emit(b, {0, 0, 0}, col);
        if (gViz.ghosts && gWorld.box.periodic && !isStatic(b)) {
            float h = (float)gWorld.box.half, r = (float)b.boundingRadius();
            float cc[3] = {(float)b.x.x, (float)b.x.y, (float)b.x.z};
            for (int a = 0; a < 3; ++a) {
                if (cc[a] > h - r) { Vec3 s{0, 0, 0}; (&s.x)[a] = -edge; emit(b, s, col * 0.5f); }
                if (cc[a] < -h + r) { Vec3 s{0, 0, 0}; (&s.x)[a] = edge; emit(b, s, col * 0.5f); }
            }
        }
    }
    if (gChar.active) {   // the WASD character capsule (not a world body)
        Vec3 cp{(float)gChar.pos.x, (float)gChar.pos.y, (float)gChar.pos.z};
        pushShape(Shape::Capsule, cp, Mat4::identity(), (float)gChar.r, (float)gChar.hh, {},
                  gChar.grounded ? Vec3{0.95f, 0.85f, 0.30f} : Vec3{0.95f, 0.55f, 0.30f});
    }
    drawMesh(gSphere, 0); drawMesh(gCyl, 0); drawMesh(gBox, 0);

    buildDebugLines();
    if (!gLines.empty()) {
        glUseProgram(gLineProg);
        glUniformMatrix4fv(gLUVP, 1, GL_FALSE, VP.m);
        glBindVertexArray(gLineVao);
        glBindBuffer(GL_ARRAY_BUFFER, gLineVbo);
        glBufferData(GL_ARRAY_BUFFER, gLines.size() * sizeof(float), gLines.data(), GL_DYNAMIC_DRAW);
        glDrawArrays(GL_LINES, 0, (GLsizei)(gLines.size() / 6));
        glBindVertexArray(0);
    }

    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glutSwapBuffers();
}

// ------------------------------------------------------------------ input
static void reshape(int w, int h) { gWinW = w; gWinH = h; glViewport(0, 0, w, h); ImGui_ImplGLUT_ReshapeFunc(w, h); }
static void timer(int) { glutPostRedisplay(); glutTimerFunc(16, timer, 0); }

static void key(unsigned char k, int x, int y) {
    ImGui_ImplGLUT_KeyboardFunc(k, x, y);
    if (ImGui::GetIO().WantCaptureKeyboard) return;
    gKey[k] = true;                          // for the WASD character
    switch (k) {
        case 27: glutLeaveMainLoop(); break;
        case ' ': gPaused = !gPaused; break;
        case 'r': buildScene(gScene); break;
        case 'n': buildScene(gScene % kNumScenes + 1); break;
        case 'd': dropSphere(); break;
        case ',': if (gStepsPerFrame > 1) --gStepsPerFrame; break;
        case '.': if (gStepsPerFrame < 16) ++gStepsPerFrame; break;
        case 'h': gDist = gWorld.box.half * (gShowCell ? 2.8f : 1.4f);
                  gAz = 0.7f; gEl = 0.45f; gCenter = {0, gShowCell ? 0.0f : 5.0f, 0}; break;
        default: if (k >= '1' && k <= '9') buildScene(k - '0');
    }
}
static void keyUp(unsigned char k, int x, int y) { gKey[k] = false; ImGui_ImplGLUT_KeyboardUpFunc(k, x, y); }
static void special(int k, int x, int y) { ImGui_ImplGLUT_SpecialFunc(k, x, y); }
static void specialUp(int k, int x, int y) { ImGui_ImplGLUT_SpecialUpFunc(k, x, y); }

static void mouse(int b, int st, int x, int y) {
    ImGui_ImplGLUT_MouseFunc(b, st, x, y);
    if (ImGui::GetIO().WantCaptureMouse) { gDrag = false; return; }
    if (b == GLUT_LEFT_BUTTON) {
        gDrag = (st == GLUT_DOWN); gLastX = x; gLastY = y;
        if (st == GLUT_DOWN) pick(x, y);
    }
}
static void wheel(int b, int dir, int x, int y) {
    ImGui_ImplGLUT_MouseWheelFunc(b, dir, x, y);
    if (ImGui::GetIO().WantCaptureMouse) return;
    gDist *= (dir > 0) ? 0.9f : 1.1f;
}
static void motion(int x, int y) {
    ImGui_ImplGLUT_MotionFunc(x, y);
    if (ImGui::GetIO().WantCaptureMouse) return;
    if (gDrag) { gAz -= (x - gLastX) * 0.01f; gEl += (y - gLastY) * 0.01f;
        if (gEl > 1.5f) gEl = 1.5f; if (gEl < -1.5f) gEl = -1.5f; gLastX = x; gLastY = y; }
}
static void passiveMotion(int x, int y) { ImGui_ImplGLUT_MotionFunc(x, y); }

int main(int argc, char** argv) {
    int scene = 1;
    if (argc >= 2) scene = std::atoi(argv[1]);
    if (scene < 1 || scene > kNumScenes) scene = 1;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutSetOption(GLUT_MULTISAMPLE, 8);
    glutInitContextVersion(3, 3); glutInitContextProfile(GLUT_CORE_PROFILE);
    glutInitWindowSize(gWinW, gWinH);
    glutCreateWindow("NativeEngine Visual Debugger");
    if (!loadGL()) { std::fprintf(stderr, "failed to load OpenGL 3.3\n"); return 1; }
    glEnable(GL_DEPTH_TEST); glEnable(GL_MULTISAMPLE);

    gProg = link(kVert, kFrag);
    gUVP = glGetUniformLocation(gProg, "uVP"); gUCam = glGetUniformLocation(gProg, "uCam");
    gULight = glGetUniformLocation(gProg, "uLight"); gUUnlit = glGetUniformLocation(gProg, "uUnlit");
    gLineProg = link(kLineVert, kLineFrag); gLUVP = glGetUniformLocation(gLineProg, "uVP");
    glGenVertexArrays(1, &gLineVao); glBindVertexArray(gLineVao);
    glGenBuffers(1, &gLineVbo); glBindBuffer(GL_ARRAY_BUFFER, gLineVbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0); glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float))); glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    makeSphere(gSphere); makeCylinder(gCyl); makeBox(gBox); makeWireBox(gWire);
    buildScene(scene);

    IMGUI_CHECKVERSION(); ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui::GetIO().FontGlobalScale = 1.15f;
    ImGui_ImplGLUT_Init();
    ImGui_ImplOpenGL3_Init("#version 330");

    glutDisplayFunc(display); glutReshapeFunc(reshape);
    glutKeyboardFunc(key); glutKeyboardUpFunc(keyUp);
    glutSpecialFunc(special); glutSpecialUpFunc(specialUp);
    glutMouseFunc(mouse); glutMotionFunc(motion); glutPassiveMotionFunc(passiveMotion);
    glutMouseWheelFunc(wheel);
    glutTimerFunc(16, timer, 0);

    std::printf("NativeEngine Visual Debugger -- live, interactive. Scenes 1..15 (10-15: ragdoll, "
                "terrain, kinematic platform, hinge bridge, compound, WASD character). Click to select.\n");
    glutMainLoop();

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGLUT_Shutdown();
    ImGui::DestroyContext();
    return 0;
}
