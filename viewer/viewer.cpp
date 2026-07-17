// NativeEngine LIVE demo viewer. Unlike a recording player, this steps the
// physics engine in real time each frame and renders the current state -- a real
// interactive demo (PhysX-snippet style): switch scenes, pause, reset, drop
// bodies, orbit. Modern OpenGL 3.3 core (instanced Blinn-Phong, MSAA), right-
// handed / Y-up so the view matches NativeEngine.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "gl_core.hpp"
#include "glmath.hpp"

#include "world.hpp"
#include "body.hpp"

using namespace nv;
using namespace ne;

// ------------------------------------------------------------------ sim state
static World gWorld;
static int   gScene = 1;
static int   gStepsPerFrame = 4;    // dt = 1/240 * 4 ~ 1/60 s per frame -> real time
static bool  gPaused = false;
static std::mt19937_64 gRng(12345);

// camera / window
static float gAz = 0.7f, gEl = 0.45f, gDist = 40.0f;
static Vec3  gCenter{0, 0, 0};
static int   gLastX = 0, gLastY = 0; static bool gDrag = false;
static int   gWinW = 1280, gWinH = 800;

static GLuint gProg = 0;
static GLint  gUVP = -1, gUCam = -1, gULight = -1, gUUnlit = -1;

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
    // 6 faces, unit half-extent [-1,1], outward normals.
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
static Vec3 palette(int id, Shape s, bool sleeping, bool isStatic) {
    if (isStatic) return {0.30f, 0.31f, 0.34f};
    float t = (float)((id * 2654435761u) % 997) / 997.0f;
    Vec3 c;
    if (s == Shape::Sphere) c = {0.90f, 0.48f + 0.22f * t, 0.26f};
    else if (s == Shape::Box) c = {0.35f + 0.2f * t, 0.62f, 0.42f + 0.2f * t};
    else c = {0.28f + 0.15f * t, 0.55f + 0.2f * t, 0.80f};
    if (sleeping) c = c * 0.55f;   // dim sleeping bodies
    return c;
}

// ------------------------------------------------------------------ scenes
// Modelled on the classic PhysX snippets: objects on an open GROUND plane, not
// inside a reflective box. The box is used only for the "elastic gas" (a
// genuine container) and the periodic cell for the PBC feature demos.
static int gNextId = 0;
static bool gShowCell = false;     // draw the domain wireframe (boxed/periodic only)
static const int kNumScenes = 9;

static double urand(double a, double b) { std::uniform_real_distribution<double> d(a, b); return d(gRng); }

// A large static ground slab (open scenes rest on this and can roll off edges).
static void addGround(double half = 30.0) {
    Body f = makeBox(gNextId++, {0, -1.0, 0}, Q{1, 0, 0, 0}, {half, 1.0, half}, 1.0);
    f.invMass = 0; f.invInertiaBody = {0, 0, 0}; f.dynamic = false;
    gWorld.bodies.push_back(f);
}
static void openGravityWorld() {
    gWorld.openBoundary = true; gWorld.gravity = {0, -10, 0};
    gWorld.friction = 0.6; gWorld.restitution = 0.0; gWorld.sleepEnabled = true;
    gWorld.box.half = 40;   // only affects the (hidden) domain + camera fit
}

static void buildScene(int id) {
    gWorld = World{};
    gWorld.dt = 1.0 / 240.0;
    gNextId = 0; gScene = id; gShowCell = false;
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
            openGravityWorld();
            gWorld.sleepEnabled = false;       // keep the wall live for the hit
            addGround();
            const int cols = 8, rows = 6;
            for (int r = 0; r < rows; ++r)
                for (int c = 0; c < cols; ++c) {
                    double off = (r % 2) ? 1.0 : 0.0;   // running bond
                    gWorld.bodies.push_back(makeBox(gNextId++,
                        {(c - cols / 2.0) * 2.0 + off, 1.0 + r * 2.0, 0},
                        Q{1, 0, 0, 0}, {1.0, 1.0, 2.0}, 1.0));
                }
            Body ball = makeSphere(gNextId++, {0, 6, 22}, 1.6, 6.0);   // heavy, fast
            ball.v = {0, 0, -35};
            gWorld.bodies.push_back(ball);
            break;
        }
        case 3: {  // dominoes: tip the first, watch the chain topple
            openGravityWorld(); addGround();
            const int N = 16;
            for (int i = 0; i < N; ++i)
                gWorld.bodies.push_back(makeBox(gNextId++,
                    {(i - N / 2.0) * 1.5, 2.0, 0}, Q{1, 0, 0, 0}, {0.25, 2.0, 1.2}, 1.0));
            gWorld.bodies[1].w = {0, 0, -3.0};   // tip the first domino
            gWorld.bodies[1].v = {2.0, 0, 0};
            break;
        }
        case 4: {  // spheres poured into a heap on the ground
            openGravityWorld(); gWorld.restitution = 0.1; gWorld.friction = 0.4;
            addGround();
            for (int i = 0; i < 140; ++i)
                gWorld.bodies.push_back(makeSphere(gNextId++, {urand(-6, 6), urand(3, 20), urand(-6, 6)}, 1.0, 1.0));
            break;
        }
        case 5: {  // cylinders tumbling into a heap on the ground
            openGravityWorld(); gWorld.restitution = 0.05;
            addGround();
            for (int i = 0; i < 60; ++i) {
                Q q{urand(-1, 1), urand(-1, 1), urand(-1, 1), urand(-1, 1)}; q = q.normalized();
                gWorld.bodies.push_back(makeCylinder(gNextId++, {urand(-6, 6), urand(3, 20), urand(-6, 6)}, q, 1.0, 4.0, 1.0));
            }
            break;
        }
        case 6: {  // hanging rope/chain of spheres (rigid joints), released horizontal
            gWorld.openBoundary = true; gWorld.gravity = {0, -10, 0}; gWorld.box.half = 24;
            Body anchor = makeSphere(gNextId++, {0, 14, 0}, 0.3, 1.0);
            anchor.invMass = 0; anchor.invInertiaBody = {0, 0, 0}; anchor.dynamic = false;
            gWorld.bodies.push_back(anchor);
            const int N = 12; const double L = 1.4;
            for (int i = 1; i <= N; ++i) gWorld.bodies.push_back(makeSphere(gNextId++, {i * L, 14, 0}, 0.5, 1.0));
            for (int i = 0; i < N; ++i) { DistanceJoint j; j.a = i; j.b = i + 1; j.rest = L; gWorld.distanceJoints.push_back(j); }
            break;
        }
        case 7: {  // elastic sphere gas in a reflective BOX (a genuine container)
            gWorld.box.half = 10; gWorld.restitution = 1.0; gShowCell = true;
            for (int i = 0; i < 130; ++i) {
                Body b = makeSphere(gNextId++, {urand(-8, 8), urand(-8, 8), urand(-8, 8)}, 0.8, 1.0);
                b.v = {urand(-5, 5), urand(-5, 5), urand(-5, 5)}; gWorld.bodies.push_back(b);
            }
            break;
        }
        case 8: {  // PERIODIC gas -- the distinctive feature (bodies wrap faces)
            gWorld.box.half = 9; gWorld.box.periodic = true; gWorld.restitution = 1.0; gShowCell = true;
            for (int i = 0; i < 160; ++i) {
                Body b = makeSphere(gNextId++, {urand(-8.5, 8.5), urand(-8.5, 8.5), urand(-8.5, 8.5)}, 0.8, 1.0);
                b.v = {urand(-4, 4), urand(-4, 4), urand(-4, 4)}; gWorld.bodies.push_back(b);
            }
            break;
        }
        case 9: {  // two spheres colliding ACROSS the periodic boundary
            gWorld.box.half = 6; gWorld.box.periodic = true; gWorld.restitution = 1.0; gShowCell = true;
            Body a = makeSphere(gNextId++, {5.2, 0, 0}, 1.0, 1.0); a.v = {2, 0, 0};
            Body b = makeSphere(gNextId++, {-5.2, 0.3, 0}, 1.0, 1.0); b.v = {-2, 0, 0};
            gWorld.bodies.push_back(a); gWorld.bodies.push_back(b);
            break;
        }
    }
    gDist = gWorld.box.half * (gShowCell ? 2.8f : 1.4f);
    gCenter = {0, gShowCell ? 0.0f : 5.0f, 0};
}
static const char* sceneName(int id) {
    switch (id) {
        case 1: return "box pyramid"; case 2: return "brick wall + projectile"; case 3: return "dominoes";
        case 4: return "sphere pile"; case 5: return "cylinder pile"; case 6: return "rope";
        case 7: return "elastic gas (box)"; case 8: return "PERIODIC gas"; case 9: return "PERIODIC pair"; }
    return "?";
}
static void dropSphere() {   // 'd' -- drop a body onto the scene, live
    Body b = makeSphere(gNextId++, {urand(-4, 4), 20.0, urand(-4, 4)}, 1.2, 1.0);
    gWorld.bodies.push_back(b);
}

// ------------------------------------------------------------------ render
static void display() {
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

    if (gShowCell) {   // the reflective box / periodic cell -- not for open scenes
        float R = (float)gWorld.box.half;
        gWire.instData.clear(); gWire.instCount = 0;
        pushInstance(gWire, Mat4::scale({R, R, R}), {0.40f, 0.45f, 0.55f});
        drawMesh(gWire, 1);
    }

    gSphere.instData.clear(); gSphere.instCount = 0;
    gCyl.instData.clear(); gCyl.instCount = 0;
    gBox.instData.clear(); gBox.instCount = 0;
    for (const Body& b : gWorld.bodies) {
        Mat4 T = Mat4::translate({(float)b.x.x, (float)b.x.y, (float)b.x.z});
        Mat4 Rm = Mat4::fromQuat((float)b.q.w, (float)b.q.x, (float)b.q.y, (float)b.q.z);
        bool isStatic = (!b.dynamic && b.invMass == 0.0);
        Vec3 col = palette(b.id, b.shape, b.sleeping, isStatic);
        if (b.shape == Shape::Sphere)
            pushInstance(gSphere, T * Rm * Mat4::scale({(float)b.radius, (float)b.radius, (float)b.radius}), col);
        else if (b.shape == Shape::Box)
            pushInstance(gBox, T * Rm * Mat4::scale({(float)b.halfExtents.x, (float)b.halfExtents.y, (float)b.halfExtents.z}), col);
        else
            pushInstance(gCyl, T * Rm * Mat4::scale({(float)b.radius, (float)b.halfHeight, (float)b.radius}), col);
    }
    drawMesh(gSphere, 0); drawMesh(gCyl, 0); drawMesh(gBox, 0);
    glutSwapBuffers();
}

static void reshape(int w, int h) { gWinW = w; gWinH = h; glViewport(0, 0, w, h); }
static void timer(int) {
    if (!gPaused) for (int i = 0; i < gStepsPerFrame; ++i) gWorld.step();
    glutPostRedisplay();
    glutTimerFunc(16, timer, 0);
}
static void key(unsigned char k, int, int) {
    switch (k) {
        case 27: glutLeaveMainLoop(); break;
        case ' ': gPaused = !gPaused; break;
        case 'r': buildScene(gScene); break;                       // reset current scene
        case 'n': buildScene(gScene % kNumScenes + 1); break;      // next scene
        case 'd': dropSphere(); break;                             // drop a body (live)
        case ',': if (gStepsPerFrame > 1) --gStepsPerFrame; break; // slower
        case '.': ++gStepsPerFrame; break;                         // faster
        case 'h': gDist = gWorld.box.half * (gShowCell ? 2.8f : 1.4f);
                  gAz = 0.7f; gEl = 0.45f; gCenter = {0, gShowCell ? 0.0f : 5.0f, 0}; break;
        default:
            if (k >= '1' && k <= '9') buildScene(k - '0');
    }
    std::printf("scene %d (%s)  bodies=%zu  %s  stepsPerFrame=%d\n",
                gScene, sceneName(gScene), gWorld.bodies.size(), gPaused ? "PAUSED" : "running", gStepsPerFrame);
    glutPostRedisplay();
}
static void mouse(int b, int st, int x, int y) {
    if (b == GLUT_LEFT_BUTTON) { gDrag = (st == GLUT_DOWN); gLastX = x; gLastY = y; }
    if (b == 3) gDist *= 0.9f; if (b == 4) gDist *= 1.1f;
    glutPostRedisplay();
}
static void motion(int x, int y) {
    if (gDrag) { gAz -= (x - gLastX) * 0.01f; gEl += (y - gLastY) * 0.01f;
        if (gEl > 1.5f) gEl = 1.5f; if (gEl < -1.5f) gEl = -1.5f; gLastX = x; gLastY = y; glutPostRedisplay(); }
}

int main(int argc, char** argv) {
    int scene = 1;
    if (argc >= 2) scene = std::atoi(argv[1]);
    if (scene < 1 || scene > 7) scene = 1;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutSetOption(GLUT_MULTISAMPLE, 8);
    glutInitContextVersion(3, 3); glutInitContextProfile(GLUT_CORE_PROFILE);
    glutInitWindowSize(gWinW, gWinH);
    glutCreateWindow("NativeEngine LIVE  [1-7 scene  n next  r reset  d drop  space pause  , . speed  h cam  drag orbit]");
    if (!loadGL()) { std::fprintf(stderr, "failed to load OpenGL 3.3\n"); return 1; }
    glEnable(GL_DEPTH_TEST); glEnable(GL_MULTISAMPLE);
    gProg = link(kVert, kFrag);
    gUVP = glGetUniformLocation(gProg, "uVP"); gUCam = glGetUniformLocation(gProg, "uCam");
    gULight = glGetUniformLocation(gProg, "uLight"); gUUnlit = glGetUniformLocation(gProg, "uUnlit");
    makeSphere(gSphere); makeCylinder(gCyl); makeBox(gBox); makeWireBox(gWire);
    buildScene(scene);

    std::printf("NativeEngine live viewer. Scenes: 1 box-stack  2 pile  3 cylinders  4 gas  5 PERIODIC-gas  6 pendulum  7 PERIODIC-pair\n");
    glutDisplayFunc(display); glutReshapeFunc(reshape);
    glutKeyboardFunc(key); glutMouseFunc(mouse); glutMotionFunc(motion);
    glutTimerFunc(16, timer, 0);
    glutMainLoop();
    return 0;
}
