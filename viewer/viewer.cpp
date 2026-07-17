// NativeEngine viewer -- a modern (OpenGL 3.3 core) instanced renderer for .pxrf
// recordings. Blinn-Phong + hemisphere ambient, MSAA, an orbit camera, a ground
// plane and a domain wireframe (the periodic cell). Right-handed / Y-up so the
// view matches NativeEngine exactly (no mirroring).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "gl_core.hpp"
#include "glmath.hpp"
#include "recording.hpp"

using namespace nv;

// ------------------------------------------------------------------ globals
static Recording gRec;
static int   gFrame = 0;
static bool  gPaused = false;
static float gSpeed = 1.0f;         // frames advanced per tick
static float gAccum = 0.0f;
static int   gWinW = 1280, gWinH = 800;

// camera (orbit)
static float gAz = 0.7f, gEl = 0.5f, gDist = 40.0f;
static Vec3  gCenter{0, 0, 0};
static int   gLastX = 0, gLastY = 0; static bool gDrag = false;

static GLuint gProg = 0;
static GLint  gUVP = -1, gUCam = -1, gULight = -1, gUUnlit = -1;

struct Mesh {
    GLuint vao = 0, vbo = 0, ebo = 0, inst = 0;
    GLsizei indexCount = 0;
    GLenum  prim = GL_TRIANGLES;
    std::vector<float> instData;    // 16 (model) + 3 (color) per instance
    int instCount = 0;
};
static Mesh gSphere, gCyl, gGround, gWire;

// ------------------------------------------------------------------ shaders
static const char* kVert =
    "#version 330 core\n"
    "layout(location=0) in vec3 aPos;\n"
    "layout(location=1) in vec3 aNormal;\n"
    "layout(location=3) in vec4 m0;\n"
    "layout(location=4) in vec4 m1;\n"
    "layout(location=5) in vec4 m2;\n"
    "layout(location=6) in vec4 m3;\n"
    "layout(location=7) in vec3 aColor;\n"
    "uniform mat4 uVP;\n"
    "out vec3 vN; out vec3 vW; out vec3 vC;\n"
    "void main(){ mat4 M=mat4(m0,m1,m2,m3); vec4 w=M*vec4(aPos,1.0);\n"
    "  vW=w.xyz; vN=mat3(M)*aNormal; vC=aColor; gl_Position=uVP*w; }\n";

static const char* kFrag =
    "#version 330 core\n"
    "in vec3 vN; in vec3 vW; in vec3 vC; out vec4 o;\n"
    "uniform vec3 uCam; uniform vec3 uLight; uniform int uUnlit;\n"
    "void main(){ if(uUnlit==1){ o=vec4(vC,1.0); return; }\n"
    "  vec3 N=normalize(vN); vec3 L=normalize(-uLight); vec3 V=normalize(uCam-vW);\n"
    "  vec3 H=normalize(L+V);\n"
    "  float diff=max(dot(N,L),0.0); float spec=pow(max(dot(N,H),0.0),48.0);\n"
    "  float hemi=0.5+0.5*N.y;\n"
    "  vec3 amb=mix(vec3(0.16,0.17,0.20),vec3(0.34,0.36,0.40),hemi);\n"
    "  vec3 c=vC*(amb+diff*0.9)+vec3(1.0)*spec*0.35;\n"
    "  c=c/(c+vec3(1.0)); c=pow(c,vec3(1.0/2.2)); o=vec4(c,1.0); }\n";

static GLuint compile(GLenum t, const char* src) {
    GLuint s = glCreateShader(t);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) { char log[1024]; glGetShaderInfoLog(s, 1024, nullptr, log); std::fprintf(stderr, "shader: %s\n", log); }
    return s;
}
static GLuint link(const char* v, const char* f) {
    GLuint p = glCreateProgram();
    GLuint vs = compile(GL_VERTEX_SHADER, v), fs = compile(GL_FRAGMENT_SHADER, f);
    glAttachShader(p, vs); glAttachShader(p, fs); glLinkProgram(p);
    GLint ok = 0; glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) { char log[1024]; glGetProgramInfoLog(p, 1024, nullptr, log); std::fprintf(stderr, "link: %s\n", log); }
    glDeleteShader(vs); glDeleteShader(fs);
    return p;
}

// ------------------------------------------------------------------ meshes
static void upload(Mesh& m, const std::vector<float>& verts, const std::vector<unsigned>& idx, GLenum prim) {
    m.prim = prim; m.indexCount = (GLsizei)idx.size();
    glGenVertexArrays(1, &m.vao); glBindVertexArray(m.vao);
    glGenBuffers(1, &m.vbo); glBindBuffer(GL_ARRAY_BUFFER, m.vbo);
    glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_STATIC_DRAW);
    glGenBuffers(1, &m.ebo); glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(unsigned), idx.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // instance buffer: mat4 (loc 3..6) + color (loc 7)
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

static void makeSphere(Mesh& m, int stacks = 24, int slices = 32) {
    std::vector<float> v; std::vector<unsigned> idx;
    const float PI = 3.14159265f;
    for (int i = 0; i <= stacks; ++i) {
        float phi = i / (float)stacks * PI;
        for (int j = 0; j <= slices; ++j) {
            float th = j / (float)slices * 2 * PI;
            float x = std::sin(phi) * std::cos(th), y = std::cos(phi), z = std::sin(phi) * std::sin(th);
            v.insert(v.end(), {x, y, z, x, y, z});
        }
    }
    int row = slices + 1;
    for (int i = 0; i < stacks; ++i)
        for (int j = 0; j < slices; ++j) {
            unsigned a = i * row + j, b = a + row;
            idx.insert(idx.end(), {a, b, a + 1u, a + 1u, b, b + 1u});
        }
    upload(m, v, idx, GL_TRIANGLES);
}

static void makeCylinder(Mesh& m, int slices = 32) {
    std::vector<float> v; std::vector<unsigned> idx;
    const float PI = 3.14159265f;
    // side
    for (int j = 0; j <= slices; ++j) {
        float th = j / (float)slices * 2 * PI, x = std::cos(th), z = std::sin(th);
        v.insert(v.end(), {x, -1, z, x, 0, z});
        v.insert(v.end(), {x,  1, z, x, 0, z});
    }
    for (int j = 0; j < slices; ++j) {
        unsigned a = j * 2u; idx.insert(idx.end(), {a, a + 1u, a + 2u, a + 2u, a + 1u, a + 3u});
    }
    // caps
    for (int sgn = -1; sgn <= 1; sgn += 2) {
        unsigned base = (unsigned)(v.size() / 6);
        float y = (float)sgn;
        v.insert(v.end(), {0, y, 0, 0, y, 0});           // center
        for (int j = 0; j <= slices; ++j) {
            float th = j / (float)slices * 2 * PI;
            v.insert(v.end(), {std::cos(th), y, std::sin(th), 0, y, 0});
        }
        for (int j = 0; j < slices; ++j) {
            if (sgn < 0) idx.insert(idx.end(), {base, base + 1u + j, base + 2u + j});
            else         idx.insert(idx.end(), {base, base + 2u + j, base + 1u + j});
        }
    }
    upload(m, v, idx, GL_TRIANGLES);
}

static void makeGround(Mesh& m) {
    std::vector<float> v = {
        -1, 0, -1, 0, 1, 0,  1, 0, -1, 0, 1, 0,  1, 0, 1, 0, 1, 0,  -1, 0, 1, 0, 1, 0};
    std::vector<unsigned> idx = {0, 2, 1, 0, 3, 2};
    upload(m, v, idx, GL_TRIANGLES);
}

static void makeWireBox(Mesh& m) {
    std::vector<float> v;
    for (int i = 0; i < 8; ++i) {
        float x = (i & 1) ? 1 : -1, y = (i & 2) ? 1 : -1, z = (i & 4) ? 1 : -1;
        v.insert(v.end(), {x, y, z, 0, 1, 0});
    }
    std::vector<unsigned> idx = {0,1, 1,3, 3,2, 2,0, 4,5, 5,7, 7,6, 6,4, 0,4, 1,5, 2,6, 3,7};
    upload(m, v, idx, GL_LINES);
}

// ------------------------------------------------------------------ instances
static Vec3 palette(int id, bool sphere) {
    // pleasant warm (spheres) / cool (cylinders) with per-id hue jitter
    float t = (float)((id * 2654435761u) % 997) / 997.0f;
    if (sphere) return {0.85f, 0.45f + 0.25f * t, 0.25f + 0.15f * t};
    return {0.25f + 0.15f * t, 0.55f + 0.2f * t, 0.75f};
}

static void pushInstance(Mesh& m, const Mat4& model, const Vec3& col) {
    for (int i = 0; i < 16; ++i) m.instData.push_back(model.m[i]);
    m.instData.push_back(col.x); m.instData.push_back(col.y); m.instData.push_back(col.z);
    ++m.instCount;
}
static void drawMesh(Mesh& m, int unlit) {
    if (m.instCount == 0) return;
    glBindVertexArray(m.vao);
    glBindBuffer(GL_ARRAY_BUFFER, m.inst);
    glBufferData(GL_ARRAY_BUFFER, m.instData.size() * sizeof(float), m.instData.data(), GL_DYNAMIC_DRAW);
    glUniform1i(gUUnlit, unlit);
    glDrawElementsInstanced(m.prim, m.indexCount, GL_UNSIGNED_INT, 0, m.instCount);
    glBindVertexArray(0);
}

// ------------------------------------------------------------------ frame build
static void buildFrame() {
    gSphere.instData.clear(); gSphere.instCount = 0;
    gCyl.instData.clear();    gCyl.instCount = 0;
    if (gRec.frames.empty()) return;
    const RecordedFrame& fr = gRec.frames[gFrame % gRec.numFrames()];
    for (size_t i = 0; i < gRec.actorMeta.size() && i < fr.snapshots.size(); ++i) {
        const ActorMeta& am = gRec.actorMeta[i];
        const ActorSnapshot& s = fr.snapshots[i];
        Mat4 R = Mat4::fromQuat(s.qw, s.qx, s.qy, s.qz);
        Mat4 T = Mat4::translate({s.px, s.py, s.pz});
        if (am.isSphere) {
            Mat4 M = T * R * Mat4::scale({am.radius, am.radius, am.radius});
            pushInstance(gSphere, M, palette(am.id, true));
        } else {
            Mat4 M = T * R * Mat4::scale({am.radius, am.halfHeight, am.radius});
            pushInstance(gCyl, M, palette(am.id, false));
        }
    }
}

// ------------------------------------------------------------------ GLUT cbs
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

    float R = gRec.roomSize > 0 ? gRec.roomSize : 12.0f;

    // ground (dark) at the box floor
    gGround.instData.clear(); gGround.instCount = 0;
    pushInstance(gGround, Mat4::translate({0, -R, 0}) * Mat4::scale({R * 2.5f, 1, R * 2.5f}),
                 {0.14f, 0.15f, 0.17f});
    drawMesh(gGround, 0);

    // domain wireframe (the periodic cell / reflective box)
    gWire.instData.clear(); gWire.instCount = 0;
    pushInstance(gWire, Mat4::scale({R, R, R}), {0.35f, 0.40f, 0.48f});
    drawMesh(gWire, 1);

    // bodies
    buildFrame();
    drawMesh(gSphere, 0);
    drawMesh(gCyl, 0);

    glutSwapBuffers();
}

static void reshape(int w, int h) { gWinW = w; gWinH = h; glViewport(0, 0, w, h); }

static void timer(int) {
    if (!gPaused && !gRec.frames.empty()) {
        gAccum += gSpeed;
        while (gAccum >= 1.0f) { gFrame = (gFrame + 1) % gRec.numFrames(); gAccum -= 1.0f; }
    }
    glutPostRedisplay();
    glutTimerFunc(16, timer, 0);
}

static void fitCamera() {
    float R = gRec.roomSize > 0 ? gRec.roomSize : 12.0f;
    gCenter = {0, 0, 0}; gDist = R * 2.8f; gAz = 0.7f; gEl = 0.45f;
}

static void key(unsigned char k, int, int) {
    switch (k) {
        case 27: glutLeaveMainLoop(); break;
        case ' ': gPaused = !gPaused; break;
        case '[': gFrame = (gFrame + gRec.numFrames() - 1) % gRec.numFrames(); gPaused = true; break;
        case ']': gFrame = (gFrame + 1) % gRec.numFrames(); gPaused = true; break;
        case ',': gSpeed = gSpeed > 0.1f ? gSpeed * 0.7f : gSpeed; break;
        case '.': gSpeed *= 1.4f; break;
        case 'r': gFrame = 0; break;
        case 'h': fitCamera(); break;
    }
    glutPostRedisplay();
}
static void mouse(int b, int st, int x, int y) {
    if (b == GLUT_LEFT_BUTTON) { gDrag = (st == GLUT_DOWN); gLastX = x; gLastY = y; }
    if (b == 3) gDist *= 0.9f;   // wheel up
    if (b == 4) gDist *= 1.1f;   // wheel down
    glutPostRedisplay();
}
static void motion(int x, int y) {
    if (gDrag) {
        gAz -= (x - gLastX) * 0.01f;
        gEl += (y - gLastY) * 0.01f;
        if (gEl > 1.5f) gEl = 1.5f; if (gEl < -1.5f) gEl = -1.5f;
        gLastX = x; gLastY = y; glutPostRedisplay();
    }
}

int main(int argc, char** argv) {
    if (argc < 2) { std::printf("usage: NativeViewer <file.pxrf>\n"); return 1; }
    if (!gRec.load(argv[1])) { std::fprintf(stderr, "failed to load %s\n", argv[1]); return 1; }
    std::printf("loaded %s: %u actors, %u frames, roomSize %.1f\n",
                argv[1], gRec.numActors(), gRec.numFrames(), gRec.roomSize);

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutSetOption(GLUT_MULTISAMPLE, 8);
    glutInitContextVersion(3, 3);
    glutInitContextProfile(GLUT_CORE_PROFILE);
    glutInitWindowSize(gWinW, gWinH);
    glutCreateWindow("NativeEngine Viewer  [space pause  [ ] step  , . speed  h reset-cam  drag orbit]");
    if (!loadGL()) { std::fprintf(stderr, "failed to load OpenGL 3.3 functions\n"); return 1; }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);
    gProg = link(kVert, kFrag);
    gUVP = glGetUniformLocation(gProg, "uVP");
    gUCam = glGetUniformLocation(gProg, "uCam");
    gULight = glGetUniformLocation(gProg, "uLight");
    gUUnlit = glGetUniformLocation(gProg, "uUnlit");
    makeSphere(gSphere); makeCylinder(gCyl); makeGround(gGround); makeWireBox(gWire);
    fitCamera();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(key);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutTimerFunc(16, timer, 0);
    glutMainLoop();
    return 0;
}
