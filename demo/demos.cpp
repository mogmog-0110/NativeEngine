// Runnable demos, PhysX-snippet style: each builds a scene, steps it while
// recording a .pxrf for the OpenGL viewer, and prints a numerical soundness
// report (no NaN / no escape / energy behaviour / no tunneling) so breakdowns
// are caught even without watching. Shapes are spheres and cylinders (what the
// existing viewer draws).
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>

#include "world.hpp"
#include "record.hpp"

namespace ne {
namespace {

bool finiteAll(const World& w) {
    for (const Body& b : w.bodies)
        if (!(std::isfinite(b.x.x) && std::isfinite(b.x.y) && std::isfinite(b.x.z) &&
              std::isfinite(b.v.x) && std::isfinite(b.q.w)))
            return false;
    return true;
}
double maxSpeed(const World& w) {
    double m = 0; for (const Body& b : w.bodies) m = std::max(m, b.v.norm()); return m;
}
bool allInside(const World& w, double lim) {
    for (const Body& b : w.bodies)
        if (std::fabs(b.x.x) > lim || std::fabs(b.x.y) > lim || std::fabs(b.x.z) > lim)
            return false;
    return true;
}
void line() { std::printf("  ------------------------------------------------------------\n"); }

// -- Soundness demos -------------------------------------------------------

// 1. Spheres poured into a box settle into a pile without exploding or sinking.
void demo_pile(const std::string& out) {
    World w; w.box.half = 10; w.box.periodic = false;
    w.gravity = {0, -10, 0}; w.restitution = 0.2; w.friction = 0.5; w.dt = 1.0 / 240.0;
    std::mt19937_64 rng(1);
    std::uniform_real_distribution<double> p(-7, 7), h(2, 9);
    for (int i = 0; i < 80; ++i) w.bodies.push_back(makeSphere(i, {p(rng), h(rng), p(rng)}, 1.0, 1.0));
    Recorder rec; rec.init(w, 8, w.box.half);
    for (int s = 0; s <= 4000; ++s) { w.step(); rec.maybeCapture(w, s); }
    rec.save(out);
    std::printf("  pile: 80 spheres poured into a reflective box\n");
    std::printf("    finite=%d  inside=%d  maxSpeed=%.3f (settled)\n",
                finiteAll(w), allInside(w, w.box.half + 0.05), maxSpeed(w));
}

// 2. Cylinders dropped settle into a heap (exercises rotation + convex contact).
void demo_cylinders(const std::string& out) {
    World w; w.box.half = 10; w.box.periodic = false;
    w.gravity = {0, -10, 0}; w.restitution = 0.1; w.friction = 0.6; w.dt = 1.0 / 240.0;
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> p(-6, 6), h(2, 9), a(-1, 1);
    for (int i = 0; i < 40; ++i) {
        Q q{a(rng), a(rng), a(rng), a(rng)}; q = q.normalized();
        w.bodies.push_back(makeCylinder(i, {p(rng), h(rng), p(rng)}, q, 1.0, 4.0, 1.0));
    }
    Recorder rec; rec.init(w, 8, w.box.half);
    for (int s = 0; s <= 4000; ++s) { w.step(); rec.maybeCapture(w, s); }
    rec.save(out);
    std::printf("  cylinders: 40 cylinders dropped into a reflective box\n");
    std::printf("    finite=%d  inside=%d  maxSpeed=%.3f (settled)\n",
                finiteAll(w), allInside(w, w.box.half + 0.5), maxSpeed(w));
}

// 3. A hanging chain of spheres on rigid distance joints swings like a rope.
void demo_pendulum(const std::string& out) {
    World w; w.box.half = 30; w.box.periodic = false;
    w.gravity = {0, -10, 0}; w.dt = 1.0 / 240.0;
    Body anchor = makeSphere(0, {0, 12, 0}, 0.3, 1.0);
    anchor.invMass = 0; anchor.invInertiaBody = {0, 0, 0}; anchor.dynamic = false;
    w.bodies.push_back(anchor);
    const int N = 8; const double L = 1.4;
    for (int i = 1; i <= N; ++i)
        w.bodies.push_back(makeSphere(i, {double(i) * L, 12, 0}, 0.5, 1.0));  // start horizontal
    for (int i = 0; i < N; ++i) {
        DistanceJoint j; j.a = i; j.b = i + 1; j.rest = L; j.stiffness = 0;  // rigid
        w.distanceJoints.push_back(j);
    }
    Recorder rec; rec.init(w, 6, w.box.half);
    double maxErr = 0;
    for (int s = 0; s <= 4000; ++s) {
        w.step(); rec.maybeCapture(w, s);
        for (const DistanceJoint& j : w.distanceJoints)
            maxErr = std::max(maxErr, std::fabs((w.bodies[j.a].x - w.bodies[j.b].x).norm() - j.rest));
    }
    rec.save(out);
    std::printf("  pendulum: 8-link chain on rigid joints, released horizontal\n");
    std::printf("    finite=%d  max link-length error=%.4f\n", finiteAll(w), maxErr);
}

// 4. An elastic sphere gas bounces in a reflective box, conserving energy.
void demo_gas(const std::string& out) {
    World w; w.box.half = 10; w.box.periodic = false;
    w.restitution = 1.0; w.dt = 1.0 / 120.0;
    std::mt19937_64 rng(3);
    std::uniform_real_distribution<double> p(-8, 8), v(-4, 4);
    for (int i = 0; i < 100; ++i) {
        Body b = makeSphere(i, {p(rng), p(rng), p(rng)}, 0.8, 1.0);
        b.v = {v(rng), v(rng), v(rng)};
        w.bodies.push_back(b);
    }
    double e0 = w.totalKinetic();
    Recorder rec; rec.init(w, 4, w.box.half);
    for (int s = 0; s <= 3000; ++s) { w.step(); rec.maybeCapture(w, s); }
    rec.save(out);
    std::printf("  gas: 100-sphere elastic gas in a reflective box\n");
    std::printf("    finite=%d  inside=%d  energy %.2f -> %.2f (drift %.1f%%)\n",
                finiteAll(w), allInside(w, w.box.half + 0.05), e0, w.totalKinetic(),
                100.0 * std::fabs(w.totalKinetic() - e0) / e0);
}

// -- PBC demos -------------------------------------------------------------

// 5. A periodic sphere gas: bodies cross faces seamlessly, momentum + energy
//    conserved (no walls).
void demo_pbc_gas(const std::string& out) {
    World w; w.box.half = 8; w.box.periodic = true;
    w.restitution = 1.0; w.dt = 1.0 / 120.0;
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> p(-7.5, 7.5), v(-3, 3);
    for (int i = 0; i < 120; ++i) {
        Body b = makeSphere(i, {p(rng), p(rng), p(rng)}, 0.8, 1.0);
        b.v = {v(rng), v(rng), v(rng)};
        w.bodies.push_back(b);
    }
    V3 p0 = w.totalLinearMomentum(); double e0 = w.totalKinetic();
    Recorder rec; rec.init(w, 4, w.box.half);
    for (int s = 0; s <= 3000; ++s) { w.step(); rec.maybeCapture(w, s); }
    rec.save(out);
    std::printf("  pbc_gas: 120-sphere elastic gas in a PERIODIC box (wraps faces)\n");
    std::printf("    finite=%d  |dP|=%.2e  energy drift %.1f%%\n",
                finiteAll(w), (w.totalLinearMomentum() - p0).norm(),
                100.0 * std::fabs(w.totalKinetic() - e0) / e0);
}

// 6. Two spheres collide THROUGH the periodic boundary (they start near opposite
//    faces, closest approach is across the wrap).
void demo_pbc_pair(const std::string& out) {
    World w; w.box.half = 6; w.box.periodic = true;
    w.restitution = 1.0; w.dt = 1.0 / 240.0;
    Body a = makeSphere(0, {5.2, 0, 0}, 1.0, 1.0); a.v = {2, 0, 0};   // moving toward +x face
    Body b = makeSphere(1, {-5.2, 0.3, 0}, 1.0, 1.0); b.v = {-2, 0, 0}; // toward -x face
    w.bodies = {a, b};
    V3 p0 = w.totalLinearMomentum(); double e0 = w.totalKinetic();
    Recorder rec; rec.init(w, 2, w.box.half);
    double minSep = 1e9;
    for (int s = 0; s <= 2000; ++s) {
        w.step(); rec.maybeCapture(w, s);
        minSep = std::min(minSep, w.box.minImage(w.bodies[0].x - w.bodies[1].x).norm());
    }
    rec.save(out);
    std::printf("  pbc_pair: two spheres collide ACROSS the periodic boundary\n");
    std::printf("    finite=%d  |dP|=%.2e  energy drift %.1f%%  minSep=%.2f (touch~2, no tunnel)\n",
                finiteAll(w), (w.totalLinearMomentum() - p0).norm(),
                100.0 * std::fabs(w.totalKinetic() - e0) / e0, minSep);
}

}  // namespace

// Load a .pxrf back through the viewer's own loader (validates magic, version,
// and every per-frame Adler-32 CRC), confirming the file the demos wrote is one
// the OpenGL player will accept.
int verifyPxrf(const std::string& path) {
    Recording r;
    bool ok = r.load(path);
    std::printf("verify %s: load=%s  actors=%u  frames=%u\n", path.c_str(),
                ok ? "OK" : "FAILED", r.numActors(), r.numFrames());
    return (ok && r.numFrames() > 0) ? 0 : 1;
}

int runDemo(const std::string& name, const std::string& out) {
    std::printf("=== NativeEngine demo: %s -> %s ===\n", name.c_str(), out.c_str());
    line();
    if (name == "pile") demo_pile(out);
    else if (name == "cylinders") demo_cylinders(out);
    else if (name == "pendulum") demo_pendulum(out);
    else if (name == "gas") demo_gas(out);
    else if (name == "pbc_gas") demo_pbc_gas(out);
    else if (name == "pbc_pair") demo_pbc_pair(out);
    else if (name == "all") {
        demo_pile("demo_pile.pxrf");        line();
        demo_cylinders("demo_cylinders.pxrf"); line();
        demo_pendulum("demo_pendulum.pxrf"); line();
        demo_gas("demo_gas.pxrf");            line();
        demo_pbc_gas("demo_pbc_gas.pxrf");    line();
        demo_pbc_pair("demo_pbc_pair.pxrf");
    } else {
        std::printf("  unknown demo '%s'. Available: pile cylinders pendulum gas "
                    "pbc_gas pbc_pair all\n", name.c_str());
        return 1;
    }
    line();
    std::printf("  Play back with: x64\\ReleaseRender\\PhysxRender.exe %s\n", out.c_str());
    return 0;
}

}  // namespace ne
