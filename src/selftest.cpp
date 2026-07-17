#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "world.hpp"
#include "contact.hpp"

namespace ne {
namespace {

int g_pass = 0;
int g_fail = 0;

void check(bool ok, const char* name) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (ok) ++g_pass; else ++g_fail;
}

bool close(double a, double b, double tol) { return std::fabs(a - b) <= tol; }
bool closeV(const V3& a, const V3& b, double tol) {
    return (a - b).norm() <= tol;
}

// ---------------------------------------------------------------------------
// 1. Minimum-image convention: separations wrap to the nearest image.
// ---------------------------------------------------------------------------
void test_min_image() {
    Box box;
    box.half = 5.0;          // edge 10
    box.periodic = true;
    // Two points near opposite faces: raw dx = 9, min-image should be -1.
    V3 d = box.minImage(V3{9.0, 0, 0});
    check(closeV(d, V3{-1.0, 0, 0}, 1e-12), "min-image wraps 9 -> -1 (edge 10)");
    // A point just past +half wraps to just past -half.
    V3 p = box.wrap(V3{6.0, -6.0, 0.0});
    check(closeV(p, V3{-4.0, 4.0, 0.0}, 1e-12), "wrap 6 -> -4, -6 -> 4");
    // Reflective box leaves both untouched.
    box.periodic = false;
    check(closeV(box.minImage(V3{9, 0, 0}), V3{9, 0, 0}, 1e-12),
          "reflective min-image is identity");
}

// ---------------------------------------------------------------------------
// 2. Head-on equal-mass elastic collision: velocities exchange; p, E conserved.
// ---------------------------------------------------------------------------
void test_elastic_headon() {
    World w;
    w.box.half = 100.0;
    w.box.periodic = false;
    w.restitution = 1.0;
    w.dt = 1e-3;
    w.contactBeta = 0.2;
    Body a = makeSphere(0, V3{-2.0, 0, 0}, 1.0, 1.0);
    Body b = makeSphere(1, V3{+2.0, 0, 0}, 1.0, 1.0);
    a.v = V3{+1.0, 0, 0};
    b.v = V3{-1.0, 0, 0};
    w.bodies = {a, b};

    V3 p0 = w.totalLinearMomentum();
    double e0 = w.totalKinetic();
    for (int s = 0; s < 20000; ++s) w.step();

    V3 p1 = w.totalLinearMomentum();
    double e1 = w.totalKinetic();
    check(closeV(p0, p1, 1e-9), "head-on: linear momentum conserved");
    check(close(e0, e1, 1e-6), "head-on: kinetic energy conserved (e=1)");
    // After an elastic exchange the two are separating.
    double sep = (w.bodies[0].v.x < w.bodies[1].v.x);
    check(sep, "head-on: spheres separate after collision");
    check(w.bodies[0].v.x < 0 && w.bodies[1].v.x > 0,
          "head-on: equal masses exchange velocity sign");
}

// ---------------------------------------------------------------------------
// 3. Periodic transit: a lone sphere crosses a face and reappears, velocity
//    unchanged, trajectory continuous under the minimum image.
// ---------------------------------------------------------------------------
void test_periodic_transit() {
    World w;
    w.box.half = 5.0;
    w.box.periodic = true;
    w.dt = 1e-2;
    Body a = makeSphere(0, V3{4.9, 0, 0}, 0.5, 1.0);
    a.v = V3{1.0, 0, 0};
    w.bodies = {a};

    V3 v0 = w.bodies[0].v;
    // Track unwrapped position by accumulating min-image steps.
    V3 prev = w.bodies[0].x;
    V3 unwrapped = prev;
    for (int s = 0; s < 300; ++s) {
        w.step();
        V3 cur = w.bodies[0].x;
        unwrapped += w.box.minImage(cur - prev);
        prev = cur;
    }
    check(closeV(w.bodies[0].v, v0, 1e-12), "periodic: velocity unchanged crossing face");
    // 300 steps * 0.01 * v=1 = 3.0 of travel from x=4.9.
    check(close(unwrapped.x, 4.9 + 3.0, 1e-9), "periodic: unwrapped trajectory continuous");
    check(std::fabs(w.bodies[0].x.x) <= 5.0 + 1e-9, "periodic: stays in primary cell");
}

// ---------------------------------------------------------------------------
// 4. Torque-free spin: quaternion stays unit; isotropic body conserves omega and
//    angular momentum exactly.
// ---------------------------------------------------------------------------
void test_free_spin() {
    World w;
    w.box.half = 100.0;
    w.box.periodic = false;
    w.dt = 1e-3;
    Body a = makeSphere(0, V3{0, 0, 0}, 1.0, 1.0);
    a.w = V3{0.3, -0.7, 1.1};
    w.bodies = {a};

    V3 L0 = w.bodies[0].angularMomentum();
    for (int s = 0; s < 5000; ++s) w.step();
    double qn = w.bodies[0].q.norm();
    check(close(qn, 1.0, 1e-9), "free spin: quaternion stays unit norm");
    check(closeV(w.bodies[0].angularMomentum(), L0, 1e-9),
          "free spin: angular momentum conserved (isotropic)");
    check(closeV(w.bodies[0].w, V3{0.3, -0.7, 1.1}, 1e-9),
          "free spin: omega constant for a sphere");
}

// ---------------------------------------------------------------------------
// 5. Reflective wall: an N-body elastic gas stays inside and conserves energy.
// ---------------------------------------------------------------------------
void test_wall_gas() {
    World w;
    w.box.half = 10.0;
    w.box.periodic = false;
    w.restitution = 1.0;
    w.dt = 5e-3;
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> pos(-8.0, 8.0), vel(-3.0, 3.0);
    for (int i = 0; i < 60; ++i) {
        Body b = makeSphere(i, V3{pos(rng), pos(rng), pos(rng)}, 1.0, 1.0);
        b.v = V3{vel(rng), vel(rng), vel(rng)};
        w.bodies.push_back(b);
    }
    double e0 = w.totalKinetic();
    bool inside = true;
    for (int s = 0; s < 4000; ++s) {
        w.step();
        for (const Body& b : w.bodies)
            if (std::fabs(b.x.x) > 10.0 + 1e-6 || std::fabs(b.x.y) > 10.0 + 1e-6 ||
                std::fabs(b.x.z) > 10.0 + 1e-6)
                inside = false;
    }
    double e1 = w.totalKinetic();
    check(inside, "wall gas: all spheres remain inside the box");
    // Impulsive walls + positional correction leak a little; allow 1%.
    check(std::fabs(e1 - e0) / e0 < 0.01, "wall gas: kinetic energy conserved to 1%");
}

// ---------------------------------------------------------------------------
// 6. Determinism: the same periodic gas run twice is bit-for-bit identical.
// ---------------------------------------------------------------------------
World makeGas(bool periodic) {
    World w;
    w.box.half = 8.0;
    w.box.periodic = periodic;
    w.restitution = 1.0;
    w.dt = 5e-3;
    std::mt19937_64 rng(999);
    std::uniform_real_distribution<double> pos(-6.0, 6.0), vel(-2.0, 2.0);
    for (int i = 0; i < 50; ++i) {
        Body b = makeSphere(i, V3{pos(rng), pos(rng), pos(rng)}, 1.0, 1.0);
        b.v = V3{vel(rng), vel(rng), vel(rng)};
        b.w = V3{vel(rng), vel(rng), vel(rng)};
        w.bodies.push_back(b);
    }
    return w;
}

void test_determinism() {
    World a = makeGas(true), b = makeGas(true);
    for (int s = 0; s < 3000; ++s) { a.step(); b.step(); }
    bool identical = true;
    for (size_t i = 0; i < a.bodies.size(); ++i) {
        const Body& ba = a.bodies[i];
        const Body& bb = b.bodies[i];
        if (ba.x.x != bb.x.x || ba.x.y != bb.x.y || ba.x.z != bb.x.z ||
            ba.v.x != bb.v.x || ba.q.w != bb.q.w) {
            identical = false;
            break;
        }
    }
    check(identical, "determinism: identical runs are bit-for-bit equal");
}

// ---------------------------------------------------------------------------
// 7. Periodic momentum: a periodic gas has no walls, so linear momentum and
//    kinetic energy are conserved (to a small integration tolerance).
// ---------------------------------------------------------------------------
void test_periodic_conservation() {
    World w = makeGas(true);
    V3 p0 = w.totalLinearMomentum();
    double e0 = w.totalKinetic();
    for (int s = 0; s < 3000; ++s) w.step();
    check(closeV(w.totalLinearMomentum(), p0, 1e-6),
          "periodic gas: linear momentum conserved");
    check(std::fabs(w.totalKinetic() - e0) / e0 < 0.02,
          "periodic gas: kinetic energy conserved to 2%");
}

// ---------------------------------------------------------------------------
// 8. Off-centre contact impulse: conserves linear momentum, angular momentum
//    (about the origin), and energy at e=1, AND actually imparts spin. This is
//    the foundation for flat-cap cylinder contacts.
// ---------------------------------------------------------------------------
void test_offcentre_impulse() {
    // A cylinder struck by a fast sphere near one cap: the sphere hits off the
    // cylinder's centre of mass, so the cylinder must both translate and spin.
    Body cyl = makeCylinder(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
    Body sph = makeSphere(1, V3{2.0, 1.8, 0}, 1.0, 1.0);   // near the +Y cap, off axis
    sph.v = V3{-1.0, 0, 0};                                 // moving into the cylinder side

    auto totalP = [&] { return cyl.v * (1.0 / cyl.invMass) + sph.v * (1.0 / sph.invMass); };
    auto totalL = [&] {
        return cyl.x.cross(cyl.v * (1.0 / cyl.invMass)) + cyl.angularMomentum()
             + sph.x.cross(sph.v * (1.0 / sph.invMass)) + sph.angularMomentum();
    };
    auto totalE = [&] { return cyl.kinetic() + sph.kinetic(); };

    V3 p0 = totalP(), l0 = totalL();
    double e0 = totalE();

    // Contact on the cylinder's +x side at y=1.8. The resolver's normal points
    // from b (sphere) to a (cylinder), i.e. -x; the impulse then pushes the
    // cylinder away from the sphere.
    V3 n{-1, 0, 0};
    V3 contact{1.0, 1.8, 0.0};                 // on the cylinder surface
    V3 ra = contact - cyl.x;                    // off-centre -> lever arm
    V3 rb = contact - sph.x;
    resolveContact(cyl, sph, ra, rb, n, 1.0);

    check(closeV(totalP(), p0, 1e-12), "off-centre: linear momentum conserved");
    check(closeV(totalL(), l0, 1e-12), "off-centre: angular momentum conserved (about origin)");
    check(close(totalE(), e0, 1e-9), "off-centre: kinetic energy conserved (e=1)");
    check(cyl.w.norm() > 1e-6, "off-centre: cylinder acquires spin");
    // Impulse J = j*(-x) applied at r=(1,1.8,0): torque r x J = (0,0,1.8 j), so
    // the strike spins the cylinder about +z.
    check(std::fabs(cyl.w.z) > std::fabs(cyl.w.x) && std::fabs(cyl.w.z) > std::fabs(cyl.w.y),
          "off-centre: spin is about the expected (z) axis");
}

}  // namespace

int runSelftest() {
    std::printf("=== NativeEngine physics selftest ===\n");
    test_min_image();
    test_elastic_headon();
    test_periodic_transit();
    test_free_spin();
    test_wall_gas();
    test_determinism();
    test_periodic_conservation();
    test_offcentre_impulse();
    std::printf("\n=== Summary ===\n  Passed: %d\n  Failed: %d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}

}  // namespace ne
