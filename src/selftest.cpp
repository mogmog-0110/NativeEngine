#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "world.hpp"
#include "contact.hpp"
#include "narrowphase.hpp"
#include "gjk.hpp"

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

// ---------------------------------------------------------------------------
// 9. Sphere-vs-flat-cap-cylinder narrow phase: the three exterior regions
//    (cap face, side, rim) give the right normal and overlap.
// ---------------------------------------------------------------------------
void test_sphere_cylinder_geometry() {
    Box box; box.periodic = false; box.half = 100.0;
    Body cyl = makeCylinder(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0); // r=1, h=2

    // Cap: sphere above the +Y cap centre, just overlapping.
    Body s1 = makeSphere(1, V3{0, 2.5, 0}, 1.0, 1.0);
    Contact c1 = sphereVsCylinder(s1, cyl, box);
    check(c1.hit && closeV(c1.normal, V3{0, 1, 0}, 1e-9) && close(c1.overlap, 0.5, 1e-9),
          "sphere-cyl: cap-face contact (normal +Y, overlap 0.5)");

    // Side: sphere beside the barrel.
    Body s2 = makeSphere(2, V3{1.5, 0, 0}, 1.0, 1.0);
    Contact c2 = sphereVsCylinder(s2, cyl, box);
    check(c2.hit && closeV(c2.normal, V3{1, 0, 0}, 1e-9) && close(c2.overlap, 0.5, 1e-9),
          "sphere-cyl: side contact (normal +X, overlap 0.5)");

    // Rim: sphere off the +X/+Y edge; closest point is the rim (1,2,0).
    Body s3 = makeSphere(3, V3{1.5, 2.5, 0}, 1.0, 1.0);
    Contact c3 = sphereVsCylinder(s3, cyl, box);
    double expDist = std::sqrt(0.25 + 0.25);
    check(c3.hit && closeV(c3.point, V3{1, 2, 0}, 1e-9) &&
              close(c3.overlap, 1.0 - expDist, 1e-9),
          "sphere-cyl: rim contact (closest point on the cap edge)");

    // No contact: clearly outside.
    Body s4 = makeSphere(4, V3{0, 4.0, 0}, 1.0, 1.0);
    check(!sphereVsCylinder(s4, cyl, box).hit, "sphere-cyl: no contact when separated");
}

// ---------------------------------------------------------------------------
// 10. Sphere-vs-cylinder dynamics through the world: an off-axis strike
//     conserves momentum + angular momentum + energy and spins the cylinder.
// ---------------------------------------------------------------------------
void test_sphere_cylinder_dynamics() {
    World w;
    w.box.half = 100.0; w.box.periodic = false;
    w.restitution = 1.0; w.dt = 1e-3;
    Body cyl = makeCylinder(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
    Body sph = makeSphere(1, V3{3.0, 1.5, 0}, 1.0, 1.0);   // off the COM, beside barrel
    sph.v = V3{-2.0, 0, 0};
    w.bodies = {cyl, sph};

    auto P = [&] {
        V3 p; for (auto& b : w.bodies) p += b.v * (1.0 / b.invMass); return p; };
    auto Lm = [&] {
        V3 l; for (auto& b : w.bodies) { l += b.x.cross(b.v * (1.0 / b.invMass));
                                         l += b.angularMomentum(); } return l; };
    auto E = [&] { double e = 0; for (auto& b : w.bodies) e += b.kinetic(); return e; };

    V3 p0 = P(), l0 = Lm(); double e0 = E();
    for (int s = 0; s < 8000; ++s) w.step();

    check(closeV(P(), p0, 1e-7), "sphere-cyl dynamics: linear momentum conserved");
    check(closeV(Lm(), l0, 1e-6), "sphere-cyl dynamics: angular momentum conserved");
    check(std::fabs(E() - e0) / e0 < 1e-4, "sphere-cyl dynamics: energy conserved (e=1)");
    check(w.bodies[0].w.norm() > 1e-4, "sphere-cyl dynamics: cylinder acquires spin");
}

// ---------------------------------------------------------------------------
// 11. GJK overlap: support functions + simplex evolution. Cross-checked against
//     the independent analytic sphere-cylinder narrow phase, plus a periodic
//     across-boundary overlap.
// ---------------------------------------------------------------------------
void test_gjk_overlap() {
    Box box; box.periodic = false; box.half = 100.0;
    Body cyl = makeCylinder(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);

    Body cB = makeCylinder(1, V3{1.5, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
    check(gjkOverlap(cyl, cB, V3{0, 0, 0}), "gjk: parallel cylinders 1.5 apart overlap");
    Body cC = makeCylinder(2, V3{3.0, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
    check(!gjkOverlap(cyl, cC, V3{0, 0, 0}), "gjk: parallel cylinders 3.0 apart separated");

    // Cross-check GJK against the analytic sphere-cylinder hit test.
    struct { V3 p; bool hit; const char* what; } cases[] = {
        {{1.5, 0, 0}, true, "side"}, {{0, 2.5, 0}, true, "cap"},
        {{1.5, 2.5, 0}, true, "rim"}, {{0, 4, 0}, false, "far"},
        {{2.1, 0, 0}, false, "just outside side"},
    };
    bool agree = true;
    for (auto& c : cases) {
        Body s = makeSphere(9, c.p, 1.0, 1.0);
        bool g = gjkOverlap(s, cyl, V3{0, 0, 0});
        bool a = sphereVsCylinder(s, cyl, box).hit;
        if (g != c.hit || a != c.hit) agree = false;
    }
    check(agree, "gjk: agrees with analytic sphere-cylinder on side/cap/rim/far");

    // Periodic: two cylinders near opposite faces overlap through the boundary.
    Box pbox; pbox.periodic = true; pbox.half = 5.0;   // edge 10
    Body pA = makeCylinder(0, V3{4.5, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
    Body pB = makeCylinder(1, V3{-4.5, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
    check(gjkOverlap(pA, pB, minImageShift(pA, pB, pbox)),
          "gjk: cylinders overlap across the periodic boundary");
    Box rbox; rbox.periodic = false; rbox.half = 5.0;
    check(!gjkOverlap(pA, pB, minImageShift(pA, pB, rbox)),
          "gjk: same pair does NOT overlap without PBC (9 apart)");
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
    test_sphere_cylinder_geometry();
    test_sphere_cylinder_dynamics();
    test_gjk_overlap();
    std::printf("\n=== Summary ===\n  Passed: %d\n  Failed: %d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}

}  // namespace ne
