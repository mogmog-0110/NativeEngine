#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "world.hpp"
#include "contact.hpp"
#include "narrowphase.hpp"
#include "gjk.hpp"
#include "epa.hpp"
#include "physics_world.hpp"
#include "character.hpp"

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
// 7b. Translation equivariance: a periodic world has NO special origin. Shift
//     every body by an arbitrary vector (wrapped) and the dynamics must be
//     identical -- positions shift by the same vector, velocities unchanged.
//     This is the defining symmetry of PBC; a position-dependent bug (e.g. a
//     wrap/min-image mismatch) breaks it even when momentum is still conserved.
//     Brute-force pairing so the constraint order is position-independent, so
//     the check holds to floating-point tolerance rather than being masked.
// ---------------------------------------------------------------------------
void test_periodic_translation_invariance() {
    const V3 s{3.1, -2.7, 5.3};
    World w = makeGas(true);
    w.forceBruteForce = true;
    World w2 = w;
    for (Body& b : w2.bodies) b.x = w2.box.wrap(b.x + s);
    double maxErr = 0.0;
    for (int step = 0; step < 40; ++step) {
        w.step();
        w2.step();
        for (size_t i = 0; i < w.bodies.size(); ++i) {
            V3 expected = w.box.wrap(w.bodies[i].x + s);
            maxErr = std::max(maxErr, w.box.minImage(w2.bodies[i].x - expected).norm());
            maxErr = std::max(maxErr, (w2.bodies[i].v - w.bodies[i].v).norm());
        }
    }
    check(maxErr < 1e-6, "periodic: dynamics invariant under a global shift");
}

// ---------------------------------------------------------------------------
// 7c. A collision straddling the boundary must be IDENTICAL to the same
//     collision in free space -- the minimum-image contact is real physics,
//     not an approximation. Two equal spheres meet head-on across the +x/-x
//     faces; e=1 exchanges their velocities exactly, matching a free-space run.
// ---------------------------------------------------------------------------
void test_cross_boundary_collision() {
    auto run = [](bool periodic) {
        World w;
        w.dt = 1.0 / 240.0; w.restitution = 1.0; w.box.periodic = periodic;
        w.box.half = periodic ? 5.0 : 50.0;
        Body a = makeSphere(0, periodic ? V3{4.4, 0, 0} : V3{-0.6, 0, 0}, 0.5, 1.0);
        Body b = makeSphere(1, periodic ? V3{-4.4, 0, 0} : V3{0.6, 0, 0}, 0.5, 1.0);
        a.v = V3{1.0, 0, 0};
        b.v = V3{-1.0, 0, 0};
        w.bodies = {a, b};
        for (int step = 0; step < 400; ++step) w.step();
        return w;
    };
    World fp = run(false), pr = run(true);
    check(close(pr.bodies[0].v.x, fp.bodies[0].v.x, 1e-9) &&
          close(pr.bodies[1].v.x, fp.bodies[1].v.x, 1e-9),
          "periodic: cross-boundary collision equals free space");
    check(close(pr.bodies[0].v.x, -1.0, 1e-6) && close(pr.bodies[1].v.x, 1.0, 1e-6),
          "periodic: cross-boundary head-on exchanges velocity");
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

// ---------------------------------------------------------------------------
// 12. EPA penetration: normal + depth cross-checked against the INDEPENDENT
//     analytic sphere-cylinder oracle, plus contact-point sanity.
// ---------------------------------------------------------------------------
void test_epa_contact() {
    Box box; box.periodic = false; box.half = 100.0;
    Body cyl = makeCylinder(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);

    // Sphere overlapping the side/cap: EPA (sphere=A, cyl=B) normal points from
    // B (cyl) toward A (sphere) -- outward from the cylinder -- matching the
    // analytic oracle.
    struct { V3 p; V3 n; double ov; const char* what; } cs[] = {
        {{1.6, 0, 0},   {1, 0, 0}, 0.4, "side"},
        {{0, 2.7, 0},   {0, 1, 0}, 0.3, "cap"},
    };
    bool ok = true;
    for (auto& t : cs) {
        Body s = makeSphere(9, t.p, 1.0, 1.0);
        Contact e = convexContact(s, cyl, box);       // A=sphere, normal B->A = outward
        Contact a = sphereVsCylinder(s, cyl, box);
        if (!e.hit || !a.hit) { ok = false; continue; }
        if ((e.normal - t.n).norm() > 1e-3) ok = false;
        if (std::fabs(e.overlap - t.ov) > 1e-3) ok = false;
        if (std::fabs(e.overlap - a.overlap) > 1e-3) ok = false;
    }
    check(ok, "epa: normal+depth match analytic sphere-cylinder (side, cap)");

    // Two parallel cylinders overlapping radially by 0.5: normal ~ +-x, depth 0.5.
    Body c2 = makeCylinder(1, V3{1.5, 0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
    Contact cc = convexContact(cyl, c2, box);          // A=cyl@0, B=cyl@1.5 -> n = -x
    check(cc.hit && std::fabs(std::fabs(cc.normal.x) - 1.0) < 1e-3 &&
              std::fabs(cc.overlap - 0.5) < 1e-3,
          "epa: parallel cylinders -- normal along x, depth 0.5");
    check(cc.hit && std::fabs(cc.point.x - 0.75) < 0.15,
          "epa: cylinder-cylinder contact point between the two");
}

// ---------------------------------------------------------------------------
// 13. Cylinder-cylinder dynamics through the world: conserves momentum +
//     angular momentum + energy, across a reflective and a periodic boundary.
// ---------------------------------------------------------------------------
void test_cylinder_cylinder_dynamics() {
    for (bool periodic : {false, true}) {
        World w;
        w.box.half = periodic ? 6.0 : 100.0;
        w.box.periodic = periodic;
        w.restitution = 1.0; w.dt = 1e-3;
        Body a = makeCylinder(0, V3{-2.5, 0.3, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
        Body b = makeCylinder(1, V3{ 2.5, 0.0, 0}, Q{1, 0, 0, 0}, 1.0, 4.0, 1.0);
        a.v = V3{ 1.5, 0, 0};
        b.v = V3{-1.5, 0, 0};
        w.bodies = {a, b};

        auto P = [&] { V3 p; for (auto& x : w.bodies) p += x.v * (1.0 / x.invMass); return p; };
        auto Lm = [&] { V3 l; for (auto& x : w.bodies) {
            l += x.x.cross(x.v * (1.0 / x.invMass)); l += x.angularMomentum(); } return l; };
        auto E = [&] { double e = 0; for (auto& x : w.bodies) e += x.kinetic(); return e; };

        V3 p0 = P(); double e0 = E();
        V3 l0 = Lm();
        double minSep = 1e9;
        for (int s = 0; s < 6000; ++s) {
            w.step();
            double sep = w.box.minImage(w.bodies[0].x - w.bodies[1].x).norm();
            if (sep < minSep) minSep = sep;
        }

        check(closeV(P(), p0, 1e-6),
              periodic ? "cyl-cyl: momentum conserved (periodic)"
                       : "cyl-cyl: momentum conserved (reflective)");
        check(std::fabs(E() - e0) / e0 < 1e-3,
              periodic ? "cyl-cyl: energy conserved (periodic)"
                       : "cyl-cyl: energy conserved (reflective)");
        // The collision MUST actually happen and MUST NOT tunnel. Two parallel
        // r=1 barrels touch at centre separation ~2; a wrong-sign normal would
        // let them pass through, driving minSep toward 0. So a real elastic
        // collision keeps minSep in a band just above 2, never near 0.
        check(minSep > 1.7 && minSep < 2.3,
              periodic ? "cyl-cyl: real collision, no tunneling (periodic)"
                       : "cyl-cyl: real collision, no tunneling (reflective)");
        if (!periodic)
            check(closeV(Lm(), l0, 1e-5), "cyl-cyl: angular momentum conserved (reflective)");
    }
}

// ---------------------------------------------------------------------------
// 14. Box shape: EPA geometry for an axis-aligned overlap, and dynamic box-box
//     + sphere-box collisions conserve momentum + energy and don't tunnel.
// ---------------------------------------------------------------------------
void test_box() {
    Box box; box.periodic = false; box.half = 100.0;

    // Two unit boxes (half-extent 1) overlapping along x by 0.5: normal ~x, depth 0.5.
    Body b1 = makeBox(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    Body b2 = makeBox(1, V3{1.5, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    Contact bb = convexContact(b1, b2, box);            // A=b1, normal B->A = -x
    check(bb.hit && std::fabs(std::fabs(bb.normal.x) - 1.0) < 1e-3 &&
              std::fabs(bb.overlap - 0.5) < 1e-3,
          "box: axis-aligned overlap -- normal along x, depth 0.5");

    // Sphere overlapping a box face.
    Body sB = makeBox(2, V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    Body sS = makeSphere(3, V3{1.6, 0, 0}, 1.0, 1.0);   // face at x=1, sphere reaches 0.6
    Contact sb = convexContact(sS, sB, box);
    check(sb.hit && (sb.normal - V3{1, 0, 0}).norm() < 1e-2 &&
              std::fabs(sb.overlap - 0.4) < 1e-2,
          "box: sphere-box face contact (normal +x, depth 0.4)");

    // Dynamic box-box head-on: conserve p + E, actually collide (no tunnel).
    World w;
    w.box.half = 100.0; w.box.periodic = false; w.restitution = 1.0; w.dt = 1e-3;
    Body a = makeBox(0, V3{-2.5, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    Body b = makeBox(1, V3{ 2.5, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    a.v = V3{ 1.5, 0, 0}; b.v = V3{-1.5, 0, 0};
    w.bodies = {a, b};
    auto P = [&] { V3 p; for (auto& x : w.bodies) p += x.v * (1.0 / x.invMass); return p; };
    auto E = [&] { double e = 0; for (auto& x : w.bodies) e += x.kinetic(); return e; };
    V3 p0 = P(); double e0 = E(); double minSep = 1e9;
    for (int s = 0; s < 6000; ++s) {
        w.step();
        double sep = (w.bodies[0].x - w.bodies[1].x).norm();
        if (sep < minSep) minSep = sep;
    }
    check(closeV(P(), p0, 1e-6), "box-box: linear momentum conserved");
    check(std::fabs(E() - e0) / e0 < 1e-3, "box-box: energy conserved");
    check(minSep > 1.7 && minSep < 2.3, "box-box: real collision, no tunneling");
}

// ---------------------------------------------------------------------------
// 15. Coulomb friction at a contact: a glancing hit loses tangential speed,
//     gains spin, dissipates energy, and stays inside the friction cone;
//     mu = 0 leaves tangential motion untouched.
// ---------------------------------------------------------------------------
void test_friction() {
    // Ball moving +x and slightly down onto a static floor (infinite mass box).
    auto setup = [] {
        Body ball = makeSphere(0, V3{0, 1.0, 0}, 1.0, 1.0);
        ball.v = V3{2.0, -1.0, 0};   // sliding +x while pressing down
        return ball;
    };
    Body floor = makeBox(1, V3{0, -1.0, 0}, Q{1, 0, 0, 0}, V3{10, 1, 10}, 1.0);
    floor.invMass = 0.0; floor.invInertiaBody = {0, 0, 0};  // static

    // Contact at the ball's bottom, normal +y (floor -> ball).
    V3 contact{0, 0, 0};
    V3 n{0, 1, 0};

    Body frictionless = setup();
    resolveContact(frictionless, floor, contact - frictionless.x, contact - floor.x, n, 0.0, 0.0);
    check(std::fabs(frictionless.v.x - 2.0) < 1e-9 && frictionless.w.norm() < 1e-9,
          "friction: mu=0 leaves tangential velocity and spin untouched");

    Body rough = setup();
    double vx0 = rough.v.x, ke0 = rough.kinetic();
    resolveContact(rough, floor, contact - rough.x, contact - floor.x, n, 0.0, 0.5);
    check(rough.v.x < vx0 - 1e-6, "friction: tangential velocity is reduced");
    check(rough.w.norm() > 1e-6, "friction: contact imparts spin (rolling)");
    check(rough.kinetic() <= ke0 + 1e-12, "friction: does not add energy");
    // Cone: tangential impulse magnitude <= mu * normal impulse. The ball's
    // normal impulse stopped v.y from -1 to 0 (dp_n = m*1). Tangential dp_x =
    // m*(vx0 - v.x) must be <= mu * dp_n.
    double m = 1.0 / rough.invMass;
    double dpN = m * 1.0, dpT = m * (vx0 - rough.v.x);
    check(dpT <= 0.5 * dpN + 1e-6, "friction: stays inside the Coulomb cone");
}

// ---------------------------------------------------------------------------
// 16. Gravity: a dynamic body accelerates at g; a static body does not.
// ---------------------------------------------------------------------------
void test_gravity() {
    World w;
    w.box.half = 1e6; w.box.periodic = false;
    w.gravity = V3{0, -10.0, 0};
    w.dt = 1e-3;
    Body a = makeSphere(0, V3{0, 0, 0}, 1.0, 1.0);
    w.bodies = {a};
    for (int s = 0; s < 1000; ++s) w.step();          // 1 second
    check(std::fabs(w.bodies[0].v.y - (-10.0)) < 1e-6, "gravity: v.y = -g*t after 1 s");
    check(w.bodies[0].x.y < -4.0 && w.bodies[0].x.y > -6.0, "gravity: fell ~5 m in 1 s");
}

// ---------------------------------------------------------------------------
// 17. Handle facade: creation, handle stability across removal, pose get/set,
//     impulse, and the contact callback.
// ---------------------------------------------------------------------------
void test_physics_world() {
    PhysicsWorld pw;
    pw.setBox(100.0, false);
    pw.setTimestep(1e-3);

    BodyId s1 = pw.addSphere(V3{-5, 0, 0}, 1.0, 1.0);
    BodyId s2 = pw.addSphere(V3{0, 0, 0}, 1.0, 1.0);
    BodyId s3 = pw.addSphere(V3{5, 0, 0}, 1.0, 1.0);
    check(s1 && s2 && s3 && s1 != s2 && s2 != s3, "facade: distinct valid handles");
    check(pw.bodyCount() == 3, "facade: three bodies");

    // Remove the middle one; the others' handles must stay valid and correct.
    pw.remove(s2);
    check(!pw.valid(s2) && pw.valid(s1) && pw.valid(s3), "facade: removal keeps other handles valid");
    check(closeV(pw.position(s1), V3{-5, 0, 0}, 1e-12) &&
              closeV(pw.position(s3), V3{5, 0, 0}, 1e-12),
          "facade: swap-and-pop preserves poses by handle");
    check(pw.bodyCount() == 2, "facade: two bodies after removal");

    // Pose set/get + impulse.
    pw.setPosition(s1, V3{1, 2, 3});
    check(closeV(pw.position(s1), V3{1, 2, 3}, 1e-12), "facade: setPosition/position round-trip");
    pw.applyImpulse(s3, V3{0, 0, 0});   // no-op sanity
    pw.applyImpulse(s1, V3{4.18879, 0, 0});   // mass = 4/3 pi ~ 4.18879 -> v.x = 1
    check(std::fabs(pw.velocity(s1).x - 1.0) < 1e-4, "facade: applyImpulse gives v = J/m");

    // Contact callback fires for an overlapping pair.
    PhysicsWorld cw;
    cw.setBox(100.0, false);
    cw.setTimestep(1e-3);
    BodyId a = cw.addSphere(V3{-0.6, 0, 0}, 1.0, 1.0);   // overlapping (centres 1.2 < 2)
    BodyId b = cw.addSphere(V3{ 0.6, 0, 0}, 1.0, 1.0);
    int hits = 0; BodyId ga = 0, gb = 0;
    cw.setContactCallback([&](const ContactInfo& ci) { ++hits; ga = ci.a; gb = ci.b; });
    cw.step();
    check(hits >= 1 && ((ga == a && gb == b) || (ga == b && gb == a)),
          "facade: contact callback fires for the overlapping pair");
}

// ---------------------------------------------------------------------------
// 18. Resting contact: a sphere dropped onto a static floor under gravity comes
//     to rest at radius height and stops bouncing (restitution slop + position
//     correction + the iterative solver).
// ---------------------------------------------------------------------------
void test_resting() {
    PhysicsWorld pw;
    pw.setBox(1e6, false);
    pw.setTimestep(1.0 / 240.0);
    pw.setGravity(V3{0, -10.0, 0});
    pw.setRestitution(0.0);   // inelastic -> should settle
    BodyId floor = pw.makeStatic(
        pw.addBox(V3{0, -1.0, 0}, Q{1, 0, 0, 0}, V3{20, 1, 20}, 1.0));  // top at y=0
    BodyId ball = pw.addSphere(V3{0, 3.0, 0}, 1.0, 1.0);
    (void)floor;
    for (int s = 0; s < 3000; ++s) pw.step();   // 12.5 s
    double y = pw.position(ball).y;
    double vy = pw.velocity(ball).y;
    check(y > 0.9 && y < 1.15, "resting: sphere settles at ~radius above the floor");
    check(std::fabs(vy) < 0.2, "resting: sphere stops bouncing (small residual velocity)");
    check(y > 0.5, "resting: sphere does not sink through the static floor");
}

// ---------------------------------------------------------------------------
// 19. Box-box manifold: a box rests FLAT on a static floor (does not tip) and a
//     two-box stack stays stacked. Impossible with single-point contact.
// ---------------------------------------------------------------------------
void test_stacking() {
    // A box dropped flat onto a static floor should settle flat, not tip.
    {
        PhysicsWorld pw;
        pw.setBox(1e6, false);
        pw.setTimestep(1.0 / 240.0);
        pw.setGravity(V3{0, -10, 0});
        pw.setRestitution(0.0);
        pw.setFriction(0.5);
        pw.makeStatic(pw.addBox(V3{0, -1, 0}, Q{1, 0, 0, 0}, V3{20, 1, 20}, 1.0)); // top y=0
        BodyId bx = pw.addBox(V3{0, 1.6, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);      // falls to y=1
        for (int s = 0; s < 3000; ++s) pw.step();
        double y = pw.position(bx).y;
        double spin = pw.angularVelocity(bx).norm();
        check(y > 0.9 && y < 1.15, "stacking: box settles flat at rest height");
        check(spin < 0.1, "stacking: box does NOT tip (angular velocity ~0)");
    }
    // Two boxes stacked stay stacked (upper near y=3, both quiet).
    {
        PhysicsWorld pw;
        pw.setBox(1e6, false);
        pw.setTimestep(1.0 / 240.0);
        pw.setGravity(V3{0, -10, 0});
        pw.setRestitution(0.0);
        pw.setFriction(0.6);
        pw.makeStatic(pw.addBox(V3{0, -1, 0}, Q{1, 0, 0, 0}, V3{20, 1, 20}, 1.0)); // top y=0
        BodyId lo = pw.addBox(V3{0, 1.05, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
        BodyId hi = pw.addBox(V3{0, 3.10, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
        for (int s = 0; s < 4000; ++s) pw.step();
        double yl = pw.position(lo).y, yh = pw.position(hi).y;
        check(yl > 0.85 && yl < 1.2, "stacking: lower box rests on floor (~y=1)");
        check(yh > 2.7 && yh < 3.3, "stacking: upper box rests on lower (~y=3)");
        check(pw.velocity(hi).norm() < 0.3, "stacking: stack is at rest");
    }
}

// ---------------------------------------------------------------------------
// 20. Sleeping: a settled body sleeps and stops moving; an impulse wakes it; a
//     falling body wakes a sleeper it lands on.
// ---------------------------------------------------------------------------
void test_sleeping() {
    PhysicsWorld pw;
    pw.setBox(1e6, false);
    pw.setTimestep(1.0 / 240.0);
    pw.setGravity(V3{0, -10, 0});
    pw.setRestitution(0.0);
    pw.setSleepEnabled(true);
    pw.makeStatic(pw.addBox(V3{0, -1, 0}, Q{1, 0, 0, 0}, V3{20, 1, 20}, 1.0));
    BodyId ball = pw.addSphere(V3{0, 3.0, 0}, 1.0, 1.0);

    for (int s = 0; s < 2000; ++s) pw.step();   // settle + sleep
    check(pw.isSleeping(ball), "sleeping: a settled body falls asleep");
    V3 restPos = pw.position(ball);
    for (int s = 0; s < 2000; ++s) pw.step();   // stays put while asleep
    check((pw.position(ball) - restPos).norm() < 1e-6, "sleeping: an asleep body does not move");

    // An impulse wakes it.
    pw.applyImpulse(ball, V3{5, 0, 0});
    check(!pw.isSleeping(ball), "sleeping: applyImpulse wakes the body");
    pw.step();
    check(pw.velocity(ball).x > 0.1, "sleeping: woken body responds to the impulse");

    // A falling body wakes a sleeper it lands on.
    PhysicsWorld pw2;
    pw2.setBox(1e6, false); pw2.setTimestep(1.0 / 240.0);
    pw2.setGravity(V3{0, -10, 0}); pw2.setRestitution(0.0); pw2.setSleepEnabled(true);
    pw2.makeStatic(pw2.addBox(V3{0, -1, 0}, Q{1, 0, 0, 0}, V3{20, 1, 20}, 1.0));
    BodyId low = pw2.addBox(V3{0, 1.0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    for (int s = 0; s < 2000; ++s) pw2.step();
    check(pw2.isSleeping(low), "sleeping: the lower box sleeps");
    BodyId drop = pw2.addBox(V3{0, 6.0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    (void)drop;
    bool woke = false;
    for (int s = 0; s < 2000; ++s) { pw2.step(); if (!pw2.isSleeping(low)) { woke = true; break; } }
    check(woke, "sleeping: a falling box wakes the sleeper it lands on");
}

// ---------------------------------------------------------------------------
// 21. Broadphase: the uniform grid gives BIT-IDENTICAL results to the O(N^2)
//     scan (it only prunes pairs that cannot touch), in a periodic box.
// ---------------------------------------------------------------------------
void test_broadphase() {
    auto build = [](bool brute) {
        World w;
        w.box.half = 12.0; w.box.periodic = true;
        w.restitution = 1.0; w.dt = 5e-3;
        w.forceBruteForce = brute;
        std::mt19937_64 rng(2024);
        std::uniform_real_distribution<double> pos(-11.0, 11.0), vel(-2.0, 2.0);
        for (int i = 0; i < 200; ++i) {          // > 64 -> grid path when !brute
            Body b = makeSphere(i, V3{pos(rng), pos(rng), pos(rng)}, 1.0, 1.0);
            b.v = V3{vel(rng), vel(rng), vel(rng)};
            b.w = V3{vel(rng), vel(rng), vel(rng)};
            w.bodies.push_back(b);
        }
        return w;
    };
    World grid = build(false), brute = build(true);
    for (int s = 0; s < 1500; ++s) { grid.step(); brute.step(); }

    bool identical = true;
    double e = 0;
    for (size_t i = 0; i < grid.bodies.size(); ++i) {
        const Body& g = grid.bodies[i];
        const Body& b = brute.bodies[i];
        if (g.x.x != b.x.x || g.x.y != b.x.y || g.x.z != b.x.z ||
            g.v.x != b.v.x || g.q.w != b.q.w) { identical = false; }
        e = std::max(e, (g.x - b.x).norm());
    }
    check(identical, "broadphase: grid result is bit-identical to brute force");
    check(e == 0.0, "broadphase: max position difference is exactly zero");
}

// ---------------------------------------------------------------------------
// 22. Distance joints: a rigid rod keeps a pendulum's length constant while it
//     swings; a spring-damper settles at the stretched equilibrium.
// ---------------------------------------------------------------------------
void test_joints() {
    // Rigid pendulum: bob hangs from a static anchor by a rigid rod of length 3.
    {
        PhysicsWorld pw;
        pw.setBox(1e6, false);
        pw.setTimestep(1.0 / 240.0);
        pw.setGravity(V3{0, -10, 0});
        BodyId anchor = pw.makeStatic(pw.addSphere(V3{0, 5, 0}, 0.1, 1.0));
        BodyId bob = pw.addSphere(V3{3, 5, 0}, 1.0, 1.0);   // 3 to the side
        pw.addDistanceJoint(anchor, bob, V3{0, 0, 0}, V3{0, 0, 0}, 3.0);  // rigid

        double maxErr = 0, minY = 5.0;
        for (int s = 0; s < 2000; ++s) {
            pw.step();
            double len = (pw.position(bob) - pw.position(anchor)).norm();
            maxErr = std::max(maxErr, std::fabs(len - 3.0));
            minY = std::min(minY, pw.position(bob).y);
        }
        check(maxErr < 0.15, "joint: rigid rod keeps pendulum length ~constant");
        check(minY < 3.0, "joint: pendulum actually swings down under gravity");
    }
    // Spring: bob hangs on a spring; settles where k*(dist-rest) = m*g.
    {
        PhysicsWorld pw;
        pw.setBox(1e6, false);
        pw.setTimestep(1.0 / 240.0);
        pw.setGravity(V3{0, -10, 0});
        BodyId anchor = pw.makeStatic(pw.addSphere(V3{0, 5, 0}, 0.1, 1.0));
        BodyId bob = pw.addSphere(V3{0, 3, 0}, 1.0, 1.0);   // 2 below
        double m = 1.0 / 1.0;   // density 1, r=1 -> mass 4/3 pi; use actual
        (void)m;
        double k = 200.0;
        pw.addDistanceJoint(anchor, bob, V3{0, 0, 0}, V3{0, 0, 0}, 2.0, k, 20.0);

        for (int s = 0; s < 6000; ++s) pw.step();
        double len = (pw.position(bob) - pw.position(anchor)).norm();
        double mass = 1.0 / 1.0;   // recompute exact below
        // exact bob mass = density * 4/3 pi r^3 = 4.18879
        mass = 4.18879;
        double expectedStretch = mass * 10.0 / k;    // k*x = m g
        check(std::fabs((len - 2.0) - expectedStretch) < 0.05,
              "joint: spring settles at the gravity-balanced stretch");
        check(pw.velocity(bob).norm() < 0.1, "joint: damped spring comes to rest");
    }
}

// ---------------------------------------------------------------------------
// Game-backend body properties (Tier 1): per-body material + combine modes,
// damping, collision layers/masks, and userData routed through contact events.
// ---------------------------------------------------------------------------
void test_per_body_material() {
    // Head-on of two equal spheres; return their separating speed after impact.
    auto headon = [](double ea, double eb, World::Combine mode) {
        World w;
        w.dt = 1.0 / 240.0; w.restitution = 0.5; w.box.half = 100.0;
        w.restitutionCombine = mode;
        Body a = makeSphere(0, V3{-1.0, 0, 0}, 0.5, 1.0); a.v = V3{1, 0, 0}; a.restitution = ea;
        Body b = makeSphere(1, V3{1.0, 0, 0}, 0.5, 1.0);  b.v = V3{-1, 0, 0}; b.restitution = eb;
        w.bodies = {a, b};
        for (int s = 0; s < 400; ++s) w.step();
        return w.bodies[1].v.x - w.bodies[0].v.x;   // >0 means receding
    };
    check(close(headon(1.0, 1.0, World::Combine::Max), 2.0, 1e-3),
          "material: both e=1 -> elastic, separates at the approach speed");
    check(close(headon(0.0, 1.0, World::Combine::Max), 2.0, 1e-3),
          "material: restitution combine Max picks the bouncier body");
    check(std::fabs(headon(0.0, 1.0, World::Combine::Min)) < 0.1,
          "material: restitution combine Min -> inelastic (no separation)");
}

void test_damping() {
    World w;
    w.dt = 1.0 / 120.0; w.box.half = 1000.0;
    Body a = makeSphere(0, V3{0, 0, 0}, 1.0, 1.0);  a.v = V3{10, 0, 0}; a.linearDamping = 1.0;
    Body b = makeSphere(1, V3{0, 50, 0}, 1.0, 1.0); b.v = V3{10, 0, 0};   // undamped, far away
    w.bodies = {a, b};
    for (int s = 0; s < 120; ++s) w.step();          // 1 second
    check(w.bodies[0].v.x < 0.5 * w.bodies[1].v.x, "damping: damped body loses speed");
    check(close(w.bodies[1].v.x, 10.0, 1e-9), "damping: zero damping leaves velocity unchanged");
}

void test_collision_layers() {
    // two spheres approaching head-on; filtered onto non-matching layers they
    // must pass through, otherwise (default masks) they must collide.
    auto finalVx = [](bool filter) {
        World w;
        w.dt = 1.0 / 240.0; w.restitution = 1.0; w.box.half = 100.0;
        Body a = makeSphere(0, V3{-2, 0, 0}, 0.5, 1.0); a.v = V3{2, 0, 0};
        Body b = makeSphere(1, V3{2, 0, 0}, 0.5, 1.0);  b.v = V3{-2, 0, 0};
        if (filter) { a.layer = 0; a.mask = 1u << 0; b.layer = 1; b.mask = 1u << 1; }
        w.bodies = {a, b};
        for (int s = 0; s < 300; ++s) w.step();
        return w.bodies[0].v.x;
    };
    check(close(finalVx(true), 2.0, 1e-9), "layers: non-matching layers pass through");
    check(finalVx(false) < 0.0, "layers: default masks collide (velocity reverses)");
}

void test_userdata() {
    PhysicsWorld pw;
    pw.setBox(100.0, false); pw.setTimestep(1.0 / 240.0); pw.setRestitution(1.0);
    BodyId a = pw.addSphere(V3{-2, 0, 0}, 0.5, 1.0);
    BodyId b = pw.addSphere(V3{2, 0, 0}, 0.5, 1.0);
    pw.setVelocity(a, V3{2, 0, 0}); pw.setVelocity(b, V3{-2, 0, 0});
    pw.setUserData(a, 111); pw.setUserData(b, 222);
    bool fired = false; std::uint64_t ua = 0, ub = 0;
    pw.setContactCallback([&](const ContactInfo& ci) { fired = true; ua = ci.aUser; ub = ci.bUser; });
    for (int s = 0; s < 300; ++s) pw.step();
    check(fired, "userData: contact event fired");
    check((ua == 111 && ub == 222) || (ua == 222 && ub == 111),
          "userData: contact reports each body's userData");
}

void test_kinematic() {
    World w;
    w.dt = 1.0 / 240.0; w.box.half = 100.0; w.restitution = 0.0;
    // kinematic sphere sweeping +x into a resting dynamic sphere ahead of it.
    Body k = makeSphere(0, V3{-2, 0, 0}, 1.0, 1.0);
    k.invMass = 0; k.invInertiaBody = {0, 0, 0}; k.dynamic = false; k.kinematic = true;
    k.v = V3{2, 0, 0};
    Body d = makeSphere(1, V3{0.5, 0, 0}, 0.5, 1.0);   // dynamic, at rest, ahead
    w.bodies = {k, d};
    for (int s = 0; s < 240; ++s) w.step();             // 1 second
    check(close(w.bodies[0].x.x, 0.0, 1e-6), "kinematic: pose advances by its own velocity");
    check(close(w.bodies[0].v.x, 2.0, 1e-9), "kinematic: unaffected by the collision it causes");
    check(w.bodies[1].v.x > 0.5, "kinematic: imparts velocity to the dynamic body");
    check(w.bodies[1].x.x > 1.0, "kinematic: pushes the dynamic body along");
}

void test_contact_events() {
    // (a) resting ball: one begin when it lands, many stays, no end.
    {
        PhysicsWorld pw;
        pw.setBox(100.0, false); pw.setTimestep(1.0 / 240.0);
        pw.setGravity(V3{0, -10, 0}); pw.setRestitution(0.0); pw.setFriction(0.5);
        BodyId floor = pw.addBox(V3{0, -1, 0}, Q{1, 0, 0, 0}, V3{10, 1, 10}, 1.0);
        pw.makeStatic(floor);
        pw.addSphere(V3{0, 3, 0}, 0.5, 1.0);
        int begins = 0, stays = 0, ends = 0;
        pw.setContactBeginCallback([&](const ContactInfo&) { ++begins; });
        pw.setContactStayCallback([&](const ContactInfo&) { ++stays; });
        pw.setContactEndCallback([&](const ContactInfo&) { ++ends; });
        for (int s = 0; s < 500; ++s) pw.step();
        check(begins == 1, "contact events: begin fires once when the ball lands");
        check(stays > 10, "contact events: stay fires while resting");
        check(ends == 0, "contact events: no end while still resting");
    }
    // (b) a clean bounce fires exactly one begin and one end.
    {
        PhysicsWorld pw;
        pw.setBox(100.0, false); pw.setTimestep(1.0 / 240.0); pw.setRestitution(1.0);
        BodyId a = pw.addSphere(V3{-2, 0, 0}, 0.5, 1.0);
        BodyId b = pw.addSphere(V3{2, 0, 0}, 0.5, 1.0);
        pw.setVelocity(a, V3{2, 0, 0}); pw.setVelocity(b, V3{-2, 0, 0});
        int begins = 0, ends = 0;
        pw.setContactBeginCallback([&](const ContactInfo&) { ++begins; });
        pw.setContactEndCallback([&](const ContactInfo&) { ++ends; });
        for (int s = 0; s < 400; ++s) pw.step();
        check(begins == 1, "contact events: begin once on the bounce");
        check(ends == 1, "contact events: end once after separation");
    }
}

void test_triggers() {
    PhysicsWorld pw;
    pw.setBox(100.0, false); pw.setTimestep(1.0 / 240.0);
    BodyId trig = pw.addBox(V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
    pw.makeStatic(trig); pw.setSensor(trig, true);
    BodyId ball = pw.addSphere(V3{-5, 0, 0}, 0.5, 1.0);
    pw.setVelocity(ball, V3{5, 0, 0});
    int enter = 0, stay = 0, exit = 0;
    pw.setTriggerEnterCallback([&](const ContactInfo&) { ++enter; });
    pw.setTriggerStayCallback([&](const ContactInfo&) { ++stay; });
    pw.setTriggerExitCallback([&](const ContactInfo&) { ++exit; });
    for (int s = 0; s < 400; ++s) pw.step();
    check(enter == 1, "trigger: enter fires once");
    check(exit == 1, "trigger: exit fires once");
    check(stay >= 1, "trigger: stay fires while inside");
    check(close(pw.velocity(ball).x, 5.0, 1e-9), "trigger: no physical response (passes through)");
    check(pw.position(ball).x > 1.5, "trigger: ball continued past the sensor");
}

void test_capsule() {
    // (a) capsule-sphere: a sphere resting on the cylindrical flank of a
    //     horizontal capsule ends up ~ (capR + sphR) above the axis.
    {
        World w;
        w.dt = 1.0 / 240.0; w.box.half = 100.0; w.gravity = {0, -10, 0};
        w.restitution = 0.0; w.friction = 0.5;
        // capsule lying along X (rotate body +Y to world +X: 90 deg about Z)
        double s = std::sin(PI / 4), cth = std::cos(PI / 4);
        Body cap = makeCapsule(0, V3{0, 0, 0}, Q{cth, 0, 0, s}, 0.5, 1.5, 1.0);
        cap.invMass = 0; cap.invInertiaBody = {0, 0, 0}; cap.dynamic = false;
        Body ball = makeSphere(1, V3{0, 3, 0}, 0.5, 1.0);
        w.bodies = {cap, ball};
        for (int st = 0; st < 1500; ++st) w.step();
        check(close(w.bodies[1].x.y, 1.0, 0.05), "capsule: sphere rests on the flank at capR+sphR");
        check(std::fabs(w.bodies[1].x.x) < 0.6 && std::fabs(w.bodies[1].x.z) < 0.6,
              "capsule: sphere stays over the capsule");
    }
    // (b) capsule-capsule: two crossed capsules resolve to a real separation.
    {
        World w;
        w.dt = 1.0 / 240.0; w.box.half = 100.0; w.restitution = 0.0;
        double s = std::sin(PI / 4), cth = std::cos(PI / 4);
        Body A = makeCapsule(0, V3{0, 0, 0}, Q{cth, 0, 0, s}, 0.5, 2.0, 1.0);   // along X
        Body B = makeCapsule(1, V3{0, 0.6, 0}, Q{cth, s, 0, 0}, 0.5, 2.0, 1.0); // along Z, overlapping
        w.bodies = {A, B};
        double gap0 = w.bodies[1].x.y - w.bodies[0].x.y;
        for (int st = 0; st < 600; ++st) w.step();
        double gap1 = w.bodies[1].x.y - w.bodies[0].x.y;
        check(gap1 > gap0 + 0.2 && gap1 > 0.9, "capsule: crossed capsules push apart to ~2r");
    }
    // (c) capsule vs box through GJK/EPA (the capsule support function): a
    //     vertical capsule overlapping a static box is pushed out to touching.
    {
        World w;
        w.dt = 1.0 / 240.0; w.box.half = 100.0; w.restitution = 0.0;
        Body bx = makeBox(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
        bx.invMass = 0; bx.invInertiaBody = {0, 0, 0}; bx.dynamic = false;
        Body cap = makeCapsule(1, V3{0, 1.2, 0}, Q{1, 0, 0, 0}, 0.5, 1.0, 1.0);  // vertical, overlapping
        w.bodies = {bx, cap};
        for (int st = 0; st < 500; ++st) w.step();
        check(w.bodies[1].x.y > 2.0, "capsule vs box: EPA pushes the capsule out to touching");
    }
}

void test_queries() {
    PhysicsWorld pw;
    pw.setBox(1000.0, false); pw.setTimestep(1.0 / 240.0);
    BodyId s1 = pw.addSphere(V3{10, 0, 0}, 1.0, 1.0); pw.makeStatic(s1);
    BodyId bx = pw.addBox(V3{0, -5, 0}, Q{1, 0, 0, 0}, V3{2, 1, 2}, 1.0); pw.makeStatic(bx);
    (void)bx;

    RayHit h = pw.raycast(V3{0, 0, 0}, V3{1, 0, 0}, 100.0);
    check(h.hit && h.body == s1, "raycast: hits the sphere ahead");
    check(close(h.distance, 9.0, 1e-6), "raycast: distance to the sphere surface");
    check(close(h.normal.x, -1.0, 1e-6), "raycast: outward normal points back along the ray");

    check(!pw.raycast(V3{0, 0, 0}, V3{-1, 0, 0}, 100.0).hit, "raycast: misses when pointing away");

    pw.setLayerMask(s1, 3, 0xFFFFFFFFu);
    check(!pw.raycast(V3{0, 0, 0}, V3{1, 0, 0}, 100.0, ~(1u << 3)).hit,
          "raycast: layer mask excludes the sphere");
    pw.setLayerMask(s1, 0, 0xFFFFFFFFu);

    auto ov = pw.overlapSphere(V3{10.5, 0, 0}, 1.5);
    check(ov.size() == 1 && ov[0] == s1, "overlapSphere: finds the nearby sphere");
    check(pw.overlapSphere(V3{100, 100, 100}, 1.0).empty(), "overlapSphere: empty far away");

    RayHit sc = pw.sphereCast(V3{0, 0, 0}, 1.0, V3{1, 0, 0}, 100.0);
    check(sc.hit && sc.body == s1, "sphereCast: sweep hits the sphere");
    check(close(sc.distance, 8.0, 1e-6), "sphereCast: TOI accounts for the swept radius");
}

void test_general_joints() {
    // (a) Fixed joint: a body held rigidly beside a static anchor stays put.
    {
        PhysicsWorld pw;
        pw.setBox(1000.0, false); pw.setTimestep(1.0 / 240.0); pw.setGravity(V3{0, -10, 0});
        BodyId anchor = pw.addBox(V3{0, 5, 0}, Q{1, 0, 0, 0}, V3{0.5, 0.5, 0.5}, 1.0);
        pw.makeStatic(anchor);
        BodyId body = pw.addBox(V3{1, 5, 0}, Q{1, 0, 0, 0}, V3{0.5, 0.5, 0.5}, 1.0);
        pw.addFixedJoint(anchor, body);
        for (int s = 0; s < 600; ++s) pw.step();
        V3 p = pw.position(body);
        check(close(p.x, 1.0, 0.08) && close(p.y, 5.0, 0.15) && close(p.z, 0.0, 0.08),
              "fixed joint: holds the relative pose under gravity");
    }
    // (b) Ball joint: a bob hangs below a static anchor, pinned at the anchor.
    {
        PhysicsWorld pw;
        pw.setBox(1000.0, false); pw.setTimestep(1.0 / 240.0); pw.setGravity(V3{0, -10, 0});
        BodyId anchor = pw.addSphere(V3{0, 5, 0}, 0.2, 1.0); pw.makeStatic(anchor);
        BodyId bob = pw.addSphere(V3{0, 3, 0}, 0.5, 1.0);
        pw.addBallJoint(anchor, bob, V3{0, 5, 0});
        for (int s = 0; s < 800; ++s) pw.step();
        V3 p = pw.position(bob);
        check(close((p - V3{0, 5, 0}).norm(), 2.0, 0.05), "ball joint: anchor distance preserved");
        check(p.y < 5.0, "ball joint: hangs below the anchor");
    }
    // (c) Hinge with motor: drives rotation about its axis only.
    {
        PhysicsWorld pw;
        pw.setBox(1000.0, false); pw.setTimestep(1.0 / 240.0);
        BodyId anchor = pw.addBox(V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{0.2, 0.2, 0.2}, 1.0);
        pw.makeStatic(anchor);
        BodyId arm = pw.addBox(V3{1, 0, 0}, Q{1, 0, 0, 0}, V3{1, 0.1, 0.1}, 1.0);
        PhysicsWorld::JointId h = pw.addHingeJoint(anchor, arm, V3{0, 0, 0}, V3{0, 1, 0});
        pw.setJointMotor(h, 2.0, 200.0);
        for (int s = 0; s < 240; ++s) pw.step();
        check(std::fabs(pw.angularVelocity(arm).y) > 1.0, "hinge motor: drives rotation about the axis");
        check(std::fabs(pw.position(arm).z) > 0.2, "hinge motor: arm swept around the axis");
    }
    // (d) Breakable: enough load severs the joint and the body falls.
    {
        PhysicsWorld pw;
        pw.setBox(1000.0, false); pw.setTimestep(1.0 / 240.0); pw.setGravity(V3{0, -30, 0});
        BodyId anchor = pw.addSphere(V3{0, 5, 0}, 0.2, 1.0); pw.makeStatic(anchor);
        BodyId heavy = pw.addBox(V3{0, 3, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 50.0);
        PhysicsWorld::JointId j = pw.addBallJoint(anchor, heavy, V3{0, 5, 0});
        pw.setJointBreakable(j, 5.0);
        for (int s = 0; s < 300; ++s) pw.step();
        check(pw.jointBroken(j), "breakable joint: breaks under heavy load");
        check(pw.position(heavy).y < 2.0, "breakable joint: body falls after breaking");
    }
    // (e) Slider: free to translate along the axis, locked perpendicular.
    {
        PhysicsWorld pw;
        pw.setBox(1000.0, false); pw.setTimestep(1.0 / 240.0);
        BodyId anchor = pw.addBox(V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{0.2, 0.2, 0.2}, 1.0);
        pw.makeStatic(anchor);
        BodyId slider = pw.addBox(V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{0.5, 0.5, 0.5}, 1.0);
        pw.addSliderJoint(anchor, slider, V3{0, 0, 0}, V3{1, 0, 0});   // slide along X
        pw.setVelocity(slider, V3{2, 3, 1});                          // pushed every way
        for (int s = 0; s < 240; ++s) pw.step();
        V3 p = pw.position(slider);
        check(std::fabs(p.x) > 0.5, "slider: translates along the axis");
        check(std::fabs(p.y) < 0.05 && std::fabs(p.z) < 0.05, "slider: locked perpendicular to the axis");
    }
}

void test_plane_and_convex() {
    // (a) Plane: a sphere dropped onto an infinite ground plane rests on it.
    {
        World w;
        w.dt = 1.0 / 240.0; w.box.half = 1000.0; w.gravity = {0, -10, 0};
        w.restitution = 0.0; w.friction = 0.5;
        w.bodies.push_back(makePlane(0, V3{0, 0, 0}, V3{0, 1, 0}));
        w.bodies.push_back(makeSphere(1, V3{0, 5, 0}, 0.5, 1.0));
        for (int s = 0; s < 1200; ++s) w.step();
        check(close(w.bodies[1].x.y, 0.5, 0.03), "plane: sphere rests on the plane at its radius");
        check(w.bodies[1].v.norm() < 0.1, "plane: sphere comes to rest");
    }
    // (b) Convex hull (octahedron) dropped onto a plane settles above it.
    {
        World w;
        w.dt = 1.0 / 240.0; w.box.half = 1000.0; w.gravity = {0, -10, 0};
        w.restitution = 0.0; w.friction = 0.6;
        w.bodies.push_back(makePlane(0, V3{0, 0, 0}, V3{0, 1, 0}));
        std::vector<V3> oct = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
        w.bodies.push_back(makeConvex(1, V3{0, 4, 0}, Q{1, 0, 0, 0}, oct, 1.0));
        for (int s = 0; s < 1500; ++s) w.step();
        double y = w.bodies[1].x.y;
        check(y > 0.3 && y < 1.2, "convex: octahedron settles resting on the plane");
        check(w.bodies[1].v.norm() < 0.2, "convex: comes to rest");
    }
    // (c) Convex vs box via GJK/EPA: an overlapping convex is pushed out.
    {
        World w;
        w.dt = 1.0 / 240.0; w.box.half = 1000.0; w.restitution = 0.0;
        Body bx = makeBox(0, V3{0, 0, 0}, Q{1, 0, 0, 0}, V3{1, 1, 1}, 1.0);
        bx.invMass = 0; bx.invInertiaBody = {0, 0, 0}; bx.dynamic = false;
        std::vector<V3> tet = {{0.8, 0.8, 0.8}, {-0.8, -0.8, 0.8}, {-0.8, 0.8, -0.8}, {0.8, -0.8, -0.8}};
        w.bodies.push_back(bx);
        w.bodies.push_back(makeConvex(1, V3{0, 1.2, 0}, Q{1, 0, 0, 0}, tet, 1.0));
        double y0 = w.bodies[1].x.y;
        for (int s = 0; s < 400; ++s) w.step();
        check(w.bodies[1].x.y > y0 + 0.05, "convex vs box: EPA separates the overlap");
    }
}

void test_character() {
    PhysicsWorld pw;
    pw.setBox(1000.0, false); pw.setTimestep(1.0 / 240.0);
    BodyId floor = pw.addBox(V3{0, -1, 0}, Q{1, 0, 0, 0}, V3{20, 1, 20}, 1.0); pw.makeStatic(floor);
    BodyId wall = pw.addBox(V3{5, 2, 0}, Q{1, 0, 0, 0}, V3{0.5, 3, 10}, 1.0); pw.makeStatic(wall);
    (void)floor; (void)wall;

    CharacterController cc(&pw, 0.5, 0.5);      // radius 0.5, half-height 0.5
    cc.setPosition(V3{0, 1.02, 0});             // standing just above the floor
    cc.move(V3{0, 0, 0});
    check(cc.grounded(), "character: grounded on the floor");

    for (int i = 0; i < 200; ++i) cc.move(V3{0.1, 0, 0});   // walk +x into the wall
    V3 p = cc.position();
    check(p.x < 4.5, "character: stopped by the wall (no tunneling)");
    check(p.x > 3.0, "character: reached the wall");
    check(close(p.y, 1.02, 0.1), "character: stays at floor height while walking");
    check(cc.grounded(), "character: still grounded after walking");
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
    test_periodic_translation_invariance();
    test_cross_boundary_collision();
    test_offcentre_impulse();
    test_sphere_cylinder_geometry();
    test_sphere_cylinder_dynamics();
    test_gjk_overlap();
    test_epa_contact();
    test_cylinder_cylinder_dynamics();
    test_box();
    test_friction();
    test_gravity();
    test_physics_world();
    test_resting();
    test_stacking();
    test_sleeping();
    test_broadphase();
    test_joints();
    test_per_body_material();
    test_damping();
    test_collision_layers();
    test_userdata();
    test_kinematic();
    test_contact_events();
    test_triggers();
    test_capsule();
    test_queries();
    test_general_joints();
    test_plane_and_convex();
    test_character();
    std::printf("\n=== Summary ===\n  Passed: %d\n  Failed: %d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}

}  // namespace ne
