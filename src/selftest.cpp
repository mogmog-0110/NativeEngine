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
    std::printf("\n=== Summary ===\n  Passed: %d\n  Failed: %d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}

}  // namespace ne
