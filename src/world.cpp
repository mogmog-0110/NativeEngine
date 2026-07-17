#include "world.hpp"

#include "contact.hpp"
#include "narrowphase.hpp"
#include "epa.hpp"
#include "detect.hpp"

#include <vector>

namespace ne {

// A one-point contact constraint for the iterative solver. Single point per pair
// for now (box-box face manifolds via clipping come next); the accumulated
// impulses + restitution bias make it a proper sequential-impulse solver that
// handles several simultaneous contacts on one body.
namespace {
struct Constraint {
    std::size_t i, j;
    V3 n;              // b -> a
    V3 rA, rB;         // lever arms (minimum-imaged)
    double overlap;
    double restBias;   // target separating speed from restitution
    double nImp = 0.0, tImp = 0.0;
};

// Below this approach speed, restitution is suppressed so resting contacts do
// not jitter or gain energy.
constexpr double kRestitutionSlop = 0.5;
}  // namespace

void World::step() {
    integrate();
    collide();
    if (box.periodic)
        wrapPositions();
    else
        applyWalls();
}

void World::integrate() {
    for (Body& b : bodies) {
        if (b.invMass <= 0.0) {          // infinite mass: never moves
            b.force = {};
            b.torque = {};
            continue;
        }
        // Linear: symplectic Euler. Gravity is an acceleration (mass-independent);
        // applied forces divide by mass.
        b.v += (gravity + b.force * b.invMass) * dt;
        b.x += b.v * dt;

        // Rotational: integrate ANGULAR MOMENTUM, not angular velocity. For an
        // asymmetric body (the cylinder, Ixx != Iyy) a torque-free w is NOT
        // constant -- only L is -- so advancing w directly would leak angular
        // momentum. Advance L by the torque, rotate by the current w, then
        // recompute w = I_world(q_new)^-1 L so L is conserved exactly at zero
        // torque. (For the isotropic sphere this reduces to constant w.)
        V3 L = b.angularMomentum();      // I_world(q) w, at the current pose
        L += b.torque * dt;
        b.q = integrateOrientation(b.q, b.w, dt);
        b.w = b.applyInvInertiaWorld(L); // w consistent with L at the new pose
        b.force = {};
        b.torque = {};
    }
}

// All-pairs hard contact. O(N^2) with minimum-image separation; a cell list
// replaces this once correctness is established. Iterating i<j in index order
// keeps the impulse application deterministic. Cylinder-cylinder narrow phase
// (GJK/EPA) is added in the next increment; here we handle sphere-sphere and
// sphere-cylinder.
void World::collide() {
    const size_t n = bodies.size();

    // 1. Detect: build a constraint per overlapping pair. Deterministic i<j order.
    std::vector<Constraint> cons;
    for (size_t i = 0; i < n; ++i) {
        Body& a = bodies[i];
        for (size_t j = i + 1; j < n; ++j) {
            Body& b = bodies[j];
            if (a.invMass + b.invMass <= 0.0) continue;
            Contact c = detectContact(a, b, box);
            if (!c.hit) continue;
            Constraint k;
            k.i = i; k.j = j; k.n = c.normal; k.overlap = c.overlap;
            k.rA = box.minImage(c.point - a.x);
            k.rB = box.minImage(c.point - b.x);
            // Restitution bias from the approach speed at contact.
            V3 vc = (a.v + a.w.cross(k.rA)) - (b.v + b.w.cross(k.rB));
            double vn = vc.dot(k.n);
            k.restBias = (vn < -kRestitutionSlop) ? -restitution * vn : 0.0;
            cons.push_back(k);
        }
    }

    // 2. Velocity solve: sequential impulse with accumulated clamping. Several
    // iterations so simultaneous contacts on one body settle consistently.
    const int kIters = 8;
    for (int it = 0; it < kIters; ++it) {
        for (Constraint& k : cons) {
            Body& a = bodies[k.i];
            Body& b = bodies[k.j];
            // Normal impulse toward the target separating speed (restBias).
            V3 vc = (a.v + a.w.cross(k.rA)) - (b.v + b.w.cross(k.rB));
            double vn = vc.dot(k.n);
            double Kn = effMass(a, b, k.rA, k.rB, k.n);
            if (Kn > 1e-18) {
                double dj = (k.restBias - vn) / Kn;
                double old = k.nImp;
                k.nImp = (old + dj > 0.0) ? old + dj : 0.0;   // clamp >= 0
                dj = k.nImp - old;
                V3 J = k.n * dj;
                a.v += J * a.invMass; a.w += a.applyInvInertiaWorld(k.rA.cross(J));
                b.v -= J * b.invMass; b.w -= b.applyInvInertiaWorld(k.rB.cross(J));
            }
            // Friction: oppose tangential velocity, clamped to the Coulomb cone.
            if (friction > 0.0 && k.nImp > 0.0) {
                V3 vc2 = (a.v + a.w.cross(k.rA)) - (b.v + b.w.cross(k.rB));
                V3 vt = vc2 - k.n * vc2.dot(k.n);
                double vtl = vt.norm();
                if (vtl > 1e-12) {
                    V3 t = vt / vtl;
                    double Kt = effMass(a, b, k.rA, k.rB, t);
                    if (Kt > 1e-18) {
                        double djt = -vtl / Kt;
                        double maxF = friction * k.nImp;
                        double old = k.tImp;
                        double sum = old + djt;
                        k.tImp = (sum < -maxF) ? -maxF : (sum > maxF ? maxF : sum);
                        djt = k.tImp - old;
                        V3 Jt = t * djt;
                        a.v += Jt * a.invMass; a.w += a.applyInvInertiaWorld(k.rA.cross(Jt));
                        b.v -= Jt * b.invMass; b.w -= b.applyInvInertiaWorld(k.rB.cross(Jt));
                    }
                }
            }
        }
    }

    // 3. Positional correction (geometric, energy-neutral).
    for (const Constraint& k : cons)
        correctPenetration(bodies[k.i], bodies[k.j], k.n, k.overlap, contactBeta);
}

// Reflective walls at +-half. A sphere whose surface crosses a wall has that
// velocity component reflected (scaled by restitution) and its position projected
// back inside. Energy-conserving at e = 1.
void World::applyWalls() {
    const double h = box.half;
    for (Body& b : bodies) {
        if (b.invMass <= 0.0 || !b.isSphere()) continue;
        double r = b.radius;
        double* px[3] = {&b.x.x, &b.x.y, &b.x.z};
        double* pv[3] = {&b.v.x, &b.v.y, &b.v.z};
        for (int k = 0; k < 3; ++k) {
            double lim = h - r;
            if (*px[k] > lim) {
                *px[k] = 2.0 * lim - *px[k];
                if (*pv[k] > 0) *pv[k] = -restitution * *pv[k];
            } else if (*px[k] < -lim) {
                *px[k] = -2.0 * lim - *px[k];
                if (*pv[k] < 0) *pv[k] = -restitution * *pv[k];
            }
        }
    }
}

void World::wrapPositions() {
    for (Body& b : bodies)
        if (b.invMass > 0.0)
            b.x = box.wrap(b.x);
}

V3 World::totalLinearMomentum() const {
    V3 p;
    for (const Body& b : bodies)
        if (b.invMass > 0.0) p += b.v * (1.0 / b.invMass);
    return p;
}

V3 World::totalAngularMomentum() const {
    // L_origin = sum( r x m v + I_world w ). Reflective only -- under PBC the
    // origin-referenced L is not conserved (wrapping shifts r by a lattice vector).
    V3 L;
    for (const Body& b : bodies) {
        if (b.invMass <= 0.0) continue;
        double m = 1.0 / b.invMass;
        L += b.x.cross(b.v * m);
        L += b.angularMomentum();
    }
    return L;
}

double World::totalKinetic() const {
    double e = 0.0;
    for (const Body& b : bodies)
        if (b.invMass > 0.0) e += b.kinetic();
    return e;
}

}  // namespace ne
