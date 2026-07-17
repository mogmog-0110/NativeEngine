#include "world.hpp"

#include "contact.hpp"
#include "narrowphase.hpp"
#include "epa.hpp"

namespace ne {

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
        // Linear: symplectic Euler.
        b.v += b.force * (b.invMass * dt);
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
    for (size_t i = 0; i < n; ++i) {
        Body& a = bodies[i];
        for (size_t j = i + 1; j < n; ++j) {
            Body& b = bodies[j];
            if (a.invMass + b.invMass <= 0.0) continue;

            const bool sphereCyl =
                (a.isSphere() && b.shape == Shape::Cylinder) ||
                (b.isSphere() && a.shape == Shape::Cylinder);

            if (a.isSphere() && b.isSphere()) {
                // Sphere-sphere: fast analytic path (contact on the line of
                // centres, no spin), through the same general resolver.
                V3 d = box.minImage(a.x - b.x);           // b -> a, nearest image
                double sumR = a.radius + b.radius;
                double dist2 = d.norm2();
                if (dist2 >= sumR * sumR || dist2 < 1e-18) continue;
                double dist = std::sqrt(dist2);
                V3 normal = d / dist;                      // b -> a
                resolveContact(a, b, normal * (-a.radius), normal * (b.radius),
                               normal, restitution, friction);
                correctPenetration(a, b, normal, sumR - dist, contactBeta);
            } else if (sphereCyl) {
                // Sphere-cylinder: exact analytic narrow phase (sphere = A so the
                // normal points from the cylinder toward the sphere).
                Body& sph = a.isSphere() ? a : b;
                Body& cyl = a.isSphere() ? b : a;
                Contact c = sphereVsCylinder(sph, cyl, box);
                if (!c.hit) continue;
                resolveContact(sph, cyl, c.point - sph.x, box.minImage(c.point - cyl.x),
                               c.normal, restitution, friction);
                double invSum = sph.invMass + cyl.invMass;
                if (invSum > 0.0) {
                    V3 corr = c.normal * (contactBeta * c.overlap / invSum);
                    sph.x += corr * sph.invMass;
                    cyl.x -= corr * cyl.invMass;
                }
            } else {
                // Everything else (box-*, cylinder-cylinder): general convex
                // contact via GJK + EPA. Normal points from B (=bodies[j]) to A.
                Contact c = convexContact(a, b, box);
                if (!c.hit) continue;
                resolveContact(a, b, box.minImage(c.point - a.x),
                               box.minImage(c.point - b.x), c.normal, restitution, friction);
                double invSum = a.invMass + b.invMass;
                if (invSum > 0.0) {
                    V3 corr = c.normal * (contactBeta * c.overlap / invSum);
                    a.x += corr * a.invMass;
                    b.x -= corr * b.invMass;
                }
            }
        }
    }
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
