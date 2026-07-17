#include "world.hpp"

#include "contact.hpp"
#include "narrowphase.hpp"

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
        // Symplectic Euler: velocities first, then positions.
        b.v += b.force * (b.invMass * dt);
        b.w += b.applyInvInertiaWorld(b.torque) * dt;
        b.x += b.v * dt;
        b.q = integrateOrientation(b.q, b.w, dt);
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

            if (a.isSphere() && b.isSphere()) {
                V3 d = box.minImage(a.x - b.x);           // b -> a, nearest image
                double sumR = a.radius + b.radius;
                double dist2 = d.norm2();
                if (dist2 >= sumR * sumR || dist2 < 1e-18) continue;
                double dist = std::sqrt(dist2);
                V3 normal = d / dist;                      // b -> a
                // Contact on the line of centres: lever arms along +-normal, so
                // no spin, but through the same general resolver.
                resolveContact(a, b, normal * (-a.radius), normal * (b.radius),
                               normal, restitution);
                correctPenetration(a, b, normal, sumR - dist, contactBeta);
            } else if (a.isSphere() != b.isSphere()) {
                // One sphere, one cylinder. Resolve as (sphere = A, cylinder = B)
                // so the contact normal points from the cylinder toward the sphere.
                Body& sph = a.isSphere() ? a : b;
                Body& cyl = a.isSphere() ? b : a;
                Contact c = sphereVsCylinder(sph, cyl, box);
                if (!c.hit) continue;
                V3 rSph = c.point - sph.x;
                V3 rCyl = box.minImage(c.point - cyl.x);
                resolveContact(sph, cyl, rSph, rCyl, c.normal, restitution);
                // Positional correction along the (off-centre) contact normal.
                double invSum = sph.invMass + cyl.invMass;
                if (invSum > 0.0) {
                    V3 corr = c.normal * (contactBeta * c.overlap / invSum);
                    sph.x += corr * sph.invMass;
                    cyl.x -= corr * cyl.invMass;
                }
            }
            // (cylinder-cylinder handled in the GJK/EPA increment)
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
