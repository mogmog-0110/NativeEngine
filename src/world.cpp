#include "world.hpp"

namespace ne {

void World::step() {
    integrate();
    collideSpheres();
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

// Sphere-sphere hard contact. O(N^2) with minimum-image separation; a cell list
// replaces this once correctness is established. Iterating i<j in index order
// keeps the impulse application deterministic.
void World::collideSpheres() {
    const size_t n = bodies.size();
    for (size_t i = 0; i < n; ++i) {
        Body& a = bodies[i];
        if (!a.isSphere()) continue;
        for (size_t j = i + 1; j < n; ++j) {
            Body& b = bodies[j];
            if (!b.isSphere()) continue;

            V3 d = box.minImage(a.x - b.x);           // b -> a, nearest image
            double sumR = a.radius + b.radius;
            double dist2 = d.norm2();
            if (dist2 >= sumR * sumR || dist2 < 1e-18) continue;

            double dist = std::sqrt(dist2);
            V3 normal = d / dist;                       // unit normal, b -> a
            double vrel = (a.v - b.v).dot(normal);      // closing if < 0

            double invSum = a.invMass + b.invMass;
            if (invSum <= 0.0) continue;

            if (vrel < 0.0) {
                double jimp = -(1.0 + restitution) * vrel / invSum;
                a.v += normal * (jimp * a.invMass);
                b.v -= normal * (jimp * b.invMass);
            }
            // Positional correction: split the overlap by inverse-mass share so
            // the pair stops interpenetrating without injecting velocity.
            double overlap = sumR - dist;
            V3 corr = normal * (contactBeta * overlap / invSum);
            a.x += corr * a.invMass;
            b.x -= corr * b.invMass;
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
