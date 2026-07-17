#pragma once
#ifndef NATIVEENGINE_BODY_HPP
#define NATIVEENGINE_BODY_HPP

#include <cstdint>

#include "vmath.hpp"

namespace ne {

// A rigid body. Sphere and (later) cylinder share this; the shape is a tag plus
// the parameters the narrow phase needs. Inertia is stored as the DIAGONAL of
// the inverse body-frame inertia tensor, valid because both shapes are integrated
// about their principal axes.
enum class Shape { Sphere, Cylinder, Box };

struct Body {
    // State
    V3 x;         // centre of mass position
    V3 v;         // linear velocity
    Q  q{1, 0, 0, 0};  // orientation (body -> world)
    V3 w;         // angular velocity (world frame)

    // Accumulators (cleared each step)
    V3 force;
    V3 torque;

    // Mass / inertia
    double invMass = 0.0;      // 0 => infinite mass (kinematic / wall)
    V3 invInertiaBody;         // diagonal of inverse inertia in body frame

    // Shape
    Shape shape = Shape::Sphere;
    double radius = 1.0;       // sphere radius, or cylinder radius
    double halfHeight = 0.0;   // cylinder half-height (0 for sphere)
    V3 halfExtents;            // box half-extents (Box only)

    int id = -1;               // stable id (sphere ids and cylinder ids are
                               // separate namespaces, matching the science layer)

    // Per-body material. -1 means "inherit the World's global friction /
    // restitution", so bodies that never set a material behave exactly as before.
    // At a contact the two bodies' effective values are combined (World combine
    // modes); combine(x,x)=x, so an all-inherit scene is unchanged.
    double friction = -1.0;
    double restitution = -1.0;

    // Per-body velocity damping (drag). 0 = none. Applied each step:
    // v *= 1/(1 + linearDamping*dt), w *= 1/(1 + angularDamping*dt).
    double linearDamping = 0.0;
    double angularDamping = 0.0;

    // Collision filtering. `layer` is an index (0..31); `mask` is the set of
    // layers this body collides with. Two bodies collide iff each one's mask
    // contains the other's layer. Defaults (layer 0, mask = all) collide with
    // everything, so unfiltered scenes are unchanged.
    std::uint32_t layer = 0;
    std::uint32_t mask = 0xFFFFFFFFu;

    // Opaque handle for the embedder to map a body back to its gameplay entity;
    // surfaced in contact/query results. The engine never interprets it.
    std::uint64_t userData = 0;

    // Sleeping: a settled dynamic body stops integrating (and acts as immovable
    // in the solver) until disturbed. Implemented by zeroing invMass/invInertia
    // while asleep -- the integrator and solver already treat 0 as immovable --
    // and restoring them on wake. `dynamic` distinguishes a genuine static body
    // (never sleeps/wakes) from a sleeping dynamic one (both have invMass 0).
    bool dynamic = true;
    bool sleeping = false;
    double sleepTimer = 0.0;
    double invMassStore = 0.0;
    V3 invInertiaStore;

    // Kinematic: infinite mass (invMass=0, so the solver never pushes it) BUT a
    // real velocity that DOES enter contact resolution, and a pose advanced from
    // that velocity each step. This is a moving platform / elevator / animated
    // collider: it drives dynamic bodies but is itself driven only by the user.
    // Distinct from a static body (v=0, pose fixed) and never sleeps.
    bool kinematic = false;

    bool isSphere() const { return shape == Shape::Sphere; }

    // Radius of a bounding sphere about the centre (for broadphase).
    double boundingRadius() const {
        if (shape == Shape::Sphere) return radius;
        if (shape == Shape::Box) return halfExtents.norm();
        return std::sqrt(halfHeight * halfHeight + radius * radius);  // cylinder
    }

    // World-frame inverse inertia applied to a world vector:  Iinv_world * u
    // = R * diag(invIbody) * R^T * u.
    V3 applyInvInertiaWorld(const V3& u) const {
        V3 ub = q.inverseRotate(u);                       // world -> body
        V3 sb{ub.x * invInertiaBody.x,
              ub.y * invInertiaBody.y,
              ub.z * invInertiaBody.z};                   // scale by diagonal
        return q.rotate(sb);                              // body -> world
    }

    double kinetic() const {
        double lin = 0.5 * (invMass > 0 ? 1.0 / invMass : 0.0) * v.norm2();
        // rotational: 0.5 * w . (I_world w) = 0.5 * w . L
        V3 Ibody{invInertiaBody.x > 0 ? 1.0 / invInertiaBody.x : 0.0,
                 invInertiaBody.y > 0 ? 1.0 / invInertiaBody.y : 0.0,
                 invInertiaBody.z > 0 ? 1.0 / invInertiaBody.z : 0.0};
        V3 wb = q.inverseRotate(w);
        V3 Lb{wb.x * Ibody.x, wb.y * Ibody.y, wb.z * Ibody.z};
        double rot = 0.5 * wb.dot(Lb);
        return lin + rot;
    }

    // World-frame angular momentum L = I_world w.
    V3 angularMomentum() const {
        V3 Ibody{invInertiaBody.x > 0 ? 1.0 / invInertiaBody.x : 0.0,
                 invInertiaBody.y > 0 ? 1.0 / invInertiaBody.y : 0.0,
                 invInertiaBody.z > 0 ? 1.0 / invInertiaBody.z : 0.0};
        V3 wb = q.inverseRotate(w);
        V3 Lb{wb.x * Ibody.x, wb.y * Ibody.y, wb.z * Ibody.z};
        return q.rotate(Lb);
    }
};

// Factory: a solid sphere of unit density scaled by `density`.
inline Body makeSphere(int id, const V3& x, double radius, double density) {
    Body b;
    b.shape = Shape::Sphere;
    b.id = id;
    b.x = x;
    b.radius = radius;
    b.halfHeight = 0.0;
    double m = density * (4.0 / 3.0) * ne::PI * radius * radius * radius;
    b.invMass = (m > 0) ? 1.0 / m : 0.0;
    double I = (2.0 / 5.0) * m * radius * radius;      // solid sphere
    double invI = (I > 0) ? 1.0 / I : 0.0;
    b.invInertiaBody = {invI, invI, invI};
    return b;
}

// Factory: a solid cylinder, axis along body +Y (matching the PhysX convention).
inline Body makeCylinder(int id, const V3& x, const Q& q, double radius,
                         double height, double density) {
    Body b;
    b.shape = Shape::Cylinder;
    b.id = id;
    b.x = x;
    b.q = q;
    b.radius = radius;
    b.halfHeight = 0.5 * height;
    double m = density * ne::PI * radius * radius * height;
    b.invMass = (m > 0) ? 1.0 / m : 0.0;
    double Iyy = 0.5 * m * radius * radius;                              // about axis
    double Ixx = (1.0 / 12.0) * m * (3.0 * radius * radius + height * height);
    b.invInertiaBody = {Ixx > 0 ? 1.0 / Ixx : 0.0,
                        Iyy > 0 ? 1.0 / Iyy : 0.0,
                        Ixx > 0 ? 1.0 / Ixx : 0.0};
    return b;
}

// Factory: a solid box with the given half-extents (games' bread-and-butter
// shape). Handled by the same GJK/EPA convex machinery as the cylinder.
inline Body makeBox(int id, const V3& x, const Q& q, const V3& halfExtents,
                    double density) {
    Body b;
    b.shape = Shape::Box;
    b.id = id;
    b.x = x;
    b.q = q;
    b.halfExtents = halfExtents;
    const double hx = halfExtents.x, hy = halfExtents.y, hz = halfExtents.z;
    double m = density * 8.0 * hx * hy * hz;
    b.invMass = (m > 0) ? 1.0 / m : 0.0;
    double Ix = (1.0 / 3.0) * m * (hy * hy + hz * hz);   // (1/12) m ((2hy)^2+(2hz)^2)
    double Iy = (1.0 / 3.0) * m * (hx * hx + hz * hz);
    double Iz = (1.0 / 3.0) * m * (hx * hx + hy * hy);
    b.invInertiaBody = {Ix > 0 ? 1.0 / Ix : 0.0,
                        Iy > 0 ? 1.0 / Iy : 0.0,
                        Iz > 0 ? 1.0 / Iz : 0.0};
    return b;
}

}  // namespace ne

#endif  // NATIVEENGINE_BODY_HPP
