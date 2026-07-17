#pragma once
#ifndef NATIVEENGINE_BODY_HPP
#define NATIVEENGINE_BODY_HPP

#include <algorithm>
#include <cstdint>
#include <vector>

#include "vmath.hpp"

namespace ne {

// A rigid body. Sphere and (later) cylinder share this; the shape is a tag plus
// the parameters the narrow phase needs. Inertia is stored as the DIAGONAL of
// the inverse body-frame inertia tensor, valid because both shapes are integrated
// about their principal axes.
//   Convex : an arbitrary convex hull given by `vertices` (body frame).
//   Plane  : an infinite static half-space; `halfExtents` holds the world normal
//            and the plane passes through `x`.
enum class Shape { Sphere, Cylinder, Box, Capsule, Convex, Plane };

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
    double radius = 1.0;       // sphere radius, or cylinder/capsule radius
    double halfHeight = 0.0;   // cylinder/capsule half-height (0 for sphere)
    V3 halfExtents;            // box half-extents (Box); world normal (Plane)
    std::vector<V3> vertices;  // convex-hull vertices, body frame (Convex only)

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

    // Per-body gravity multiplier (1 = full world gravity; 0 = floats; <0 = rises).
    double gravityScale = 1.0;

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

    // Sensor / trigger: detected for overlap events but produces NO physical
    // response (the solver skips any pair touching a sensor). A trigger volume.
    bool sensor = false;

    // Continuous collision: when set, a fast body is swept from its pre-step pose
    // and stopped at the first impact so it cannot tunnel through thin geometry in
    // one step. Off by default (the discrete solver is used). prevX is the pose
    // captured before integration, for the sweep.
    bool ccd = false;
    V3 prevX;
    Q prevQ{1, 0, 0, 0};   // pre-step orientation, for render interpolation

    bool isSphere() const { return shape == Shape::Sphere; }

    // Radius of a bounding sphere about the centre (for broadphase).
    double boundingRadius() const {
        if (shape == Shape::Sphere) return radius;
        if (shape == Shape::Box) return halfExtents.norm();
        if (shape == Shape::Capsule) return halfHeight + radius;   // segment + cap
        if (shape == Shape::Plane) return 1e30;                    // infinite (grid-excluded)
        if (shape == Shape::Convex) {
            double m = 0.0;
            for (const V3& vert : vertices) m = std::max(m, vert.norm());
            return m;
        }
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

    // World-frame inverse inertia as a full 3x3 matrix: R * diag(invIbody) * R^T.
    // Used by the joint solver (needs the matrix, not just its action on a vector).
    Mat3 invInertiaWorld() const {
        Mat3 R = rotationMatrix();
        return R * Mat3::diagonal(invInertiaBody) * R.transposed();
    }
    // Body-to-world rotation as a 3x3 (columns are the rotated basis vectors).
    Mat3 rotationMatrix() const {
        V3 cx = q.rotate({1, 0, 0}), cy = q.rotate({0, 1, 0}), cz = q.rotate({0, 0, 1});
        return Mat3{cx.x, cy.x, cz.x, cx.y, cy.y, cz.y, cx.z, cy.z, cz.z};
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

// Factory: a solid capsule -- a cylinder of radius r and length 2*halfHeight
// (axis along body +Y) capped by two hemispheres. The inertia is exact in both
// limits: a sphere when halfHeight -> 0, a thin rod when r -> 0.
inline Body makeCapsule(int id, const V3& x, const Q& q, double radius,
                        double halfHeight, double density) {
    Body b;
    b.shape = Shape::Capsule;
    b.id = id;
    b.x = x;
    b.q = q;
    b.radius = radius;
    b.halfHeight = halfHeight;
    const double r = radius, L = 2.0 * halfHeight;      // cylindrical length
    const double mc = density * ne::PI * r * r * L;     // cylinder part
    const double mh = density * (4.0 / 3.0) * ne::PI * r * r * r;  // two hemispheres
    const double m = mc + mh;
    b.invMass = (m > 0) ? 1.0 / m : 0.0;
    // About the axis (Y): cylinder 1/2 m r^2 + hemispheres 2/5 m r^2.
    double Iyy = 0.5 * mc * r * r + 0.4 * mh * r * r;
    // Perpendicular: cylinder m(r^2/4 + L^2/12) + caps offset ~ L/2.
    double Ixx = mc * (0.25 * r * r + L * L / 12.0) + mh * (0.4 * r * r + 0.25 * L * L);
    b.invInertiaBody = {Ixx > 0 ? 1.0 / Ixx : 0.0,
                        Iyy > 0 ? 1.0 / Iyy : 0.0,
                        Ixx > 0 ? 1.0 / Ixx : 0.0};
    return b;
}

// Factory: an arbitrary convex hull from body-frame vertices. Mass and inertia
// are approximated by the vertices' axis-aligned bounding box (adequate for
// gameplay); the exact convex path is the same GJK/EPA machinery via the support.
inline Body makeConvex(int id, const V3& x, const Q& q, std::vector<V3> verts,
                       double density) {
    Body b;
    b.shape = Shape::Convex;
    b.id = id;
    b.x = x;
    b.q = q;
    b.vertices = std::move(verts);
    V3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
    for (const V3& v : b.vertices) {
        lo = V3{std::min(lo.x, v.x), std::min(lo.y, v.y), std::min(lo.z, v.z)};
        hi = V3{std::max(hi.x, v.x), std::max(hi.y, v.y), std::max(hi.z, v.z)};
    }
    V3 he = (hi - lo) * 0.5;
    const double hx = std::max(he.x, 1e-4), hy = std::max(he.y, 1e-4), hz = std::max(he.z, 1e-4);
    double m = density * 8.0 * hx * hy * hz;
    b.invMass = (m > 0) ? 1.0 / m : 0.0;
    double Ix = (1.0 / 3.0) * m * (hy * hy + hz * hz);
    double Iy = (1.0 / 3.0) * m * (hx * hx + hz * hz);
    double Iz = (1.0 / 3.0) * m * (hx * hx + hy * hy);
    b.invInertiaBody = {Ix > 0 ? 1.0 / Ix : 0.0, Iy > 0 ? 1.0 / Iy : 0.0, Iz > 0 ? 1.0 / Iz : 0.0};
    return b;
}

// Factory: an infinite static plane (half-space) through `point` with outward
// `normal`. Always static: nothing can move it.
inline Body makePlane(int id, const V3& point, const V3& normal) {
    Body b;
    b.shape = Shape::Plane;
    b.id = id;
    b.x = point;
    b.halfExtents = normal.normalized();   // world normal stored here
    b.invMass = 0.0;
    b.invInertiaBody = {0, 0, 0};
    b.dynamic = false;
    return b;
}

}  // namespace ne

#endif  // NATIVEENGINE_BODY_HPP
