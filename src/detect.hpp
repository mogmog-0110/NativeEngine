#pragma once
#ifndef NATIVEENGINE_DETECT_HPP
#define NATIVEENGINE_DETECT_HPP

#include "body.hpp"
#include "box.hpp"
#include "narrowphase.hpp"
#include "epa.hpp"

namespace ne {

// Unified narrow-phase dispatch: returns a Contact whose `normal` points from b
// toward a (the resolver / manifold convention), `point` is on a's side, and
// `overlap` is the penetration depth. Fast exact analytic paths for
// sphere-sphere and sphere-cylinder; GJK+EPA for everything else.
// Infinite half-space (plane) vs any convex body. `plane.halfExtents` is the
// world outward normal; the plane passes through plane.x. Returns normal pointing
// from the plane toward `other` (plane -> other), point at the deepest vertex.
inline Contact planeVsConvex(const Body& other, const Body& plane, const Box& box) {
    Contact c;
    V3 n = plane.halfExtents;                          // unit world normal
    V3 sp = supportWorld(other, -n);                   // deepest point toward the plane
    V3 rel = box.minImage(sp - plane.x);
    double pen = -n.dot(rel);                          // >0 means below the surface
    if (pen <= 0.0) return c;
    c.hit = true;
    c.overlap = pen;
    c.normal = n;                                      // plane -> other
    c.point = sp + n * (0.5 * pen);                    // midway to the surface
    return c;
}

inline Contact detectContact(const Body& a, const Body& b, const Box& box) {
    Contact c;
    // Plane (infinite half-space) vs anything: analytic, not GJK.
    if (a.shape == Shape::Plane || b.shape == Shape::Plane) {
        if (a.shape == Shape::Plane && b.shape == Shape::Plane) return c;   // two planes: ignore
        const Body& plane = (a.shape == Shape::Plane) ? a : b;
        const Body& other = (a.shape == Shape::Plane) ? b : a;
        Contact t = planeVsConvex(other, plane, box);   // normal plane -> other
        if (t.hit && a.shape == Shape::Plane) t.normal = -t.normal;   // want b -> a
        return t;
    }
    if (a.isSphere() && b.isSphere()) {
        V3 d = box.minImage(a.x - b.x);
        double sumR = a.radius + b.radius, dist2 = d.norm2();
        if (dist2 >= sumR * sumR || dist2 < 1e-18) return c;
        double dist = std::sqrt(dist2);
        c.hit = true;
        c.normal = d / dist;                 // b -> a
        c.overlap = sumR - dist;
        c.point = a.x - c.normal * a.radius;  // on a's surface
        return c;
    }
    const bool aCyl = a.shape == Shape::Cylinder;
    const bool bCyl = b.shape == Shape::Cylinder;
    if (a.isSphere() && bCyl) return sphereVsCylinder(a, b, box);   // normal b->a
    if (b.isSphere() && aCyl) {                                     // normal was a->b
        Contact t = sphereVsCylinder(b, a, box);
        if (t.hit) t.normal = -t.normal;                            // -> b->a
        return t;
    }

    // Capsule analytic paths (smooth shape -> cleaner than EPA). capsule-box and
    // capsule-cylinder fall through to GJK/EPA via the capsule support function.
    const bool aCap = a.shape == Shape::Capsule;
    const bool bCap = b.shape == Shape::Capsule;
    if (aCap && bCap) return capsuleVsCapsule(a, b, box);           // normal b->a
    if (aCap && b.isSphere()) return capsuleVsSphere(a, b, box);    // A=cap, normal b->a
    if (a.isSphere() && bCap) {
        Contact t = capsuleVsSphere(b, a, box);                     // normal sph->cap
        if (t.hit) t.normal = -t.normal;                            // -> b->a
        return t;
    }

    return convexContact(a, b, box);          // normal b->a
}

}  // namespace ne

#endif  // NATIVEENGINE_DETECT_HPP
