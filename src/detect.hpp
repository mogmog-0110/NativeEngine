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
inline Contact detectContact(const Body& a, const Body& b, const Box& box) {
    Contact c;
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
    return convexContact(a, b, box);          // normal b->a
}

}  // namespace ne

#endif  // NATIVEENGINE_DETECT_HPP
