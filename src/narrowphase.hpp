#pragma once
#ifndef NATIVEENGINE_NARROWPHASE_HPP
#define NATIVEENGINE_NARROWPHASE_HPP

#include "body.hpp"
#include "box.hpp"

namespace ne {

// A resolved contact between two bodies.
//   normal  : unit, points from body B toward body A (the resolveContact
//             convention). For sphere-vs-cylinder we call with A = sphere,
//             B = cylinder, so `normal` is the OUTWARD cylinder surface normal.
//   point   : the contact point in world space (on the cylinder surface).
//   overlap : positive penetration depth (surface-to-surface).
struct Contact {
    bool hit = false;
    V3 normal;
    V3 point;
    double overlap = 0.0;
};

// Sphere vs finite FLAT-CAP cylinder (axis along the cylinder's body +Y).
//
// Exact closest-point-on-a-solid-cylinder, matching the PhysX model's flat caps
// (PhysX faceted it into a 16-gon prism; a smooth flat-cap cylinder is cleaner
// and is what "cylindrical linker" means). The three exterior regions -- side,
// cap face, and cap rim -- all fall out of clamping the local coordinates:
//   y  clamped to [-h, +h]   selects cap face vs side band
//   radial clamped to r      selects side vs axis
// giving side (radial>r, |y|<h), cap (radial<=r, |y|>h) and rim (both) uniformly.
// The deep-interior case (inside both) is pushed out along the least-penetrating
// axis so the contact stays well defined.
inline Contact sphereVsCylinder(const Body& sph, const Body& cyl, const Box& box) {
    Contact c;
    const double r = cyl.radius, h = cyl.halfHeight, sr = sph.radius;

    // Sphere centre in the cylinder's body frame, minimum-imaged.
    V3 rel = box.minImage(sph.x - cyl.x);
    V3 L = cyl.q.inverseRotate(rel);

    double rho = std::sqrt(L.x * L.x + L.z * L.z);
    bool insideAxial = std::fabs(L.y) < h;
    bool insideRadial = rho < r;

    V3 Plocal;      // closest point on the cylinder surface, body frame
    V3 nLocal;      // outward surface normal at that point, body frame

    if (insideAxial && insideRadial) {
        // Sphere centre is inside the solid: exit along the nearer face.
        double dCap = h - std::fabs(L.y);     // distance to nearer cap plane
        double dSide = r - rho;               // distance to side
        if (dSide < dCap) {
            double s = (rho > 1e-12) ? r / rho : 0.0;
            Plocal = V3{L.x * s, L.y, L.z * s};
            nLocal = (rho > 1e-12) ? V3{L.x / rho, 0, L.z / rho} : V3{1, 0, 0};
            c.overlap = sr + dSide;
        } else {
            double sy = (L.y >= 0) ? h : -h;
            Plocal = V3{L.x, sy, L.z};
            nLocal = V3{0, (L.y >= 0) ? 1.0 : -1.0, 0};
            c.overlap = sr + dCap;
        }
        c.hit = true;
    } else {
        double yc = (L.y > h) ? h : (L.y < -h ? -h : L.y);
        double rx = L.x, rz = L.z;
        if (rho > r) { double s = r / rho; rx *= s; rz *= s; }
        Plocal = V3{rx, yc, rz};
        V3 diff = L - Plocal;
        double dist = diff.norm();
        if (dist >= sr) return c;             // no contact
        c.hit = true;
        c.overlap = sr - dist;
        nLocal = (dist > 1e-12) ? diff / dist : V3{0, 1, 0};
    }

    c.normal = cyl.q.rotate(nLocal).normalized();   // outward, world frame
    c.point = cyl.x + cyl.q.rotate(Plocal);          // on the cylinder surface
    return c;
}

}  // namespace ne

#endif  // NATIVEENGINE_NARROWPHASE_HPP
