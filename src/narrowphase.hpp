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

// --- capsule helpers --------------------------------------------------------
// A capsule is a segment (endpoints c.x +- halfHeight * axis) swept by a sphere
// of `radius`. Contacts reduce to the closest point(s) between the segment(s).
inline double npClamp(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }

inline V3 closestOnSegment(const V3& p, const V3& a, const V3& b) {
    V3 ab = b - a;
    double L2 = ab.norm2();
    if (L2 < 1e-18) return a;
    return a + ab * npClamp((p - a).dot(ab) / L2, 0.0, 1.0);
}

// Closest points c1 (on p1q1) and c2 (on p2q2) between two segments (Ericson).
inline void closestSegSeg(const V3& p1, const V3& q1, const V3& p2, const V3& q2,
                          V3& c1, V3& c2) {
    V3 d1 = q1 - p1, d2 = q2 - p2, r = p1 - p2;
    double a = d1.dot(d1), e = d2.dot(d2), f = d2.dot(r);
    const double EPS = 1e-12;
    double s, t;
    if (a <= EPS && e <= EPS) { s = t = 0.0; }
    else if (a <= EPS) { s = 0.0; t = npClamp(f / e, 0.0, 1.0); }
    else {
        double c = d1.dot(r);
        if (e <= EPS) { t = 0.0; s = npClamp(-c / a, 0.0, 1.0); }
        else {
            double b = d1.dot(d2), denom = a * e - b * b;
            s = denom > EPS ? npClamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
            t = (b * s + f) / e;
            if (t < 0.0) { t = 0.0; s = npClamp(-c / a, 0.0, 1.0); }
            else if (t > 1.0) { t = 1.0; s = npClamp((b - c) / a, 0.0, 1.0); }
        }
    }
    c1 = p1 + d1 * s; c2 = p2 + d2 * t;
}

// Capsule (A) vs sphere (B). normal points from B toward A (detect convention),
// point is on the capsule surface.
inline Contact capsuleVsSphere(const Body& cap, const Body& sph, const Box& box) {
    Contact c;
    V3 sc = cap.x + box.minImage(sph.x - cap.x);          // sphere centre, capsule's image
    V3 axis = cap.q.rotate(V3{0, 1, 0});
    V3 cp = closestOnSegment(sc, cap.x - axis * cap.halfHeight, cap.x + axis * cap.halfHeight);
    V3 e = sc - cp;                                        // capsule -> sphere
    double dist = e.norm();
    double sumR = cap.radius + sph.radius;
    if (dist >= sumR) return c;
    V3 nrm = (dist > 1e-12) ? e / dist : V3{0, 1, 0};      // capsule -> sphere
    c.hit = true;
    c.overlap = sumR - dist;
    c.normal = -nrm;                                       // sphere -> capsule (B -> A)
    c.point = cp + nrm * cap.radius;                       // on the capsule surface
    return c;
}

// Capsule (A) vs capsule (B). normal from B toward A, point on A's surface.
inline Contact capsuleVsCapsule(const Body& A, const Body& B, const Box& box) {
    Contact c;
    V3 axA = A.q.rotate(V3{0, 1, 0}), axB = B.q.rotate(V3{0, 1, 0});
    V3 bc = A.x + box.minImage(B.x - A.x);                 // B centre, A's image
    V3 pA, pB;
    closestSegSeg(A.x - axA * A.halfHeight, A.x + axA * A.halfHeight,
                  bc - axB * B.halfHeight, bc + axB * B.halfHeight, pA, pB);
    V3 e = pB - pA;                                        // A -> B
    double dist = e.norm();
    double sumR = A.radius + B.radius;
    if (dist >= sumR) return c;
    V3 nrm = (dist > 1e-12) ? e / dist : V3{0, 1, 0};      // A -> B
    c.hit = true;
    c.overlap = sumR - dist;
    c.normal = -nrm;                                       // B -> A
    c.point = pA + nrm * A.radius;                         // on A's surface
    return c;
}

}  // namespace ne

#endif  // NATIVEENGINE_NARROWPHASE_HPP
