#pragma once
#ifndef NATIVEENGINE_QUERIES_HPP
#define NATIVEENGINE_QUERIES_HPP

#include <cmath>

#include "body.hpp"
#include "box.hpp"

// Scene queries: ray vs a single body, in world space. `d` must be unit. On a
// hit, `tOut` is the distance along the ray (<= maxT), `pOut` the world point,
// `nOut` the outward surface normal at the hit. Callers place the body at its
// minimum-image position first, so these are plain world-space tests.
namespace ne {

// Ray vs a sphere centred at c. Returns the near intersection (or the exit if the
// origin is inside).
inline bool rayVsSphereAt(const V3& o, const V3& d, const V3& c, double r,
                          double maxT, double& tOut, V3& nOut) {
    V3 m = o - c;
    double b = m.dot(d), cc = m.dot(m) - r * r;
    if (cc > 0.0 && b > 0.0) return false;         // origin outside, pointing away
    double disc = b * b - cc;
    if (disc < 0.0) return false;
    double t = -b - std::sqrt(disc);
    if (t < 0.0) t = 0.0;                           // origin inside
    if (t > maxT) return false;
    tOut = t;
    nOut = ((o + d * t) - c).normalized();
    return true;
}

namespace detail {
// Ray vs box in the box's local frame (centre at origin, half-extents he).
inline bool rayVsBoxLocal(const V3& ol, const V3& dl, const V3& he, double maxT,
                          double& tOut, V3& nLocal) {
    double tmin = 0.0, tmax = maxT;
    int axis = -1; double sign = 1.0;
    const double* o = &ol.x; const double* dd = &dl.x; const double* h = &he.x;
    for (int a = 0; a < 3; ++a) {
        if (std::fabs(dd[a]) < 1e-12) {
            if (o[a] < -h[a] || o[a] > h[a]) return false;   // parallel & outside
        } else {
            double inv = 1.0 / dd[a];
            double t1 = (-h[a] - o[a]) * inv, t2 = (h[a] - o[a]) * inv;
            double s = -1.0;
            if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; s = 1.0; }
            if (t1 > tmin) { tmin = t1; axis = a; sign = s; }
            if (t2 < tmax) tmax = t2;
            if (tmin > tmax) return false;
        }
    }
    if (axis < 0) return false;   // origin already inside on every axis
    tOut = tmin;
    nLocal = {axis == 0 ? sign : 0.0, axis == 1 ? sign : 0.0, axis == 2 ? sign : 0.0};
    return true;
}
}  // namespace detail

// Ray vs any body. Sphere/box/cylinder/capsule handled analytically.
inline bool rayVsBody(const Body& b, const V3& c, const V3& o, const V3& d,
                      double maxT, double& tOut, V3& pOut, V3& nOut) {
    if (b.shape == Shape::Sphere) {
        if (!rayVsSphereAt(o, d, c, b.radius, maxT, tOut, nOut)) return false;
        pOut = o + d * tOut; return true;
    }
    // Work in the body's local frame (its centre at c).
    V3 ol = b.q.inverseRotate(o - c);
    V3 dl = b.q.inverseRotate(d);

    if (b.shape == Shape::Box) {
        V3 nl;
        if (!detail::rayVsBoxLocal(ol, dl, b.halfExtents, maxT, tOut, nl)) return false;
        nOut = b.q.rotate(nl); pOut = o + d * tOut; return true;
    }

    const double r = b.radius, h = b.halfHeight;
    double best = maxT + 1.0; V3 nl;
    // Radial side: (ox+t dx)^2 + (oz+t dz)^2 = r^2.
    double a2 = dl.x * dl.x + dl.z * dl.z;
    double b2 = 2.0 * (ol.x * dl.x + ol.z * dl.z);
    double c2 = ol.x * ol.x + ol.z * ol.z - r * r;
    if (a2 > 1e-18) {
        double disc = b2 * b2 - 4.0 * a2 * c2;
        if (disc >= 0.0) {
            double t = (-b2 - std::sqrt(disc)) / (2.0 * a2);
            if (t < 0.0) t = (-b2 + std::sqrt(disc)) / (2.0 * a2);
            if (t >= 0.0 && t <= maxT) {
                double y = ol.y + t * dl.y;
                if (y >= -h && y <= h && t < best) {   // on the finite side band
                    best = t; nl = V3{(ol.x + t * dl.x) / r, 0.0, (ol.z + t * dl.z) / r};
                }
            }
        }
    }
    if (b.shape == Shape::Cylinder) {
        // Flat caps: planes y = +-h, inside the disc.
        for (int s = -1; s <= 1; s += 2) {
            double yc = s * h;
            if (std::fabs(dl.y) < 1e-12) continue;
            double t = (yc - ol.y) / dl.y;
            if (t < 0.0 || t > maxT || t >= best) continue;
            double px = ol.x + t * dl.x, pz = ol.z + t * dl.z;
            if (px * px + pz * pz <= r * r) { best = t; nl = V3{0.0, (double)s, 0.0}; }
        }
    } else {   // Capsule: hemispherical caps at the segment ends.
        if (best > maxT) {   // side missed (or y out of band): try the cap spheres
            for (int s = -1; s <= 1; s += 2) {
                double tS; V3 nS;
                if (rayVsSphereAt(ol, dl, V3{0.0, (double)s * h, 0.0}, r, maxT, tS, nS) && tS < best) {
                    best = tS; nl = nS;
                }
            }
        }
    }
    if (best > maxT) return false;
    tOut = best; nOut = b.q.rotate(nl).normalized(); pOut = o + d * tOut;
    return true;
}

}  // namespace ne

#endif  // NATIVEENGINE_QUERIES_HPP
