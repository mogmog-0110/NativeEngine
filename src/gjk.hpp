#pragma once
#ifndef NATIVEENGINE_GJK_HPP
#define NATIVEENGINE_GJK_HPP

#include <array>

#include "body.hpp"
#include "box.hpp"

namespace ne {

// Support mapping: the farthest point of a body's surface along a world
// direction d. This is all GJK/EPA need about a shape, so a sphere and a
// flat-cap cylinder are handled by the same convex machinery -- and both support
// functions are trivial. (A smooth flat-cap cylinder is what "cylindrical linker"
// means; PhysX faceted it into a 16-gon, but the exact convex is cleaner here.)
inline V3 supportWorld(const Body& b, const V3& d) {
    if (b.isSphere()) return b.x + d.normalized() * b.radius;
    // Cylinder: axis along body +Y. Support = farther cap plane along y, and the
    // rim in the radial direction.
    V3 db = b.q.inverseRotate(d);
    double sy = (db.y >= 0.0) ? b.halfHeight : -b.halfHeight;
    double rho = std::sqrt(db.x * db.x + db.z * db.z);
    V3 sb = (rho > 1e-12) ? V3{db.x / rho * b.radius, sy, db.z / rho * b.radius}
                          : V3{0.0, sy, 0.0};
    return b.x + b.q.rotate(sb);
}

// Support of the Minkowski difference A (-) B, with B shifted by `bShift` so the
// pair is evaluated in the correct periodic image.
inline V3 csoSupport(const Body& A, const Body& B, const V3& bShift, const V3& d) {
    return supportWorld(A, d) - (supportWorld(B, -d) + bShift);
}

// Boolean GJK: does the Minkowski difference contain the origin, i.e. do A and B
// overlap? Standard evolving-simplex GJK (Muratori formulation). `bShift` places
// B at its minimum-image position relative to A.
inline bool gjkOverlap(const Body& A, const Body& B, const V3& bShift) {
    V3 d{1, 0, 0};
    V3 a = csoSupport(A, B, bShift, d);
    if (a.dot(d) <= 0.0) {
        // try another direction before giving up (d may be tangent)
        d = V3{0, 1, 0};
        a = csoSupport(A, B, bShift, d);
        if (a.dot(d) <= 0.0) return false;
    }
    std::array<V3, 4> simplex{a, {}, {}, {}};
    int n = 1;
    d = -a;

    for (int iter = 0; iter < 64; ++iter) {
        V3 p = csoSupport(A, B, bShift, d);
        if (p.dot(d) <= 1e-12) return false;      // no origin past this support
        // prepend p
        for (int i = n; i > 0; --i) simplex[i] = simplex[i - 1];
        simplex[0] = p;
        ++n;

        // doSimplex: update simplex + search direction; return true if origin
        // enclosed.
        V3 A0 = simplex[0];
        V3 ao = -A0;
        if (n == 2) {
            V3 B0 = simplex[1];
            V3 ab = B0 - A0;
            if (ab.dot(ao) > 0.0) {
                d = ab.cross(ao).cross(ab);
                if (d.norm2() < 1e-24) { d = V3{ab.z, ab.x, ab.y}; }  // ao ∥ ab
            } else {
                simplex[0] = A0; n = 1; d = ao;
            }
        } else if (n == 3) {
            V3 B0 = simplex[1], C0 = simplex[2];
            V3 ab = B0 - A0, ac = C0 - A0;
            V3 abc = ab.cross(ac);
            if (abc.cross(ac).dot(ao) > 0.0) {
                if (ac.dot(ao) > 0.0) { simplex[1] = C0; n = 2; d = ac.cross(ao).cross(ac); }
                else { n = 2; d = ab.dot(ao) > 0 ? ab.cross(ao).cross(ab) : ao; if (ab.dot(ao) <= 0) n = 1; }
            } else if (ab.cross(abc).dot(ao) > 0.0) {
                n = 2; d = ab.dot(ao) > 0 ? ab.cross(ao).cross(ab) : ao; if (ab.dot(ao) <= 0) n = 1;
            } else {
                if (abc.dot(ao) > 0.0) { d = abc; }
                else { V3 t = simplex[1]; simplex[1] = simplex[2]; simplex[2] = t; d = -abc; }
            }
        } else {  // n == 4: tetrahedron -- test the three new faces
            V3 B0 = simplex[1], C0 = simplex[2], D0 = simplex[3];
            V3 ab = B0 - A0, ac = C0 - A0, ad = D0 - A0;
            V3 abc = ab.cross(ac), acd = ac.cross(ad), adb = ad.cross(ab);
            bool over = true;
            if (abc.dot(ao) > 0.0) { simplex[3] = C0; simplex[2] = B0; simplex[1] = A0; n = 3; d = abc; over = false; }
            else if (acd.dot(ao) > 0.0) { simplex[1] = C0; simplex[2] = D0; n = 3; d = acd; over = false; }
            else if (adb.dot(ao) > 0.0) { simplex[2] = D0; simplex[1] = B0; n = 3; d = adb; over = false; }
            if (over) return true;   // origin inside the tetrahedron
        }
        if (d.norm2() < 1e-24) return true;  // degenerate: treat as touching
    }
    return false;
}

// Convenience: min-image shift that places B in the image nearest A.
inline V3 minImageShift(const Body& A, const Body& B, const Box& box) {
    V3 raw = B.x - A.x;
    return box.minImage(raw) - raw;
}

}  // namespace ne

#endif  // NATIVEENGINE_GJK_HPP
