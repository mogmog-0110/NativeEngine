#pragma once
#ifndef NATIVEENGINE_EPA_HPP
#define NATIVEENGINE_EPA_HPP

#include <array>
#include <vector>

#include "gjk.hpp"
#include "narrowphase.hpp"   // Contact

namespace ne {

// A Minkowski-difference vertex that remembers which surface points on A and B
// produced it, so EPA can recover a world-space contact point.
struct SupportV {
    V3 v;   // support(A,d) - support(B,-d)   (in A's image frame)
    V3 a;   // the supporting point on A
    V3 b;   // the supporting point on B (already shifted to A's nearest image)
};

inline SupportV csoSupportV(const Body& A, const Body& B, const V3& bShift, const V3& d) {
    V3 pa = supportWorld(A, d);
    V3 pb = supportWorld(B, -d) + bShift;
    return {pa - pb, pa, pb};
}

// GJK that, on overlap, returns a tetrahedron (4 SupportV) enclosing the origin.
inline bool gjkTetra(const Body& A, const Body& B, const V3& bShift,
                     std::array<SupportV, 4>& tet) {
    V3 d{1, 0, 0};
    SupportV s = csoSupportV(A, B, bShift, d);
    if (s.v.norm2() < 1e-24) { d = V3{0, 1, 0}; s = csoSupportV(A, B, bShift, d); }
    std::vector<SupportV> W{s};
    d = -s.v;

    for (int iter = 0; iter < 64; ++iter) {
        if (d.norm2() < 1e-24) d = V3{1, 0, 0};
        SupportV p = csoSupportV(A, B, bShift, d);
        if (p.v.dot(d) <= 1e-12) return false;
        W.insert(W.begin(), p);

        const V3 a = W[0].v, ao = -a;
        if (W.size() == 2) {
            V3 ab = W[1].v - a;
            d = ab.cross(ao).cross(ab);
            if (d.norm2() < 1e-24) d = V3{ab.y, -ab.x, 0};
        } else if (W.size() == 3) {
            V3 ab = W[1].v - a, ac = W[2].v - a;
            V3 abc = ab.cross(ac);
            if (abc.dot(ao) > 0) d = abc;
            else { std::swap(W[1], W[2]); d = -abc; }
            if (d.norm2() < 1e-24) return false;
        } else {  // 4: tetra -- check which face (if any) the origin is beyond
            V3 b = W[1].v, c = W[2].v, dd = W[3].v;
            V3 ab = b - a, ac = c - a, ad = dd - a;
            V3 abc = ab.cross(ac), acd = ac.cross(ad), adb = ad.cross(ab);
            if (abc.dot(ao) > 0)      { W = {W[0], W[1], W[2]}; d = abc; }
            else if (acd.dot(ao) > 0) { W = {W[0], W[2], W[3]}; d = acd; }
            else if (adb.dot(ao) > 0) { W = {W[0], W[3], W[1]}; d = adb; }
            else { tet = {W[0], W[1], W[2], W[3]}; return true; }  // origin enclosed
        }
    }
    return false;
}

// EPA: expand the tetrahedron toward the origin's nearest face on the Minkowski
// boundary, giving penetration normal + depth + a world contact point.
// Returns a Contact with `normal` pointing from B toward A (resolveContact's
// convention), `overlap` = penetration depth, `point` = midpoint of the witness
// points on A and B.
inline Contact epaContact(const Body& A, const Body& B, const V3& bShift,
                          const std::array<SupportV, 4>& tet) {
    struct Face { int a, b, c; V3 n; double dist; };
    std::vector<SupportV> verts(tet.begin(), tet.end());

    auto makeFace = [&](int i, int j, int k) -> Face {
        V3 n = (verts[j].v - verts[i].v).cross(verts[k].v - verts[i].v);
        double L = n.norm();
        n = (L > 1e-18) ? n / L : V3{0, 0, 1};
        double dist = n.dot(verts[i].v);
        if (dist < 0) { n = -n; dist = -dist; std::swap(j, k); }  // outward
        return {i, j, k, n, dist};
    };
    std::vector<Face> faces{makeFace(0, 1, 2), makeFace(0, 2, 3),
                            makeFace(0, 3, 1), makeFace(1, 3, 2)};

    Face best{};
    for (int iter = 0; iter < 64; ++iter) {
        // closest face to the origin
        int ci = 0;
        for (int i = 1; i < (int)faces.size(); ++i)
            if (faces[i].dist < faces[ci].dist) ci = i;
        best = faces[ci];

        SupportV p = csoSupportV(A, B, bShift, best.n);
        double d = p.v.dot(best.n);
        if (d - best.dist < 1e-7) break;   // converged: face is on the boundary

        // Remove faces the new point can see; collect the horizon edges.
        int pidx = (int)verts.size();
        verts.push_back(p);
        std::vector<std::array<int, 2>> horizon;
        std::vector<Face> keep;
        for (const Face& f : faces) {
            if (f.n.dot(p.v - verts[f.a].v) > 1e-12) {
                auto edge = [&](int u, int w) {
                    for (size_t e = 0; e < horizon.size(); ++e)
                        if (horizon[e][0] == w && horizon[e][1] == u) {
                            horizon.erase(horizon.begin() + e); return;
                        }
                    horizon.push_back({u, w});
                };
                edge(f.a, f.b); edge(f.b, f.c); edge(f.c, f.a);
            } else {
                keep.push_back(f);
            }
        }
        faces = keep;
        for (auto& e : horizon) faces.push_back(makeFace(e[0], e[1], pidx));
        if (faces.empty()) break;
    }

    // Barycentric projection of the origin onto the best face -> witness points.
    const SupportV& A0 = verts[best.a];
    const SupportV& B0 = verts[best.b];
    const SupportV& C0 = verts[best.c];
    V3 proj = best.n * best.dist;
    V3 v0 = B0.v - A0.v, v1 = C0.v - A0.v, v2 = proj - A0.v;
    double d00 = v0.dot(v0), d01 = v0.dot(v1), d11 = v1.dot(v1);
    double d20 = v2.dot(v0), d21 = v2.dot(v1);
    double den = d00 * d11 - d01 * d01;
    double bv = 0, bw = 0;
    if (std::fabs(den) > 1e-18) {
        bv = (d11 * d20 - d01 * d21) / den;
        bw = (d00 * d21 - d01 * d20) / den;
    }
    double bu = 1.0 - bv - bw;
    V3 onA = A0.a * bu + B0.a * bv + C0.a * bw;
    V3 onB = A0.b * bu + B0.b * bv + C0.b * bw;

    Contact c;
    c.hit = true;
    // best.n is the outward normal of the closest face of the CSO A(-)B, i.e. it
    // points from the origin toward the boundary = from A toward B. resolveContact
    // wants the normal from B toward A (so +J separates A from B), so negate.
    c.normal = -best.n;
    c.overlap = best.dist;
    c.point = (onA + onB) * 0.5;
    return c;
}

// Convex-convex contact via GJK + EPA. Returns {hit=false} if separated.
inline Contact convexContact(const Body& A, const Body& B, const Box& box) {
    V3 shift = minImageShift(A, B, box);
    std::array<SupportV, 4> tet;
    if (!gjkTetra(A, B, shift, tet)) return Contact{};
    return epaContact(A, B, shift, tet);
}

}  // namespace ne

#endif  // NATIVEENGINE_EPA_HPP
