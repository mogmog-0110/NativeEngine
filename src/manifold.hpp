#pragma once
#ifndef NATIVEENGINE_MANIFOLD_HPP
#define NATIVEENGINE_MANIFOLD_HPP

#include <vector>

#include "body.hpp"
#include "box.hpp"
#include "gjk.hpp"   // minImageShift

namespace ne {

// A multi-point contact manifold. Single-point contact makes a box rest on one
// corner and tip; a resting box-on-box face contact needs the ~4-point patch
// this produces, so boxes settle flat and stack.
struct Manifold {
    int count = 0;
    V3 point[4];
    double depth[4];
};

namespace detail {

inline V3 axisVec(int a) {
    return a == 0 ? V3{1, 0, 0} : (a == 1 ? V3{0, 1, 0} : V3{0, 0, 1});
}
inline double comp(const V3& v, int a) { return a == 0 ? v.x : (a == 1 ? v.y : v.z); }

struct BoxFace {
    V3 n;              // outward world normal
    V3 c;              // face centre (world)
    V3 u, v;           // world tangent axes (unit)
    double eu, ev;     // half-lengths along u, v
    V3 vert[4];        // world corners
};

// Face of a box (given its orientation/extents and a world centre) along body
// axis `axis` with the given sign.
inline BoxFace boxFace(const Body& b, const V3& centreWorld, int axis, int sign) {
    const int t1 = (axis + 1) % 3, t2 = (axis + 2) % 3;
    const double h[3] = {b.halfExtents.x, b.halfExtents.y, b.halfExtents.z};
    BoxFace f;
    f.n = b.q.rotate(axisVec(axis) * double(sign));
    f.c = centreWorld + f.n * h[axis];
    f.u = b.q.rotate(axisVec(t1));  f.eu = h[t1];
    f.v = b.q.rotate(axisVec(t2));  f.ev = h[t2];
    f.vert[0] = f.c + f.u * f.eu + f.v * f.ev;
    f.vert[1] = f.c - f.u * f.eu + f.v * f.ev;
    f.vert[2] = f.c - f.u * f.eu - f.v * f.ev;
    f.vert[3] = f.c + f.u * f.eu - f.v * f.ev;
    return f;
}

// The face whose outward normal is most aligned with `dir`.
inline BoxFace bestFace(const Body& b, const V3& centreWorld, const V3& dir) {
    int bestAxis = 0, bestSign = 1;
    double best = -1e30;
    for (int axis = 0; axis < 3; ++axis)
        for (int sign = -1; sign <= 1; sign += 2) {
            double a = b.q.rotate(axisVec(axis) * double(sign)).dot(dir);
            if (a > best) { best = a; bestAxis = axis; bestSign = sign; }
        }
    return boxFace(b, centreWorld, bestAxis, bestSign);
}

// Sutherland-Hodgman: keep the part of `poly` with (p . pn) <= off.
inline void clipPlane(std::vector<V3>& poly, const V3& pn, double off) {
    std::vector<V3> out;
    const int m = (int)poly.size();
    for (int i = 0; i < m; ++i) {
        const V3& cur = poly[i];
        const V3& nxt = poly[(i + 1) % m];
        double dc = cur.dot(pn) - off, dn = nxt.dot(pn) - off;
        bool inC = dc <= 0.0, inN = dn <= 0.0;
        if (inC) out.push_back(cur);
        if (inC != inN) {
            double t = dc / (dc - dn);
            out.push_back(cur + (nxt - cur) * t);
        }
    }
    poly.swap(out);
}

}  // namespace detail

// Box-box contact manifold via reference/incident face clipping. `n` is the
// collision normal (world, b -> a) from EPA. Returns up to 4 points on the
// contact patch with their penetration depths; count 0 => fall back to the
// single EPA point (edge/corner contact).
inline Manifold boxBoxManifold(const Body& a, const Body& b, const V3& n, const Box& box) {
    using namespace detail;
    Manifold m;
    const V3 bShift = minImageShift(a, b, box);        // b in a's nearest image
    const V3 bC = b.x + bShift;

    // a faces b along -n; b faces a along +n. Reference = the flatter alignment.
    BoxFace fa = bestFace(a, a.x, -n);
    BoxFace fb = bestFace(b, bC, n);
    double alignA = fa.n.dot(-n), alignB = fb.n.dot(n);

    const BoxFace& ref = (alignA >= alignB) ? fa : fb;
    const Body& incBody = (alignA >= alignB) ? b : a;
    const V3 incCentre = (alignA >= alignB) ? bC : a.x;

    // Incident face = most anti-parallel to the reference normal.
    BoxFace inc = bestFace(incBody, incCentre, -ref.n);

    // Clip the incident polygon against the reference face's four side planes.
    std::vector<V3> poly{inc.vert[0], inc.vert[1], inc.vert[2], inc.vert[3]};
    clipPlane(poly, ref.u,  ref.c.dot(ref.u) + ref.eu);
    clipPlane(poly, ref.u * -1.0, -(ref.c.dot(ref.u) - ref.eu));
    clipPlane(poly, ref.v,  ref.c.dot(ref.v) + ref.ev);
    clipPlane(poly, ref.v * -1.0, -(ref.c.dot(ref.v) - ref.ev));

    // Keep points on/behind the reference face (penetrating); project onto it.
    const double slop = 0.01;
    for (const V3& p : poly) {
        double sd = (p - ref.c).dot(ref.n);            // signed distance to face
        if (sd > slop) continue;
        if (m.count >= 4) break;
        m.point[m.count] = p - ref.n * sd;             // on the reference face
        m.depth[m.count] = -sd;                        // >= -slop
        ++m.count;
    }
    return m;
}

}  // namespace ne

#endif  // NATIVEENGINE_MANIFOLD_HPP
