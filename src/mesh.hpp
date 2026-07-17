#pragma once
#ifndef NATIVEENGINE_MESH_HPP
#define NATIVEENGINE_MESH_HPP

#include <algorithm>
#include <vector>

#include "vmath.hpp"

namespace ne {

// Static triangle mesh with an AABB BVH (mid-phase) so a moving body only tests
// the few triangles near it, not the whole level. Built once; never mutated.
struct MeshData {
    struct Tri { int a, b, c; };
    struct Node {
        V3 lo, hi;
        int left = -1, right = -1;   // internal: child indices; leaf: both -1
        int start = 0, count = 0;    // leaf: [start,start+count) into `order`
    };

    std::vector<V3> verts;
    std::vector<Tri> tris;
    std::vector<int> order;          // triangle indices, BVH-sorted
    std::vector<Node> nodes;

    V3 v(int tri, int k) const {
        const Tri& t = tris[tri];
        return verts[k == 0 ? t.a : (k == 1 ? t.b : t.c)];
    }

    void triBounds(int tri, V3& lo, V3& hi) const {
        V3 a = v(tri, 0), b = v(tri, 1), c = v(tri, 2);
        lo = {std::min({a.x, b.x, c.x}), std::min({a.y, b.y, c.y}), std::min({a.z, b.z, c.z})};
        hi = {std::max({a.x, b.x, c.x}), std::max({a.y, b.y, c.y}), std::max({a.z, b.z, c.z})};
    }

    void build() {
        order.resize(tris.size());
        for (int i = 0; i < (int)tris.size(); ++i) order[i] = i;
        nodes.clear();
        if (!tris.empty()) buildNode(0, (int)tris.size());
    }

    // Collect triangles whose AABB overlaps [lo,hi] into `out`.
    void query(const V3& lo, const V3& hi, std::vector<int>& out) const {
        if (nodes.empty()) return;
        int stack[64]; int sp = 0; stack[sp++] = 0;
        while (sp > 0) {
            const Node& nd = nodes[stack[--sp]];
            if (lo.x > nd.hi.x || hi.x < nd.lo.x || lo.y > nd.hi.y || hi.y < nd.lo.y ||
                lo.z > nd.hi.z || hi.z < nd.lo.z) continue;
            if (nd.left < 0) {
                for (int i = 0; i < nd.count; ++i) out.push_back(order[nd.start + i]);
            } else {
                if (sp < 62) { stack[sp++] = nd.left; stack[sp++] = nd.right; }
            }
        }
    }

private:
    int buildNode(int start, int count) {
        int idx = (int)nodes.size();
        nodes.push_back({});
        V3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
        for (int i = 0; i < count; ++i) {
            V3 tlo, thi; triBounds(order[start + i], tlo, thi);
            lo = {std::min(lo.x, tlo.x), std::min(lo.y, tlo.y), std::min(lo.z, tlo.z)};
            hi = {std::max(hi.x, thi.x), std::max(hi.y, thi.y), std::max(hi.z, thi.z)};
        }
        Node nd; nd.lo = lo; nd.hi = hi;
        if (count <= 4) {                       // leaf
            nd.start = start; nd.count = count;
            nodes[idx] = nd;
            return idx;
        }
        V3 ext = hi - lo;
        int axis = (ext.x > ext.y && ext.x > ext.z) ? 0 : (ext.y > ext.z ? 1 : 2);
        int mid = start + count / 2;
        std::nth_element(order.begin() + start, order.begin() + mid, order.begin() + start + count,
                         [&](int A, int B) {
                             V3 al, ah, bl, bh; triBounds(A, al, ah); triBounds(B, bl, bh);
                             double ca = (axis == 0 ? al.x + ah.x : axis == 1 ? al.y + ah.y : al.z + ah.z);
                             double cb = (axis == 0 ? bl.x + bh.x : axis == 1 ? bl.y + bh.y : bl.z + bh.z);
                             return ca < cb;
                         });
        int l = buildNode(start, mid - start);
        int r = buildNode(mid, start + count - mid);
        nd.left = l; nd.right = r;
        nodes[idx] = nd;
        return idx;
    }
};

}  // namespace ne

#endif  // NATIVEENGINE_MESH_HPP
