#pragma once
#ifndef NATIVEENGINE_BOX_HPP
#define NATIVEENGINE_BOX_HPP

#include "vmath.hpp"

namespace ne {

// Cubic domain of half-edge `half` (so the interior is [-half, +half]^3 and the
// full edge L = 2*half, matching the PhysX project's roomSize convention).
//
// This is where "native PBC" actually lives: when periodic, pair separations use
// the minimum-image convention and positions wrap through the faces. There are no
// ghost bodies -- the toroidal topology is in the distance function itself, which
// is exactly what an off-the-shelf rigid-body engine cannot express.
struct Box {
    double half = 12.0;
    bool periodic = false;

    double edge() const { return 2.0 * half; }

    // Minimum-image displacement d, wrapped into (-half, +half] per axis when
    // periodic. Reflective boxes return d unchanged.
    V3 minImage(V3 d) const {
        if (!periodic) return d;
        const double L = edge();
        d.x -= L * std::floor(d.x / L + 0.5);
        d.y -= L * std::floor(d.y / L + 0.5);
        d.z -= L * std::floor(d.z / L + 0.5);
        return d;
    }

    // Wrap a position back into the primary cell [-half, +half). Periodic only.
    V3 wrap(V3 p) const {
        if (!periodic) return p;
        const double L = edge();
        p.x -= L * std::floor((p.x + half) / L);
        p.y -= L * std::floor((p.y + half) / L);
        p.z -= L * std::floor((p.z + half) / L);
        return p;
    }
};

}  // namespace ne

#endif  // NATIVEENGINE_BOX_HPP
