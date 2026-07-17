#pragma once
#ifndef NATIVEENGINE_WORLD_HPP
#define NATIVEENGINE_WORLD_HPP

#include <vector>

#include "body.hpp"
#include "box.hpp"

namespace ne {

// The simulation world: a set of rigid bodies in a (periodic or reflective) box,
// advanced by a fixed-timestep symplectic integrator with impulse-based hard
// contact. "Ballistic between contacts" holds because forces are zero unless a
// bond spring-damper is added (later); collisions are resolved as instantaneous
// momentum-conserving impulses.
//
// Determinism: every loop iterates bodies / pairs in index order and reduces in
// that order, so a given build produces identical output for identical input.
class World {
public:
    Box box;
    double dt = 1.0 / 120.0;
    double restitution = 1.0;   // e in [0,1]; 1 = perfectly elastic
    double contactBeta = 0.2;   // positional-correction fraction per step

    std::vector<Body> bodies;

    void step();

    // Aggregate diagnostics used by the physics selftest.
    V3 totalLinearMomentum() const;
    V3 totalAngularMomentum() const;   // about the origin
    double totalKinetic() const;

private:
    void integrate();
    void collide();        // all shape pairs
    void applyWalls();     // reflective boundary
    void wrapPositions();  // periodic boundary
};

}  // namespace ne

#endif  // NATIVEENGINE_WORLD_HPP
