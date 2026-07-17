#pragma once
#ifndef NATIVEENGINE_WORLD_HPP
#define NATIVEENGINE_WORLD_HPP

#include <utility>
#include <vector>

#include "body.hpp"
#include "box.hpp"

namespace ne {

// A distance joint between two bodies at local anchor points. rest is the target
// separation. stiffness == 0 makes it RIGID (a rod / rope link, solved as a
// bilateral velocity constraint with Baumgarte length correction); stiffness > 0
// makes it a spring-damper (applied as a force). Anchoring to a fixed point is
// done by joining to a static body.
struct DistanceJoint {
    std::size_t a = 0, b = 0;
    V3 localA, localB;
    double rest = 0.0;
    double stiffness = 0.0;   // 0 = rigid
    double damping = 0.0;
};

// A general bilateral joint solved by sequential impulse. All types pin the two
// local anchor points together (the ball-socket / point constraint); Hinge and
// Fixed add angular constraints on top; Slider frees one translation axis instead.
enum class JointType { Ball, Hinge, Fixed, Slider };

struct Joint {
    JointType type = JointType::Ball;
    std::size_t a = 0, b = 0;
    V3 localA, localB;               // anchor point in each body frame
    V3 axisA{0, 1, 0}, axisB{0, 1, 0};   // hinge / slider axis in each body frame
    Q  refRel{1, 0, 0, 0};           // reference qA * qB^-1 (captured at creation)

    bool useMotor = false;
    double motorSpeed = 0.0;         // target relative rate about the axis
    double maxMotorImpulse = 0.0;    // clamp

    bool useLimit = false;
    double lower = 0.0, upper = 0.0; // hinge: angle (rad); slider: distance

    bool breakable = false;
    double breakImpulse = 0.0;

    // runtime
    bool broken = false;
    double appliedImpulse = 0.0;     // accumulated this step (break test)
    double motorAcc = 0.0;           // motor impulse accumulator (per step)
};

// A resolved contact point, captured for debug visualization (opt-in). point and
// normal are in world space; normal points from body j toward body i.
struct ContactViz {
    V3 point;
    V3 normal;
    double depth = 0.0;
    std::size_t i = 0, j = 0;
};

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
    double restitution = 1.0;   // e in [0,1]; 1 = perfectly elastic (default material)
    double friction = 0.0;      // Coulomb mu; 0 = frictionless (default material)
    double contactBeta = 0.2;   // positional-correction fraction per step
    V3 gravity;                 // uniform acceleration (0 for the science runs)

    // How a pair's two per-body materials combine at a contact. Defaults chosen
    // so that combine(x,x)=x -- an all-default scene reduces to the old global
    // behaviour exactly (Multiply is intentionally NOT the default).
    enum class Combine { Average, Min, Max, Multiply, GeometricMean };
    Combine restitutionCombine = Combine::Max;
    Combine frictionCombine = Combine::GeometricMean;

    // A body's effective material: its own value, or the world default if unset
    // (sentinel < 0). Pair values are then run through the combine modes above.
    double effFriction(const Body& b) const { return b.friction < 0 ? friction : b.friction; }
    double effRestitution(const Body& b) const { return b.restitution < 0 ? restitution : b.restitution; }
    // Layer/mask filter: collide iff each body's mask contains the other's layer.
    static bool layersCollide(const Body& a, const Body& b) {
        return (a.mask & (1u << (b.layer & 31u))) && (b.mask & (1u << (a.layer & 31u)));
    }

    // Sleeping thresholds. A dynamic body below sleepLinVel/sleepAngVel for
    // sleepTime seconds sleeps; a contact with an awake partner faster than
    // wakeVel wakes it. sleepEnabled=false disables the whole mechanism (the
    // science runs keep every body awake).
    bool sleepEnabled = false;
    double sleepLinVel = 0.05, sleepAngVel = 0.05, sleepTime = 0.5, wakeVel = 0.5;

    // Open boundary: no reflective walls (objects rest on a ground body and can
    // fall off its edges) -- the normal game setup. Ignored when periodic.
    // Default false keeps the reflective-box behaviour the science runs rely on.
    bool openBoundary = false;

    // Broadphase: a PBC-aware uniform grid replaces the O(N^2) pair scan for
    // large scenes. Candidate pairs are sorted to (i<j) order so the constraint
    // list -- and thus the result -- is bit-identical to the brute-force scan.
    // forceBruteForce is a test hook to compare the two.
    bool forceBruteForce = false;

    std::vector<Body> bodies;
    std::vector<DistanceJoint> distanceJoints;
    std::vector<Joint> joints;

    // Debug: when captureContacts is set, collide() records each resolved contact
    // point into debugContacts (cleared each step). Off by default -- zero cost to
    // the solver; a tool (the viewer) turns it on to draw contacts.
    bool captureContacts = false;
    std::vector<ContactViz> debugContacts;

    void step();
    static void wake(Body& b);

    // Aggregate diagnostics used by the physics selftest.
    V3 totalLinearMomentum() const;
    V3 totalAngularMomentum() const;   // about the origin
    double totalKinetic() const;

private:
    void integrate();
    void collide();        // all shape pairs
    void applyJointForces();      // spring joints
    void solveRigidJoints(int iterations);   // bilateral distance constraints
    void solveJoints(int iterations);        // ball / hinge / fixed / slider
    // Candidate pairs (sorted, i<j) from the broadphase.
    std::vector<std::pair<std::size_t, std::size_t>> broadphasePairs() const;
    void updateSleep();    // put settled bodies to sleep
    void applyWalls();     // reflective boundary
    void wrapPositions();  // periodic boundary
};

}  // namespace ne

#endif  // NATIVEENGINE_WORLD_HPP
