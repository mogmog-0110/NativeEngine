#pragma once
#ifndef NATIVEENGINE_PHYSICS_WORLD_HPP
#define NATIVEENGINE_PHYSICS_WORLD_HPP

#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

#include "world.hpp"

namespace ne {

// Embeddable game-facing facade over World.
//
// Games create and destroy bodies at runtime, so they need STABLE handles, not
// vector indices that shift on removal. This mirrors the shape MitiruEngine's
// backends expose (cf. JoltPhysicsWorld): integer BodyId handles, create/remove,
// pose get/set, impulse, and a contact callback. The heavy solver stays in
// World; this layer only owns the handle bookkeeping.
using BodyId = std::uint32_t;
inline constexpr BodyId kInvalidBody = 0;

struct ContactInfo {
    BodyId a = kInvalidBody, b = kInvalidBody;
    V3 point;
    V3 normal;      // from b toward a
    double depth = 0.0;
};

class PhysicsWorld {
public:
    // --- configuration (forwarded to the solver) ---
    void setBox(double half, bool periodic) { w_.box.half = half; w_.box.periodic = periodic; }
    void setTimestep(double dt) { w_.dt = dt; }
    void setRestitution(double e) { w_.restitution = e; }
    void setFriction(double mu) { w_.friction = mu; }
    void setGravity(const V3& g) { w_.gravity = g; }
    void setSleepEnabled(bool on) { w_.sleepEnabled = on; }
    double timestep() const { return w_.dt; }

    // --- body creation (returns a stable handle; kInvalidBody on failure) ---
    BodyId addSphere(const V3& pos, double radius, double density) {
        return add(makeSphere(0, pos, radius, density));
    }
    BodyId addBox(const V3& pos, const Q& rot, const V3& halfExtents, double density) {
        return add(makeBox(0, pos, rot, halfExtents, density));
    }
    BodyId addCylinder(const V3& pos, const Q& rot, double radius, double height, double density) {
        return add(makeCylinder(0, pos, rot, radius, height, density));
    }
    // A static (infinite-mass) body: never integrates, never sleeps/wakes,
    // collides as immovable.
    BodyId makeStatic(BodyId h) {
        Body* b = body(h);
        if (b) {
            b->invMass = 0.0; b->invInertiaBody = {0, 0, 0};
            b->v = {}; b->w = {};
            b->dynamic = false;
        }
        return h;
    }

    void remove(BodyId h) {
        auto it = slot_.find(h);
        if (it == slot_.end()) return;
        size_t s = it->second, last = w_.bodies.size() - 1;
        if (s != last) {                       // swap-and-pop, keep handles valid
            w_.bodies[s] = w_.bodies[last];
            BodyId moved = handle_[last];
            handle_[s] = moved;
            slot_[moved] = s;
        }
        w_.bodies.pop_back();
        handle_.pop_back();
        slot_.erase(it);
    }

    // --- stepping ---
    void step() {
        w_.step();
        if (onContact_) emitContacts();
    }

    // --- accessors ---
    bool valid(BodyId h) const { return slot_.count(h) != 0; }
    bool isSleeping(BodyId h) const { const Body* b = body(h); return b && b->sleeping; }
    V3 position(BodyId h) const { const Body* b = body(h); return b ? b->x : V3{}; }
    Q  orientation(BodyId h) const { const Body* b = body(h); return b ? b->q : Q{}; }
    V3 velocity(BodyId h) const { const Body* b = body(h); return b ? b->v : V3{}; }
    V3 angularVelocity(BodyId h) const { const Body* b = body(h); return b ? b->w : V3{}; }

    // Setters wake the body: an externally-driven body must not stay frozen.
    void setPosition(BodyId h, const V3& p) { if (Body* b = body(h)) { World::wake(*b); b->x = p; } }
    void setOrientation(BodyId h, const Q& q) { if (Body* b = body(h)) { World::wake(*b); b->q = q; } }
    void setVelocity(BodyId h, const V3& v) { if (Body* b = body(h)) { World::wake(*b); b->v = v; } }
    void setAngularVelocity(BodyId h, const V3& w) { if (Body* b = body(h)) { World::wake(*b); b->w = w; } }

    void applyForce(BodyId h, const V3& f) { if (Body* b = body(h)) { World::wake(*b); b->force += f; } }
    void applyImpulse(BodyId h, const V3& imp) {
        if (Body* b = body(h)) { World::wake(*b); b->v += imp * b->invMass; }
    }
    // Impulse at a world point (produces torque).
    void applyImpulseAt(BodyId h, const V3& imp, const V3& point) {
        if (Body* b = body(h)) {
            World::wake(*b);
            b->v += imp * b->invMass;
            b->w += b->applyInvInertiaWorld((point - b->x).cross(imp));
        }
    }

    void setContactCallback(std::function<void(const ContactInfo&)> cb) {
        onContact_ = std::move(cb);
    }

    std::size_t bodyCount() const { return w_.bodies.size(); }
    World& world() { return w_; }   // escape hatch for advanced use

private:
    World w_;
    std::vector<BodyId> handle_;                    // parallel to w_.bodies
    std::unordered_map<BodyId, std::size_t> slot_;  // handle -> index
    BodyId next_ = 1;                               // 0 reserved for "invalid"
    std::function<void(const ContactInfo&)> onContact_;

    BodyId add(Body b) {
        BodyId h = next_++;
        slot_[h] = w_.bodies.size();
        handle_.push_back(h);
        w_.bodies.push_back(b);
        return h;
    }
    Body* body(BodyId h) {
        auto it = slot_.find(h);
        return (it != slot_.end()) ? &w_.bodies[it->second] : nullptr;
    }
    const Body* body(BodyId h) const {
        auto it = slot_.find(h);
        return (it != slot_.end()) ? &w_.bodies[it->second] : nullptr;
    }

    // Report overlapping pairs this step (positions are post-step). Simple and
    // deterministic; a persistent-manifold version comes with warm starting.
    void emitContacts();
};

}  // namespace ne

#endif  // NATIVEENGINE_PHYSICS_WORLD_HPP
