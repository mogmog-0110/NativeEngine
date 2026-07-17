#pragma once
#ifndef NATIVEENGINE_PHYSICS_WORLD_HPP
#define NATIVEENGINE_PHYSICS_WORLD_HPP

#include <cstdint>
#include <functional>
#include <map>
#include <unordered_map>
#include <utility>
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
    std::uint64_t aUser = 0, bUser = 0;   // the bodies' userData (entity ids)
    V3 point;
    V3 normal;      // from b toward a
    double depth = 0.0;
};

// Result of a scene raycast: the nearest body the ray hits, if any.
struct RayHit {
    bool hit = false;
    BodyId body = kInvalidBody;
    std::uint64_t userData = 0;
    V3 point;
    V3 normal;         // outward surface normal at the hit
    double distance = 0.0;
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
    BodyId addCapsule(const V3& pos, const Q& rot, double radius, double halfHeight, double density) {
        return add(makeCapsule(0, pos, rot, radius, halfHeight, density));
    }
    BodyId addConvex(const V3& pos, const Q& rot, const std::vector<V3>& verts, double density) {
        return add(makeConvex(0, pos, rot, verts, density));
    }
    // A compound of primitive children given relative to the centre of mass.
    BodyId addCompound(const V3& pos, const Q& rot, const std::vector<ChildShape>& children, double density) {
        return add(makeCompound(0, pos, rot, children, density));
    }
    // An infinite static plane (half-space) through `point` with outward `normal`.
    BodyId addPlane(const V3& point, const V3& normal) {
        return add(makePlane(0, point, normal));
    }
    // A static triangle mesh (level geometry): `verts` + flat `indices` (triples).
    // Builds the BVH once here.
    BodyId addMesh(const V3& pos, const Q& rot, const std::vector<V3>& verts,
                   const std::vector<int>& indices) {
        auto md = std::make_shared<MeshData>();
        md->verts = verts;
        for (std::size_t i = 0; i + 2 < indices.size(); i += 3)
            md->tris.push_back({indices[i], indices[i + 1], indices[i + 2]});
        md->build();
        return add(makeMesh(0, pos, rot, md));
    }
    // A heightfield (terrain): a rows x cols grid of heights, spacing `spacing` in
    // x/z, triangulated into a static mesh. heights is row-major [r*cols + c].
    BodyId addHeightfield(const V3& origin, double spacing, int rows, int cols,
                          const std::vector<double>& heights) {
        std::vector<V3> verts; verts.reserve((std::size_t)rows * cols);
        for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c)
                verts.push_back(origin + V3{c * spacing, heights[(std::size_t)r * cols + c], r * spacing});
        std::vector<int> idx;
        for (int r = 0; r + 1 < rows; ++r)
            for (int c = 0; c + 1 < cols; ++c) {
                int i00 = r * cols + c, i01 = i00 + 1, i10 = i00 + cols, i11 = i10 + 1;
                idx.insert(idx.end(), {i00, i10, i11, i00, i11, i01});
            }
        return addMesh(V3{0, 0, 0}, Q{1, 0, 0, 0}, verts, idx);
    }
    // A static (infinite-mass) body: never integrates, never sleeps/wakes,
    // collides as immovable.
    BodyId makeStatic(BodyId h) {
        Body* b = body(h);
        if (b) {
            b->invMass = 0.0; b->invInertiaBody = {0, 0, 0};
            b->v = {}; b->w = {};
            b->dynamic = false; b->kinematic = false;
        }
        return h;
    }
    // Kinematic: infinite mass, but a real velocity that pushes dynamic bodies and
    // advances its own pose. Drive it with setVelocity/setAngularVelocity, or with
    // setKinematicTarget to move it to a pose over the next step (the implied
    // velocity is what actually pushes contacts, so use this, not setPosition).
    BodyId makeKinematic(BodyId h) {
        Body* b = body(h);
        if (b) {
            b->invMass = 0.0; b->invInertiaBody = {0, 0, 0};
            b->dynamic = false; b->kinematic = true; b->sleeping = false;
        }
        return h;
    }
    void setKinematicTarget(BodyId h, const V3& pos) {
        if (Body* b = body(h)) b->v = (pos - b->x) * (1.0 / w_.dt);
    }
    void setKinematicTarget(BodyId h, const V3& pos, const Q& rot) {
        if (Body* b = body(h)) {
            b->v = (pos - b->x) * (1.0 / w_.dt);
            b->w = angularVelocityBetween(b->q, rot, w_.dt);
        }
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

    // --- joints ---
    // A distance joint (rigid rod / rope when stiffness==0, else a spring-damper)
    // between two bodies at local anchors. Stored by handle so it stays valid
    // across body removal; anchor to a fixed point by joining to a static body.
    void addDistanceJoint(BodyId a, BodyId b, const V3& localA, const V3& localB,
                          double rest, double stiffness = 0.0, double damping = 0.0) {
        joints_.push_back({a, b, localA, localB, rest, stiffness, damping});
    }

    // General joints: anchors/axes given in WORLD space at creation (converted to
    // each body's local frame). Returns a JointId for later configuration/queries.
    using JointId = std::uint32_t;
    JointId addBallJoint(BodyId a, BodyId b, const V3& worldAnchor) {
        return addGeneralJoint(JointType::Ball, a, b, worldAnchor, V3{0, 1, 0});
    }
    JointId addHingeJoint(BodyId a, BodyId b, const V3& worldAnchor, const V3& worldAxis) {
        return addGeneralJoint(JointType::Hinge, a, b, worldAnchor, worldAxis);
    }
    JointId addSliderJoint(BodyId a, BodyId b, const V3& worldAnchor, const V3& worldAxis) {
        return addGeneralJoint(JointType::Slider, a, b, worldAnchor, worldAxis);
    }
    JointId addFixedJoint(BodyId a, BodyId b) {
        const Body* ba = body(a); const Body* bb = body(b);
        V3 anchor = (ba && bb) ? (ba->x + bb->x) * 0.5 : V3{};
        return addGeneralJoint(JointType::Fixed, a, b, anchor, V3{0, 1, 0});
    }
    void setJointMotor(JointId id, double speed, double maxImpulse) {
        if (id < gjoints_.size()) { auto& j = gjoints_[id]; j.useMotor = true; j.motorSpeed = speed; j.maxMotorImpulse = maxImpulse; }
    }
    void setJointLimit(JointId id, double lower, double upper) {
        if (id < gjoints_.size()) { auto& j = gjoints_[id]; j.useLimit = true; j.lower = lower; j.upper = upper; }
    }
    void setJointBreakable(JointId id, double breakImpulse) {
        if (id < gjoints_.size()) { auto& j = gjoints_[id]; j.breakable = true; j.breakImpulse = breakImpulse; }
    }
    bool jointBroken(JointId id) const { return id < gjoints_.size() && gjoints_[id].broken; }

    // --- stepping ---
    void step() {
        // Resolve handle-based joints to current slot indices (indices move on
        // removal); drop joints whose bodies were removed.
        w_.distanceJoints.clear();
        for (const FacadeJoint& fj : joints_) {
            auto ia = slot_.find(fj.a), ib = slot_.find(fj.b);
            if (ia == slot_.end() || ib == slot_.end()) continue;
            DistanceJoint dj;
            dj.a = ia->second; dj.b = ib->second;
            dj.localA = fj.localA; dj.localB = fj.localB;
            dj.rest = fj.rest; dj.stiffness = fj.stiffness; dj.damping = fj.damping;
            w_.distanceJoints.push_back(dj);
        }
        // Rebuild general joints with current slot indices; drop joints whose
        // bodies were removed. Preserve the persistent runtime state (broken).
        w_.joints.clear();
        jointMap_.clear();
        for (std::size_t gi = 0; gi < gjoints_.size(); ++gi) {
            Joint jj = gjoints_[gi];   // a/b currently hold HANDLES
            auto ia = slot_.find((BodyId)jj.a), ib = slot_.find((BodyId)jj.b);
            if (ia == slot_.end() || ib == slot_.end()) continue;
            jj.a = ia->second; jj.b = ib->second;
            w_.joints.push_back(jj);
            jointMap_.push_back(gi);
        }
        w_.step();
        for (std::size_t k = 0; k < w_.joints.size(); ++k)
            gjoints_[jointMap_[k]].broken = w_.joints[k].broken;   // persist break state
        if (wantsEvents()) processContacts();
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
    void applyTorque(BodyId h, const V3& t) { if (Body* b = body(h)) { World::wake(*b); b->torque += t; } }
    // Force applied at a world point over the step (produces torque too).
    void applyForceAtPoint(BodyId h, const V3& f, const V3& point) {
        if (Body* b = body(h)) { World::wake(*b); b->force += f; b->torque += (point - b->x).cross(f); }
    }
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
    void applyAngularImpulse(BodyId h, const V3& angImp) {
        if (Body* b = body(h)) { World::wake(*b); b->w += b->applyInvInertiaWorld(angImp); }
    }
    void setGravityScale(BodyId h, double s) { if (Body* b = body(h)) b->gravityScale = s; }

    // --- render interpolation ---
    // Pose blended between the pre-step and post-step states; call with the render
    // accumulator ratio alpha in [0,1] for judder-free rendering at a fixed step.
    V3 interpolatedPosition(BodyId h, double alpha) const {
        const Body* b = body(h); return b ? lerp(b->prevX, b->x, alpha) : V3{};
    }
    Q interpolatedOrientation(BodyId h, double alpha) const {
        const Body* b = body(h); return b ? nlerp(b->prevQ, b->q, alpha) : Q{};
    }

    // --- debug draw hook (feeds the embedder's own renderer) ---
    enum DebugFlags {
        DbgColliders  = 1 << 0,
        DbgAABB       = 1 << 1,
        DbgContacts   = 1 << 2,   // requires contacts capture (enabled here)
        DbgVelocities = 1 << 3,
        DbgJoints     = 1 << 4,
        DbgCOM        = 1 << 5,
        DbgAll        = 0x3f,
    };
    using LineFn = std::function<void(const V3& a, const V3& b, const V3& color)>;
    void debugDraw(const LineFn& line, std::uint32_t flags = DbgColliders) const;

    // --- per-body properties (game backend surface) ---
    // Material: friction / restitution for this body; combined with the partner's
    // at each contact (World combine modes). A negative value inherits the world
    // default.
    void setMaterial(BodyId h, double friction, double restitution) {
        if (Body* b = body(h)) { b->friction = friction; b->restitution = restitution; }
    }
    void setDamping(BodyId h, double linear, double angular) {
        if (Body* b = body(h)) { b->linearDamping = linear; b->angularDamping = angular; }
    }
    // Collision filter: `layer` is an index 0..31; `mask` is the set of layers
    // this body collides with. Two bodies collide iff each mask holds the other's
    // layer.
    void setLayerMask(BodyId h, std::uint32_t layer, std::uint32_t mask) {
        if (Body* b = body(h)) { b->layer = layer; b->mask = mask; }
    }
    void setUserData(BodyId h, std::uint64_t data) { if (Body* b = body(h)) b->userData = data; }
    std::uint64_t userData(BodyId h) const { const Body* b = body(h); return b ? b->userData : 0; }
    // A sensor detects overlap but has no physical response, and reports through
    // the trigger enter/stay/exit callbacks instead of the contact ones.
    void setSensor(BodyId h, bool on) { if (Body* b = body(h)) b->sensor = on; }
    // Continuous collision for a fast body (won't tunnel through thin geometry).
    void setCcd(BodyId h, bool on) { if (Body* b = body(h)) b->ccd = on; }
    // Enable capture of solver contact points (needed for DbgContacts debug draw).
    void setDebugContactCapture(bool on) { w_.captureContacts = on; }

    // --- events ---
    // Legacy: fires once per solid overlapping pair every step.
    void setContactCallback(std::function<void(const ContactInfo&)> cb) { onContact_ = std::move(cb); }
    // Transitions (keyed by body pair across frames): the events gameplay wants.
    using ContactCb = std::function<void(const ContactInfo&)>;
    void setContactBeginCallback(ContactCb cb) { onBegin_ = std::move(cb); }
    void setContactStayCallback(ContactCb cb) { onStay_ = std::move(cb); }
    void setContactEndCallback(ContactCb cb) { onEnd_ = std::move(cb); }
    void setTriggerEnterCallback(ContactCb cb) { onTrigEnter_ = std::move(cb); }
    void setTriggerStayCallback(ContactCb cb) { onTrigStay_ = std::move(cb); }
    void setTriggerExitCallback(ContactCb cb) { onTrigExit_ = std::move(cb); }

    // --- scene queries (layer-filtered, PBC-aware) ---
    // Nearest body hit by the ray from origin along dir (need not be unit) up to
    // maxDist. `mask` is the set of layers to test (bit set => layer included).
    RayHit raycast(const V3& origin, const V3& dir, double maxDist,
                   std::uint32_t mask = 0xFFFFFFFFu) const;
    // All bodies overlapping a query sphere / box (e.g. area-of-effect, placement).
    std::vector<BodyId> overlapSphere(const V3& center, double radius,
                                      std::uint32_t mask = 0xFFFFFFFFu) const;
    std::vector<BodyId> overlapBox(const V3& center, const Q& rot, const V3& half,
                                   std::uint32_t mask = 0xFFFFFFFFu) const;
    // Sweep a sphere of `radius` from `origin` along `dir` for maxDist; returns the
    // first body hit (time-of-impact via conservative advancement). The building
    // block for the character controller and CCD.
    RayHit sphereCast(const V3& origin, double radius, const V3& dir, double maxDist,
                      std::uint32_t mask = 0xFFFFFFFFu) const;

    // Penetration of a query shape (not in the world) against every solid body it
    // overlaps: normal points from the world body toward the query (push-out
    // direction), depth is the overlap. Used by the character controller.
    struct OverlapHit { BodyId body = kInvalidBody; std::uint64_t userData = 0; V3 normal; double depth = 0.0; };
    std::vector<OverlapHit> overlapContacts(const Body& query, std::uint32_t mask = 0xFFFFFFFFu) const;

    std::size_t bodyCount() const { return w_.bodies.size(); }
    World& world() { return w_; }   // escape hatch for advanced use

private:
    World w_;
    std::vector<BodyId> handle_;                    // parallel to w_.bodies
    std::unordered_map<BodyId, std::size_t> slot_;  // handle -> index
    BodyId next_ = 1;                               // 0 reserved for "invalid"
    std::function<void(const ContactInfo&)> onContact_;
    ContactCb onBegin_, onStay_, onEnd_, onTrigEnter_, onTrigStay_, onTrigExit_;
    // Overlapping pairs from the previous step (normalized handle key -> info),
    // split solid vs trigger, so we can emit begin/stay/end and enter/stay/exit.
    using PairKey = std::pair<BodyId, BodyId>;
    std::map<PairKey, ContactInfo> prevContacts_, prevTriggers_;
    bool wantsEvents() const {
        return onContact_ || onBegin_ || onStay_ || onEnd_ ||
               onTrigEnter_ || onTrigStay_ || onTrigExit_;
    }

    struct FacadeJoint {
        BodyId a, b;
        V3 localA, localB;
        double rest, stiffness, damping;
    };
    std::vector<FacadeJoint> joints_;

    // General joints stored with body HANDLES in a/b; rebuilt to indices per step.
    std::vector<Joint> gjoints_;
    std::vector<std::size_t> jointMap_;   // w_.joints index -> gjoints_ index

    JointId addGeneralJoint(JointType type, BodyId a, BodyId b,
                            const V3& worldAnchor, const V3& worldAxis) {
        const Body* ba = body(a); const Body* bb = body(b);
        Joint j;
        j.type = type;
        j.a = a; j.b = b;   // handles (resolved to indices each step)
        if (ba && bb) {
            j.localA = ba->q.inverseRotate(worldAnchor - ba->x);
            j.localB = bb->q.inverseRotate(worldAnchor - bb->x);
            V3 ax = worldAxis.normalized();
            j.axisA = ba->q.inverseRotate(ax);
            j.axisB = bb->q.inverseRotate(ax);
            j.refRel = ba->q * bb->q.conjugate();   // lock current relative orientation
        }
        gjoints_.push_back(j);
        return (JointId)(gjoints_.size() - 1);
    }

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

    // Diff this step's overlapping pairs against the previous step to emit the
    // legacy per-step callback plus contact begin/stay/end and trigger
    // enter/stay/exit. Positions are post-step.
    void processContacts();
};

}  // namespace ne

#endif  // NATIVEENGINE_PHYSICS_WORLD_HPP
