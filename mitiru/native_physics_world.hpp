#pragma once
#ifndef NATIVEENGINE_MITIRU_NATIVE_PHYSICS_WORLD_HPP
#define NATIVEENGINE_MITIRU_NATIVE_PHYSICS_WORLD_HPP

// MitiruEngine-facing physics backend built on NativeEngine.
//
// Mirrors the concrete-World shape MitiruEngine's backends expose (cf.
// JoltPhysicsWorld), speaking sgc types at the boundary (right-handed, Y-up,
// quaternion (x,y,z,w) w-last, SI units). It surfaces NativeEngine's full feature
// set -- shapes (incl. capsule / static mesh / heightfield), dynamic/static/
// kinematic bodies, per-body materials/damping/layers/userData, triggers + contact
// events, scene queries, joints, a capsule character controller, render
// interpolation, debug draw, and snapshotting -- plus the one thing no other
// MitiruEngine backend has: a NATIVE PERIODIC BOUNDARY.
//
// Build with MITIRU_HAS_NATIVEPHYS defined for the real backend; without it, the
// Null stub keeps header-only / no-submodule builds compiling (physics no-ops).

#include <cstdint>
#include <functional>
#include <vector>

#include <sgc/math/Vec3.hpp>
#include <sgc/math/Quaternion.hpp>

namespace mitiru {
namespace nativephys {

using BodyId = unsigned int;
using JointId = unsigned int;
inline constexpr BodyId kInvalidBody = 0;

// --- sgc-typed results / descriptors (available with or without the backend) ---
struct RayHitS {
    bool hit = false;
    BodyId body = kInvalidBody;
    std::uint64_t userData = 0;
    sgc::Vec3f point, normal;
    float distance = 0.0f;
};
struct ContactEventS {
    BodyId a = kInvalidBody, b = kInvalidBody;
    std::uint64_t aUser = 0, bUser = 0;
    sgc::Vec3f point, normal;   // normal points from b toward a
    float depth = 0.0f;
};

// Mirrors MitiruEngine's PhysicsTrait / RigidBodyComponent so a component can be
// turned into a body in one call (see createBody). mass > 0 overrides density to
// hit that mass; friction/restitution < 0 inherit the world defaults.
struct BodyDesc {
    enum class Type { Dynamic, Static, Kinematic };
    enum class Shape { Sphere, Box, Capsule, Cylinder };
    Type type = Type::Dynamic;
    Shape shape = Shape::Box;
    sgc::Vec3f position;
    sgc::Quaternionf rotation{0, 0, 0, 1};
    sgc::Vec3f halfExtents{0.5f, 0.5f, 0.5f};   // Box
    float radius = 0.5f;                        // Sphere / Capsule / Cylinder
    float halfHeight = 0.5f;                    // Capsule / Cylinder
    float density = 1.0f;
    float mass = -1.0f;                         // > 0 overrides density
    float friction = -1.0f, restitution = -1.0f;
    float linearDamping = 0.0f, angularDamping = 0.0f;
    bool isTrigger = false, ccd = false;
    std::uint32_t layer = 0, mask = 0xFFFFFFFFu;
    std::uint64_t userData = 0;
};

using ContactFn = std::function<void(const ContactEventS&)>;
using DebugLineFn = std::function<void(const sgc::Vec3f& a, const sgc::Vec3f& b, const sgc::Vec3f& color)>;

// Debug-draw flags (mirror ne::PhysicsWorld::DebugFlags).
enum DebugFlags {
    DbgColliders = 1 << 0, DbgAABB = 1 << 1, DbgContacts = 1 << 2,
    DbgVelocities = 1 << 3, DbgJoints = 1 << 4, DbgCOM = 1 << 5, DbgAll = 0x3f
};

}  // namespace nativephys
}  // namespace mitiru

#ifdef MITIRU_HAS_NATIVEPHYS

#include "physics_world.hpp"   // ne::PhysicsWorld
#include "character.hpp"       // ne::CharacterController

namespace mitiru {
namespace nativephys {

// -- sgc <-> NativeEngine conversions (float <-> double; quat order differs) --
inline ne::V3 toNe(const sgc::Vec3f& v) { return {double(v.x), double(v.y), double(v.z)}; }
inline sgc::Vec3f toSgc(const ne::V3& v) { return sgc::Vec3f{float(v.x), float(v.y), float(v.z)}; }
inline ne::Q toNe(const sgc::Quaternionf& q) { return ne::Q{double(q.w), double(q.x), double(q.y), double(q.z)}; }
inline sgc::Quaternionf toSgc(const ne::Q& q) { return sgc::Quaternionf{float(q.x), float(q.y), float(q.z), float(q.w)}; }

class NativePhysicsWorld {
public:
    void init() {}
    void shutdown() {}

    // ---- configuration ----
    void setGravity(const sgc::Vec3f& g) { pw_.setGravity(toNe(g)); }
    void setTimestep(float dt) { pw_.setTimestep(dt); }
    void setRestitution(float e) { pw_.setRestitution(e); }
    void setFriction(float mu) { pw_.setFriction(mu); }
    void setSleepEnabled(bool on) { pw_.setSleepEnabled(on); }
    void setPeriodicBox(float half, bool periodic) { pw_.setBox(half, periodic); }   // NATIVE PBC
    void setOpenBoundary(bool on) { pw_.setOpenBoundary(on); }                       // no reflective walls

    // ---- body creation ----
    BodyId createSphere(const sgc::Vec3f& p, float r, float density) { return pw_.addSphere(toNe(p), r, density); }
    BodyId createBox(const sgc::Vec3f& p, const sgc::Quaternionf& q, const sgc::Vec3f& he, float density) { return pw_.addBox(toNe(p), toNe(q), toNe(he), density); }
    BodyId createCylinder(const sgc::Vec3f& p, const sgc::Quaternionf& q, float r, float h, float density) { return pw_.addCylinder(toNe(p), toNe(q), r, h, density); }
    BodyId createCapsule(const sgc::Vec3f& p, const sgc::Quaternionf& q, float r, float halfHeight, float density) { return pw_.addCapsule(toNe(p), toNe(q), r, halfHeight, density); }
    BodyId createPlane(const sgc::Vec3f& point, const sgc::Vec3f& normal) { return pw_.addPlane(toNe(point), toNe(normal)); }
    BodyId createStaticMesh(const std::vector<sgc::Vec3f>& verts, const std::vector<int>& indices) {
        std::vector<ne::V3> v; v.reserve(verts.size());
        for (const auto& p : verts) v.push_back(toNe(p));
        return pw_.addMesh({0, 0, 0}, ne::Q{1, 0, 0, 0}, v, indices);
    }
    BodyId createHeightfield(const sgc::Vec3f& origin, float spacing, int rows, int cols, const std::vector<double>& heights) {
        return pw_.addHeightfield(toNe(origin), spacing, rows, cols, heights);
    }
    // Component-driven creation (the PhysicsTrait / RigidBodyComponent adapter).
    BodyId createBody(const BodyDesc& d) {
        float dens = d.density;
        if (d.mass > 0.0f) { float vol = volumeOf(d); if (vol > 1e-9f) dens = d.mass / vol; }
        BodyId h = kInvalidBody;
        switch (d.shape) {
            case BodyDesc::Shape::Sphere:   h = createSphere(d.position, d.radius, dens); break;
            case BodyDesc::Shape::Box:      h = createBox(d.position, d.rotation, d.halfExtents, dens); break;
            case BodyDesc::Shape::Capsule:  h = createCapsule(d.position, d.rotation, d.radius, d.halfHeight, dens); break;
            case BodyDesc::Shape::Cylinder: h = createCylinder(d.position, d.rotation, d.radius, 2.0f * d.halfHeight, dens); break;
        }
        if (d.type == BodyDesc::Type::Static) pw_.makeStatic(h);
        else if (d.type == BodyDesc::Type::Kinematic) pw_.makeKinematic(h);
        pw_.setMaterial(h, d.friction, d.restitution);
        pw_.setDamping(h, d.linearDamping, d.angularDamping);
        pw_.setLayerMask(h, d.layer, d.mask);
        pw_.setUserData(h, d.userData);
        pw_.setSensor(h, d.isTrigger);
        pw_.setCcd(h, d.ccd);
        return h;
    }
    BodyId makeStatic(BodyId h) { return pw_.makeStatic(h); }
    BodyId makeKinematic(BodyId h) { return pw_.makeKinematic(h); }
    void setKinematicTarget(BodyId h, const sgc::Vec3f& pos) { pw_.setKinematicTarget(h, toNe(pos)); }
    void setKinematicTarget(BodyId h, const sgc::Vec3f& pos, const sgc::Quaternionf& rot) { pw_.setKinematicTarget(h, toNe(pos), toNe(rot)); }
    void removeBody(BodyId h) { pw_.remove(h); }

    // ---- per-body properties ----
    void setMaterial(BodyId h, float friction, float restitution) { pw_.setMaterial(h, friction, restitution); }
    void setDamping(BodyId h, float lin, float ang) { pw_.setDamping(h, lin, ang); }
    void setGravityScale(BodyId h, float s) { pw_.setGravityScale(h, s); }
    void setLayerMask(BodyId h, std::uint32_t layer, std::uint32_t mask) { pw_.setLayerMask(h, layer, mask); }
    void setUserData(BodyId h, std::uint64_t d) { pw_.setUserData(h, d); }
    std::uint64_t getUserData(BodyId h) const { return pw_.userData(h); }
    void setSensor(BodyId h, bool on) { pw_.setSensor(h, on); }
    void setCcd(BodyId h, bool on) { pw_.setCcd(h, on); }

    // ---- stepping (mirrors JoltPhysicsWorld::update(dt)) ----
    void update(float dt) { pw_.setTimestep(dt); pw_.step(); }

    // ---- pose / velocity ----
    sgc::Vec3f getPosition(BodyId h) const { return toSgc(pw_.position(h)); }
    sgc::Quaternionf getRotation(BodyId h) const { return toSgc(pw_.orientation(h)); }
    sgc::Vec3f getVelocity(BodyId h) const { return toSgc(pw_.velocity(h)); }
    sgc::Vec3f getAngularVelocity(BodyId h) const { return toSgc(pw_.angularVelocity(h)); }
    void setPosition(BodyId h, const sgc::Vec3f& p) { pw_.setPosition(h, toNe(p)); }
    void setRotation(BodyId h, const sgc::Quaternionf& q) { pw_.setOrientation(h, toNe(q)); }
    void setVelocity(BodyId h, const sgc::Vec3f& v) { pw_.setVelocity(h, toNe(v)); }
    void setAngularVelocity(BodyId h, const sgc::Vec3f& w) { pw_.setAngularVelocity(h, toNe(w)); }
    // Interpolated pose for smooth rendering at a fixed step (alpha in [0,1]).
    sgc::Vec3f interpolatedPosition(BodyId h, float alpha) const { return toSgc(pw_.interpolatedPosition(h, alpha)); }
    sgc::Quaternionf interpolatedRotation(BodyId h, float alpha) const { return toSgc(pw_.interpolatedOrientation(h, alpha)); }

    // ---- forces ----
    void applyForce(BodyId h, const sgc::Vec3f& f) { pw_.applyForce(h, toNe(f)); }
    void applyTorque(BodyId h, const sgc::Vec3f& t) { pw_.applyTorque(h, toNe(t)); }
    void applyForceAtPoint(BodyId h, const sgc::Vec3f& f, const sgc::Vec3f& p) { pw_.applyForceAtPoint(h, toNe(f), toNe(p)); }
    void applyImpulse(BodyId h, const sgc::Vec3f& imp) { pw_.applyImpulse(h, toNe(imp)); }
    void applyImpulseAt(BodyId h, const sgc::Vec3f& imp, const sgc::Vec3f& p) { pw_.applyImpulseAt(h, toNe(imp), toNe(p)); }
    void applyAngularImpulse(BodyId h, const sgc::Vec3f& j) { pw_.applyAngularImpulse(h, toNe(j)); }

    // ---- joints ----
    void addDistanceJoint(BodyId a, BodyId b, const sgc::Vec3f& la, const sgc::Vec3f& lb, float rest, float stiffness = 0.f, float damping = 0.f) { pw_.addDistanceJoint(a, b, toNe(la), toNe(lb), rest, stiffness, damping); }
    JointId addBallJoint(BodyId a, BodyId b, const sgc::Vec3f& anchor) { return pw_.addBallJoint(a, b, toNe(anchor)); }
    JointId addHingeJoint(BodyId a, BodyId b, const sgc::Vec3f& anchor, const sgc::Vec3f& axis) { return pw_.addHingeJoint(a, b, toNe(anchor), toNe(axis)); }
    JointId addFixedJoint(BodyId a, BodyId b) { return pw_.addFixedJoint(a, b); }
    JointId addSliderJoint(BodyId a, BodyId b, const sgc::Vec3f& anchor, const sgc::Vec3f& axis) { return pw_.addSliderJoint(a, b, toNe(anchor), toNe(axis)); }
    void setJointMotor(JointId j, float speed, float maxImpulse) { pw_.setJointMotor(j, speed, maxImpulse); }
    void setJointLimit(JointId j, float lower, float upper) { pw_.setJointLimit(j, lower, upper); }
    void setJointBreakable(JointId j, float breakImpulse) { pw_.setJointBreakable(j, breakImpulse); }
    bool jointBroken(JointId j) const { return pw_.jointBroken(j); }

    // ---- events ----
    void onContact(ContactFn cb) { onContact_ = std::move(cb); pw_.setContactCallback(wrap(&onContact_)); }
    void onContactBegin(ContactFn cb) { onBegin_ = std::move(cb); pw_.setContactBeginCallback(wrap(&onBegin_)); }
    void onContactStay(ContactFn cb) { onStay_ = std::move(cb); pw_.setContactStayCallback(wrap(&onStay_)); }
    void onContactEnd(ContactFn cb) { onEnd_ = std::move(cb); pw_.setContactEndCallback(wrap(&onEnd_)); }
    void onTriggerEnter(ContactFn cb) { onTrigEnter_ = std::move(cb); pw_.setTriggerEnterCallback(wrap(&onTrigEnter_)); }
    void onTriggerStay(ContactFn cb) { onTrigStay_ = std::move(cb); pw_.setTriggerStayCallback(wrap(&onTrigStay_)); }
    void onTriggerExit(ContactFn cb) { onTrigExit_ = std::move(cb); pw_.setTriggerExitCallback(wrap(&onTrigExit_)); }

    // ---- scene queries (PBC-aware, layer-filtered) ----
    RayHitS raycast(const sgc::Vec3f& o, const sgc::Vec3f& dir, float maxDist, std::uint32_t mask = 0xFFFFFFFFu) const { return toHit(pw_.raycast(toNe(o), toNe(dir), maxDist, mask)); }
    RayHitS sphereCast(const sgc::Vec3f& o, float radius, const sgc::Vec3f& dir, float maxDist, std::uint32_t mask = 0xFFFFFFFFu) const { return toHit(pw_.sphereCast(toNe(o), radius, toNe(dir), maxDist, mask)); }
    std::vector<BodyId> overlapSphere(const sgc::Vec3f& c, float r, std::uint32_t mask = 0xFFFFFFFFu) const { return pw_.overlapSphere(toNe(c), r, mask); }

    // ---- debug draw (feeds MitiruEngine's own renderer) ----
    void setDebugContactCapture(bool on) { pw_.setDebugContactCapture(on); }
    void debugDraw(const DebugLineFn& line, std::uint32_t flags = DbgColliders) const {
        pw_.debugDraw([&](const ne::V3& a, const ne::V3& b, const ne::V3& col) { line(toSgc(a), toSgc(b), toSgc(col)); }, flags);
    }

    // ---- snapshot (rollback / checkpoint) ----
    std::vector<char> saveState() const { return pw_.saveState(); }
    bool loadState(const std::vector<char>& data) { return pw_.loadState(data); }

    bool isSleeping(BodyId h) const { return pw_.isSleeping(h); }
    std::size_t bodyCount() const { return pw_.bodyCount(); }
    ne::PhysicsWorld& native() { return pw_; }   // escape hatch (needed by NativeCharacter)

private:
    static float volumeOf(const BodyDesc& d) {
        const float PI = 3.14159265f;
        switch (d.shape) {
            case BodyDesc::Shape::Sphere:   return (4.0f / 3.0f) * PI * d.radius * d.radius * d.radius;
            case BodyDesc::Shape::Box:      return 8.0f * d.halfExtents.x * d.halfExtents.y * d.halfExtents.z;
            case BodyDesc::Shape::Capsule:  return PI * d.radius * d.radius * (2.0f * d.halfHeight) + (4.0f / 3.0f) * PI * d.radius * d.radius * d.radius;
            case BodyDesc::Shape::Cylinder: return PI * d.radius * d.radius * (2.0f * d.halfHeight);
        }
        return 1.0f;
    }
    ContactEventS toEvent(const ne::ContactInfo& ci) const {
        return {ci.a, ci.b, ci.aUser, ci.bUser, toSgc(ci.point), toSgc(ci.normal), float(ci.depth)};
    }
    RayHitS toHit(const ne::RayHit& r) const {
        RayHitS h; h.hit = r.hit; h.body = r.body; h.userData = r.userData;
        h.point = toSgc(r.point); h.normal = toSgc(r.normal); h.distance = float(r.distance);
        return h;
    }
    // Bind a stable pointer to one of the member callbacks (they outlive the
    // solver's lambda; capturing the parameter by reference would dangle).
    std::function<void(const ne::ContactInfo&)> wrap(const ContactFn* fn) {
        return [this, fn](const ne::ContactInfo& ci) { if (*fn) (*fn)(toEvent(ci)); };
    }

    ne::PhysicsWorld pw_;
    ContactFn onContact_, onBegin_, onStay_, onEnd_, onTrigEnter_, onTrigStay_, onTrigExit_;
};

// A kinematic capsule character bound to a world (collide-and-slide).
class NativeCharacter {
public:
    NativeCharacter(NativePhysicsWorld& w, float radius, float halfHeight, std::uint32_t mask = 0xFFFFFFFFu)
        : cc_(&w.native(), radius, halfHeight, mask) {}
    void setPosition(const sgc::Vec3f& p) { cc_.setPosition(toNe(p)); }
    sgc::Vec3f position() const { return toSgc(cc_.position()); }
    void move(const sgc::Vec3f& disp) { cc_.move(toNe(disp)); }
    bool grounded() const { return cc_.grounded(); }
    void setSlopeLimit(float cosAngle) { cc_.setSlopeLimit(cosAngle); }
private:
    ne::CharacterController cc_;
};

}  // namespace nativephys
}  // namespace mitiru

#else  // !MITIRU_HAS_NATIVEPHYS -- Null stub with identical signatures

namespace mitiru {
namespace nativephys {

class NativePhysicsWorld {
public:
    void init() {}
    void shutdown() {}
    void setGravity(const sgc::Vec3f&) {}
    void setTimestep(float) {}
    void setRestitution(float) {}
    void setFriction(float) {}
    void setSleepEnabled(bool) {}
    void setPeriodicBox(float, bool) {}
    void setOpenBoundary(bool) {}
    BodyId createSphere(const sgc::Vec3f&, float, float) { return kInvalidBody; }
    BodyId createBox(const sgc::Vec3f&, const sgc::Quaternionf&, const sgc::Vec3f&, float) { return kInvalidBody; }
    BodyId createCylinder(const sgc::Vec3f&, const sgc::Quaternionf&, float, float, float) { return kInvalidBody; }
    BodyId createCapsule(const sgc::Vec3f&, const sgc::Quaternionf&, float, float, float) { return kInvalidBody; }
    BodyId createPlane(const sgc::Vec3f&, const sgc::Vec3f&) { return kInvalidBody; }
    BodyId createStaticMesh(const std::vector<sgc::Vec3f>&, const std::vector<int>&) { return kInvalidBody; }
    BodyId createHeightfield(const sgc::Vec3f&, float, int, int, const std::vector<double>&) { return kInvalidBody; }
    BodyId createBody(const BodyDesc&) { return kInvalidBody; }
    BodyId makeStatic(BodyId h) { return h; }
    BodyId makeKinematic(BodyId h) { return h; }
    void setKinematicTarget(BodyId, const sgc::Vec3f&) {}
    void setKinematicTarget(BodyId, const sgc::Vec3f&, const sgc::Quaternionf&) {}
    void removeBody(BodyId) {}
    void setMaterial(BodyId, float, float) {}
    void setDamping(BodyId, float, float) {}
    void setGravityScale(BodyId, float) {}
    void setLayerMask(BodyId, std::uint32_t, std::uint32_t) {}
    void setUserData(BodyId, std::uint64_t) {}
    std::uint64_t getUserData(BodyId) const { return 0; }
    void setSensor(BodyId, bool) {}
    void setCcd(BodyId, bool) {}
    void update(float) {}
    sgc::Vec3f getPosition(BodyId) const { return {}; }
    sgc::Quaternionf getRotation(BodyId) const { return {}; }
    sgc::Vec3f getVelocity(BodyId) const { return {}; }
    sgc::Vec3f getAngularVelocity(BodyId) const { return {}; }
    void setPosition(BodyId, const sgc::Vec3f&) {}
    void setRotation(BodyId, const sgc::Quaternionf&) {}
    void setVelocity(BodyId, const sgc::Vec3f&) {}
    void setAngularVelocity(BodyId, const sgc::Vec3f&) {}
    sgc::Vec3f interpolatedPosition(BodyId, float) const { return {}; }
    sgc::Quaternionf interpolatedRotation(BodyId, float) const { return {}; }
    void applyForce(BodyId, const sgc::Vec3f&) {}
    void applyTorque(BodyId, const sgc::Vec3f&) {}
    void applyForceAtPoint(BodyId, const sgc::Vec3f&, const sgc::Vec3f&) {}
    void applyImpulse(BodyId, const sgc::Vec3f&) {}
    void applyImpulseAt(BodyId, const sgc::Vec3f&, const sgc::Vec3f&) {}
    void applyAngularImpulse(BodyId, const sgc::Vec3f&) {}
    void addDistanceJoint(BodyId, BodyId, const sgc::Vec3f&, const sgc::Vec3f&, float, float = 0.f, float = 0.f) {}
    JointId addBallJoint(BodyId, BodyId, const sgc::Vec3f&) { return 0; }
    JointId addHingeJoint(BodyId, BodyId, const sgc::Vec3f&, const sgc::Vec3f&) { return 0; }
    JointId addFixedJoint(BodyId, BodyId) { return 0; }
    JointId addSliderJoint(BodyId, BodyId, const sgc::Vec3f&, const sgc::Vec3f&) { return 0; }
    void setJointMotor(JointId, float, float) {}
    void setJointLimit(JointId, float, float) {}
    void setJointBreakable(JointId, float) {}
    bool jointBroken(JointId) const { return false; }
    void onContact(ContactFn) {}
    void onContactBegin(ContactFn) {}
    void onContactStay(ContactFn) {}
    void onContactEnd(ContactFn) {}
    void onTriggerEnter(ContactFn) {}
    void onTriggerStay(ContactFn) {}
    void onTriggerExit(ContactFn) {}
    RayHitS raycast(const sgc::Vec3f&, const sgc::Vec3f&, float, std::uint32_t = 0xFFFFFFFFu) const { return {}; }
    RayHitS sphereCast(const sgc::Vec3f&, float, const sgc::Vec3f&, float, std::uint32_t = 0xFFFFFFFFu) const { return {}; }
    std::vector<BodyId> overlapSphere(const sgc::Vec3f&, float, std::uint32_t = 0xFFFFFFFFu) const { return {}; }
    void setDebugContactCapture(bool) {}
    void debugDraw(const DebugLineFn&, std::uint32_t = DbgColliders) const {}
    std::vector<char> saveState() const { return {}; }
    bool loadState(const std::vector<char>&) { return false; }
    bool isSleeping(BodyId) const { return false; }
    std::size_t bodyCount() const { return 0; }
};

class NativeCharacter {
public:
    NativeCharacter(NativePhysicsWorld&, float, float, std::uint32_t = 0xFFFFFFFFu) {}
    void setPosition(const sgc::Vec3f&) {}
    sgc::Vec3f position() const { return {}; }
    void move(const sgc::Vec3f&) {}
    bool grounded() const { return false; }
    void setSlopeLimit(float) {}
};

}  // namespace nativephys
}  // namespace mitiru

#endif  // MITIRU_HAS_NATIVEPHYS

#endif  // NATIVEENGINE_MITIRU_NATIVE_PHYSICS_WORLD_HPP
