#pragma once
#ifndef NATIVEENGINE_MITIRU_NATIVE_PHYSICS_WORLD_HPP
#define NATIVEENGINE_MITIRU_NATIVE_PHYSICS_WORLD_HPP

// MitiruEngine-facing physics backend built on NativeEngine.
//
// Mirrors the concrete-World shape MitiruEngine's backends expose (cf.
// JoltPhysicsWorld): integer BodyId handles, create/remove, pose get/set,
// force/impulse, contact callbacks -- but speaking sgc types at the boundary
// (right-handed, Y-up, quaternion (x,y,z,w) w-last, SI units), and adding the
// one thing no other MitiruEngine backend has: a NATIVE PERIODIC BOUNDARY.
//
// Build with MITIRU_HAS_NATIVEPHYS defined to get the real backend; without it,
// the Null stub below keeps header-only / no-submodule builds compiling.

#include <sgc/math/Vec3.hpp>
#include <sgc/math/Quaternion.hpp>

namespace mitiru {
namespace nativephys {

using BodyId = unsigned int;
inline constexpr BodyId kInvalidBody = 0;

}  // namespace nativephys
}  // namespace mitiru

#ifdef MITIRU_HAS_NATIVEPHYS

#include <functional>

#include "physics_world.hpp"   // ne::PhysicsWorld

namespace mitiru {
namespace nativephys {

// -- sgc <-> NativeEngine conversions (float <-> double; quat order differs) --
inline ne::V3 toNe(const sgc::Vec3f& v) { return {double(v.x), double(v.y), double(v.z)}; }
inline sgc::Vec3f toSgc(const ne::V3& v) {
    return sgc::Vec3f{float(v.x), float(v.y), float(v.z)};
}
inline ne::Q toNe(const sgc::Quaternionf& q) {
    return ne::Q{double(q.w), double(q.x), double(q.y), double(q.z)};  // (x,y,z,w) -> (w,x,y,z)
}
inline sgc::Quaternionf toSgc(const ne::Q& q) {
    return sgc::Quaternionf{float(q.x), float(q.y), float(q.z), float(q.w)};
}

class NativePhysicsWorld {
public:
    void init() {}
    void shutdown() {}

    // Configuration.
    void setGravity(const sgc::Vec3f& g) { pw_.setGravity(toNe(g)); }
    void setTimestep(float dt) { pw_.setTimestep(dt); }
    void setRestitution(float e) { pw_.setRestitution(e); }
    void setFriction(float mu) { pw_.setFriction(mu); }
    void setSleepEnabled(bool on) { pw_.setSleepEnabled(on); }
    // NATIVE PERIODIC BOUNDARY -- the distinctive capability. half = half-edge.
    void setPeriodicBox(float half, bool periodic) { pw_.setBox(half, periodic); }

    // Body creation.
    BodyId createSphere(const sgc::Vec3f& pos, float radius, float density) {
        return pw_.addSphere(toNe(pos), radius, density);
    }
    BodyId createBox(const sgc::Vec3f& pos, const sgc::Quaternionf& rot,
                     const sgc::Vec3f& halfExtents, float density) {
        return pw_.addBox(toNe(pos), toNe(rot), toNe(halfExtents), density);
    }
    BodyId createCylinder(const sgc::Vec3f& pos, const sgc::Quaternionf& rot,
                          float radius, float height, float density) {
        return pw_.addCylinder(toNe(pos), toNe(rot), radius, height, density);
    }
    BodyId makeStatic(BodyId h) { return pw_.makeStatic(h); }
    void removeBody(BodyId h) { pw_.remove(h); }

    // Stepping (mirrors JoltPhysicsWorld::update(dt)).
    void update(float dt) { pw_.setTimestep(dt); pw_.step(); }

    // Pose / velocity.
    sgc::Vec3f getPosition(BodyId h) const { return toSgc(pw_.position(h)); }
    sgc::Quaternionf getRotation(BodyId h) const { return toSgc(pw_.orientation(h)); }
    sgc::Vec3f getVelocity(BodyId h) const { return toSgc(pw_.velocity(h)); }
    void setPosition(BodyId h, const sgc::Vec3f& p) { pw_.setPosition(h, toNe(p)); }
    void setRotation(BodyId h, const sgc::Quaternionf& q) { pw_.setOrientation(h, toNe(q)); }
    void setVelocity(BodyId h, const sgc::Vec3f& v) { pw_.setVelocity(h, toNe(v)); }

    void applyForce(BodyId h, const sgc::Vec3f& f) { pw_.applyForce(h, toNe(f)); }
    void applyImpulse(BodyId h, const sgc::Vec3f& imp) { pw_.applyImpulse(h, toNe(imp)); }

    // Joints (rigid rod / spring-damper).
    void addDistanceJoint(BodyId a, BodyId b, const sgc::Vec3f& la, const sgc::Vec3f& lb,
                          float rest, float stiffness = 0.f, float damping = 0.f) {
        pw_.addDistanceJoint(a, b, toNe(la), toNe(lb), rest, stiffness, damping);
    }

    // Collision events: (a, b, worldPoint, worldNormal[b->a]).
    using ContactFn = std::function<void(BodyId, BodyId, sgc::Vec3f, sgc::Vec3f)>;
    void onContact(ContactFn cb) {
        cb_ = std::move(cb);
        pw_.setContactCallback([this](const ne::ContactInfo& ci) {
            if (cb_) cb_(ci.a, ci.b, toSgc(ci.point), toSgc(ci.normal));
        });
    }

    bool isSleeping(BodyId h) const { return pw_.isSleeping(h); }
    ne::PhysicsWorld& native() { return pw_; }   // escape hatch for advanced config

private:
    ne::PhysicsWorld pw_;
    ContactFn cb_;
};

}  // namespace nativephys
}  // namespace mitiru

#else  // !MITIRU_HAS_NATIVEPHYS -- Null stub with identical signatures

#include <functional>

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
    BodyId createSphere(const sgc::Vec3f&, float, float) { return kInvalidBody; }
    BodyId createBox(const sgc::Vec3f&, const sgc::Quaternionf&, const sgc::Vec3f&, float) { return kInvalidBody; }
    BodyId createCylinder(const sgc::Vec3f&, const sgc::Quaternionf&, float, float, float) { return kInvalidBody; }
    BodyId makeStatic(BodyId h) { return h; }
    void removeBody(BodyId) {}
    void update(float) {}
    sgc::Vec3f getPosition(BodyId) const { return {}; }
    sgc::Quaternionf getRotation(BodyId) const { return {}; }
    sgc::Vec3f getVelocity(BodyId) const { return {}; }
    void setPosition(BodyId, const sgc::Vec3f&) {}
    void setRotation(BodyId, const sgc::Quaternionf&) {}
    void setVelocity(BodyId, const sgc::Vec3f&) {}
    void applyForce(BodyId, const sgc::Vec3f&) {}
    void applyImpulse(BodyId, const sgc::Vec3f&) {}
    void addDistanceJoint(BodyId, BodyId, const sgc::Vec3f&, const sgc::Vec3f&, float, float = 0.f, float = 0.f) {}
    using ContactFn = std::function<void(BodyId, BodyId, sgc::Vec3f, sgc::Vec3f)>;
    void onContact(ContactFn) {}
    bool isSleeping(BodyId) const { return false; }
};

}  // namespace nativephys
}  // namespace mitiru

#endif  // MITIRU_HAS_NATIVEPHYS

#endif  // NATIVEENGINE_MITIRU_NATIVE_PHYSICS_WORLD_HPP
