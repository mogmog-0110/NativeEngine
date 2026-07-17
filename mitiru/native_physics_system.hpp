#pragma once
#ifndef NATIVEENGINE_MITIRU_NATIVE_PHYSICS_SYSTEM_HPP
#define NATIVEENGINE_MITIRU_NATIVE_PHYSICS_SYSTEM_HPP

// ECS adapter that drives NativePhysicsWorld from MitiruEngine's scene.
//
// DROP-IN for MitiruEngine (references its scene/component types, so it is NOT
// compiled in NativeEngine's standalone build -- only inside MitiruEngine).
// Modelled on the existing PhysicsSystem3D: it owns EntityId<->BodyId maps,
// creates bodies from a component, runs the fixed-step world, and writes poses
// back to TransformComponent (converting the quaternion to the Euler angles
// TransformComponent stores).
//
// Requires, from MitiruEngine:
//   scene::ISystem, scene::GameWorld, scene::EntityId, scene::TransformComponent
//   (position: sgc::Vec3f, rotation: sgc::Vec3f Euler radians, scale: sgc::Vec3f)
// and a component you define, e.g.:
//   struct NativeRigidBodyComponent {
//       enum Shape { Sphere, Box, Cylinder } shape = Sphere;
//       sgc::Vec3f halfExtents{0.5f,0.5f,0.5f};  // box; or radius in .x
//       float radius = 0.5f, height = 1.0f, density = 1.0f;
//       bool isStatic = false;
//   };
//
// This header is illustrative scaffolding: adjust the includes / component field
// names to match the game's actual definitions.

#if defined(MITIRU_HAS_NATIVEPHYS) && defined(MITIRU_ENGINE_PRESENT)

#include <cmath>
#include <unordered_map>

#include <mitiru/scene/GameWorld.hpp>
#include <mitiru/scene/SystemRunner.hpp>

#include "native_physics_world.hpp"

namespace mitiru {
namespace nativephys {

// Quaternion (x,y,z,w) -> intrinsic YXZ Euler angles, matching the engine's
// TransformComponent rotation convention. Adjust if the engine differs.
inline sgc::Vec3f quatToEulerYXZ(const sgc::Quaternionf& q) {
    const float x = q.x, y = q.y, z = q.z, w = q.w;
    float sinp = 2.f * (w * x - y * z);
    float px = std::fabs(sinp) >= 1.f ? std::copysign(3.14159265f / 2.f, sinp) : std::asin(sinp);
    float py = std::atan2(2.f * (w * y + x * z), 1.f - 2.f * (x * x + y * y));
    float pz = std::atan2(2.f * (w * z + x * y), 1.f - 2.f * (x * x + z * z));
    return sgc::Vec3f{px, py, pz};
}

class NativePhysicsSystem : public scene::ISystem {
public:
    std::string name() const override { return "NativePhysicsSystem"; }

    void update(scene::GameWorld& world, float dt) override {
        syncToPhysics(world);
        accumulator_ += dt;
        const float h = world_.native().timestep() > 0 ? float(world_.native().timestep()) : 1.f / 60.f;
        int steps = 0;
        while (accumulator_ >= h && steps < maxSubsteps_) { world_.update(h); accumulator_ -= h; ++steps; }
        syncFromPhysics(world);
    }

    NativePhysicsWorld& world() { return world_; }

private:
    NativePhysicsWorld world_;
    std::unordered_map<scene::EntityId, BodyId> toBody_;
    float accumulator_ = 0.f;
    int maxSubsteps_ = 5;

    void syncToPhysics(scene::GameWorld& world) {
        world.template forEach<NativeRigidBodyComponent>(
            [&](scene::EntityId e, NativeRigidBodyComponent& rb) {
                if (toBody_.count(e)) return;                 // already created
                auto* tr = world.template getComponent<scene::TransformComponent>(e);
                if (!tr) return;
                sgc::Quaternionf rot{};                        // identity; Euler->quat if needed
                BodyId h = kInvalidBody;
                switch (rb.shape) {
                    case NativeRigidBodyComponent::Sphere:
                        h = world_.createSphere(tr->position, rb.radius, rb.density); break;
                    case NativeRigidBodyComponent::Box:
                        h = world_.createBox(tr->position, rot, rb.halfExtents, rb.density); break;
                    case NativeRigidBodyComponent::Cylinder:
                        h = world_.createCylinder(tr->position, rot, rb.radius, rb.height, rb.density); break;
                }
                if (rb.isStatic) world_.makeStatic(h);
                toBody_[e] = h;
            });
    }

    void syncFromPhysics(scene::GameWorld& world) {
        for (auto& kv : toBody_) {
            auto* tr = world.template getComponent<scene::TransformComponent>(kv.first);
            if (!tr) continue;
            tr->position = world_.getPosition(kv.second);
            tr->rotation = quatToEulerYXZ(world_.getRotation(kv.second));  // TransformComponent stores Euler
        }
    }
};

}  // namespace nativephys
}  // namespace mitiru

#endif  // MITIRU_HAS_NATIVEPHYS && MITIRU_ENGINE_PRESENT
#endif  // NATIVEENGINE_MITIRU_NATIVE_PHYSICS_SYSTEM_HPP
