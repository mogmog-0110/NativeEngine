#pragma once
#ifndef NATIVEENGINE_PHYSX_COMPAT_H
#define NATIVEENGINE_PHYSX_COMPAT_H

// -----------------------------------------------------------------------------
// PhysX-compatibility shim.
//
// The science layer of the original project -- bond_graph.{hpp,cpp} (connected
// components, cycle/target/viability V1-V4) and metrics.{hpp,cpp} (Phi_target,
// S4, MSER) -- is engine-agnostic in substance but was written against a handful
// of PhysX symbols. An audit shows it uses ONLY:
//     physx::PxVec3   (with .cross/.dot/.magnitude and arithmetic operators)
//     physx::PxJoint  (held as an opaque pointer; PxD6Joint/PxFixedJoint are
//                      just documented subtypes, never dereferenced)
//     physx::PxPi, physx::PxClamp
//
// Providing those here, on NativeEngine's include path AHEAD of the real SDK,
// lets bond_graph.cpp and metrics.cpp compile UNCHANGED against this engine. No
// PhysX SDK, no rigid-body types -- just enough to satisfy the science layer.
//
// This is the whole point of the "science layer is portable" claim, made
// concrete: the 2000-odd lines of graph/metrics/viability logic and their unit
// tests move to any engine behind a ~60-line shim.
// -----------------------------------------------------------------------------

#include <cmath>

namespace physx {

using PxReal = float;

constexpr PxReal PxPi = 3.14159265358979323846f;

struct PxVec3 {
    float x = 0.f, y = 0.f, z = 0.f;

    PxVec3() = default;
    PxVec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
    explicit PxVec3(float a) : x(a), y(a), z(a) {}

    PxVec3 operator+(const PxVec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    PxVec3 operator-(const PxVec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    PxVec3 operator-() const { return {-x, -y, -z}; }
    PxVec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    PxVec3 operator/(float s) const { return {x / s, y / s, z / s}; }
    PxVec3& operator+=(const PxVec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    PxVec3& operator-=(const PxVec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }

    float dot(const PxVec3& o) const { return x * o.x + y * o.y + z * o.z; }
    PxVec3 cross(const PxVec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    float magnitudeSquared() const { return x * x + y * y + z * z; }
    float magnitude() const { return std::sqrt(magnitudeSquared()); }
    PxVec3 getNormalized() const {
        float m = magnitude();
        return (m > 1e-30f) ? PxVec3{x / m, y / m, z / m} : PxVec3{0.f, 0.f, 0.f};
    }
    void normalize() {
        float m = magnitude();
        if (m > 1e-30f) { x /= m; y /= m; z /= m; }
    }
};

inline PxVec3 operator*(float s, const PxVec3& v) { return v * s; }

template <typename T>
inline T PxClamp(T v, T lo, T hi) { return v < lo ? lo : (v > hi ? hi : v); }

inline float PxAbs(float v) { return std::fabs(v); }

// Bonds are stored by the graph as opaque handles only -- never dereferenced --
// so incomplete forward-declared types are sufficient.
class PxJoint;
class PxD6Joint;
class PxFixedJoint;

}  // namespace physx

#endif  // NATIVEENGINE_PHYSX_COMPAT_H
