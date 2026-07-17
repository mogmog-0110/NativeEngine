#pragma once
#ifndef NATIVEENGINE_VMATH_HPP
#define NATIVEENGINE_VMATH_HPP

// Minimal double-precision vector / quaternion math for the native engine.
//
// Double precision is deliberate: the engine's job is to reproduce PhysX's
// reflective results STATISTICALLY and to add a native periodic boundary, not to
// match PhysX bit-for-bit. Determinism here means "same input -> same output on
// this build", which double + a fixed evaluation order gives us for free.

#include <cmath>

namespace ne {

constexpr double PI = 3.14159265358979323846;

struct V3 {
    double x = 0.0, y = 0.0, z = 0.0;

    V3() = default;
    V3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    V3 operator+(const V3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    V3 operator-(const V3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    V3 operator-() const { return {-x, -y, -z}; }
    V3 operator*(double s) const { return {x * s, y * s, z * s}; }
    V3 operator/(double s) const { return {x / s, y / s, z / s}; }
    V3& operator+=(const V3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    V3& operator-=(const V3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }

    double dot(const V3& o) const { return x * o.x + y * o.y + z * o.z; }
    V3 cross(const V3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double norm2() const { return x * x + y * y + z * z; }
    double norm() const { return std::sqrt(norm2()); }
    V3 normalized() const {
        double n = norm();
        return (n > 1e-300) ? V3{x / n, y / n, z / n} : V3{0, 0, 0};
    }
};

inline V3 operator*(double s, const V3& v) { return v * s; }

// Quaternion (w, x, y, z), rotates body -> world.
struct Q {
    double w = 1.0, x = 0.0, y = 0.0, z = 0.0;

    Q() = default;
    Q(double w_, double x_, double y_, double z_) : w(w_), x(x_), y(y_), z(z_) {}

    Q operator*(const Q& o) const {
        return {
            w * o.w - x * o.x - y * o.y - z * o.z,
            w * o.x + x * o.w + y * o.z - z * o.y,
            w * o.y - x * o.z + y * o.w + z * o.x,
            w * o.z + x * o.y - y * o.x + z * o.w};
    }
    Q operator+(const Q& o) const { return {w + o.w, x + o.x, y + o.y, z + o.z}; }
    Q operator*(double s) const { return {w * s, x * s, y * s, z * s}; }

    double norm() const { return std::sqrt(w * w + x * x + y * y + z * z); }
    Q normalized() const {
        double n = norm();
        return (n > 1e-300) ? Q{w / n, x / n, y / n, z / n} : Q{1, 0, 0, 0};
    }
    Q conjugate() const { return {w, -x, -y, -z}; }

    // Rotate a world/body vector by this quaternion (v' = q v q*).
    V3 rotate(const V3& v) const {
        // t = 2 * (q_vec x v); v' = v + w*t + q_vec x t
        V3 qv{x, y, z};
        V3 t = 2.0 * qv.cross(v);
        return v + w * t + qv.cross(t);
    }
    V3 inverseRotate(const V3& v) const { return conjugate().rotate(v); }
};

// Integrate an orientation by a world-frame angular velocity for time dt,
// then renormalise. dq/dt = 0.5 * omega_quat * q.
inline Q integrateOrientation(const Q& q, const V3& omega, double dt) {
    Q omegaQ{0.0, omega.x, omega.y, omega.z};
    Q dq = (omegaQ * q) * (0.5 * dt);
    return (q + dq).normalized();
}

}  // namespace ne

#endif  // NATIVEENGINE_VMATH_HPP
