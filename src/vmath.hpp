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

// Row-major 3x3 matrix, for constraint (joint) effective-mass solves.
struct Mat3 {
    double m[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};   // m[row*3 + col]

    static Mat3 zero() { Mat3 r; for (double& v : r.m) v = 0.0; return r; }
    static Mat3 identity() { return Mat3{}; }
    static Mat3 scaled(double s) { Mat3 r = zero(); r.m[0] = r.m[4] = r.m[8] = s; return r; }
    static Mat3 diagonal(const V3& d) { Mat3 r = zero(); r.m[0] = d.x; r.m[4] = d.y; r.m[8] = d.z; return r; }
    // Skew (cross-product) matrix: skew(v) * x == v.cross(x).
    static Mat3 skew(const V3& v) {
        return Mat3{ 0, -v.z, v.y,  v.z, 0, -v.x,  -v.y, v.x, 0 };
    }
    Mat3() = default;
    Mat3(double a, double b, double c, double d, double e, double f, double g, double h, double i)
        : m{a, b, c, d, e, f, g, h, i} {}

    V3 operator*(const V3& v) const {
        return {m[0] * v.x + m[1] * v.y + m[2] * v.z,
                m[3] * v.x + m[4] * v.y + m[5] * v.z,
                m[6] * v.x + m[7] * v.y + m[8] * v.z};
    }
    Mat3 operator*(const Mat3& o) const {
        Mat3 r = zero();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                for (int k = 0; k < 3; ++k) r.m[i * 3 + j] += m[i * 3 + k] * o.m[k * 3 + j];
        return r;
    }
    Mat3 operator+(const Mat3& o) const { Mat3 r; for (int i = 0; i < 9; ++i) r.m[i] = m[i] + o.m[i]; return r; }
    Mat3 operator-(const Mat3& o) const { Mat3 r; for (int i = 0; i < 9; ++i) r.m[i] = m[i] - o.m[i]; return r; }
    Mat3 transposed() const { return Mat3{m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]}; }

    Mat3 inverse() const {
        double a = m[0], b = m[1], c = m[2], d = m[3], e = m[4], f = m[5], g = m[6], h = m[7], i = m[8];
        double A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
        double det = a * A + b * B + c * C;
        if (std::fabs(det) < 1e-18) return zero();
        double inv = 1.0 / det;
        return Mat3{A * inv, (c * h - b * i) * inv, (b * f - c * e) * inv,
                    B * inv, (a * i - c * g) * inv, (c * d - a * f) * inv,
                    C * inv, (b * g - a * h) * inv, (a * e - b * d) * inv};
    }
    // Solve M x = rhs (returns zero on a singular matrix).
    V3 solve(const V3& rhs) const { return inverse() * rhs; }
};

// Integrate an orientation by a world-frame angular velocity for time dt,
// then renormalise. dq/dt = 0.5 * omega_quat * q.
inline Q integrateOrientation(const Q& q, const V3& omega, double dt) {
    Q omegaQ{0.0, omega.x, omega.y, omega.z};
    Q dq = (omegaQ * q) * (0.5 * dt);
    return (q + dq).normalized();
}

// The world-frame angular velocity that rotates `from` to `to` over dt (inverse
// of integrateOrientation, used to drive kinematic bodies to a target pose).
inline V3 angularVelocityBetween(const Q& from, const Q& to, double dt) {
    if (dt <= 0.0) return {0, 0, 0};
    Q rel = to * from.conjugate();          // world-frame relative rotation
    if (rel.w < 0.0) rel = rel * -1.0;      // shortest arc
    V3 axis{rel.x, rel.y, rel.z};
    double s = axis.norm();
    if (s < 1e-12) return {0, 0, 0};
    double angle = 2.0 * std::atan2(s, rel.w);
    return axis * (angle / (s * dt));
}

}  // namespace ne

#endif  // NATIVEENGINE_VMATH_HPP
