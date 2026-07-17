#pragma once
#ifndef NATIVEENGINE_VIEWER_GLMATH_HPP
#define NATIVEENGINE_VIEWER_GLMATH_HPP

// Minimal float vector / column-major 4x4 matrix math for the OpenGL viewer.
// Right-handed, Y-up -- matching NativeEngine, so the render is not mirrored.
#include <cmath>

namespace nv {

struct Vec3 {
    float x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(float a, float b, float c) : x(a), y(b), z(c) {}
    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(float s) const { return {x * s, y * s, z * s}; }
    float dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    Vec3 norm() const { float n = std::sqrt(dot(*this)); return n > 1e-8f ? *this * (1 / n) : *this; }
};

// Column-major 4x4 (OpenGL layout): m[col*4 + row].
struct Mat4 {
    float m[16] = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1};

    static Mat4 identity() { return {}; }

    Mat4 operator*(const Mat4& b) const {
        Mat4 r;
        for (int c = 0; c < 4; ++c)
            for (int row = 0; row < 4; ++row) {
                float s = 0;
                for (int k = 0; k < 4; ++k) s += m[k * 4 + row] * b.m[c * 4 + k];
                r.m[c * 4 + row] = s;
            }
        return r;
    }

    static Mat4 translate(const Vec3& t) {
        Mat4 r; r.m[12] = t.x; r.m[13] = t.y; r.m[14] = t.z; return r;
    }
    static Mat4 scale(const Vec3& s) {
        Mat4 r; r.m[0] = s.x; r.m[5] = s.y; r.m[10] = s.z; return r;
    }
    // Rotation from a unit quaternion (w, x, y, z).
    static Mat4 fromQuat(float w, float x, float y, float z) {
        Mat4 r;
        r.m[0] = 1 - 2 * (y * y + z * z); r.m[1] = 2 * (x * y + z * w);     r.m[2] = 2 * (x * z - y * w);
        r.m[4] = 2 * (x * y - z * w);     r.m[5] = 1 - 2 * (x * x + z * z); r.m[6] = 2 * (y * z + x * w);
        r.m[8] = 2 * (x * z + y * w);     r.m[9] = 2 * (y * z - x * w);     r.m[10] = 1 - 2 * (x * x + y * y);
        return r;
    }

    static Mat4 perspective(float fovyRad, float aspect, float znear, float zfar) {
        float f = 1.0f / std::tan(fovyRad * 0.5f);
        Mat4 r; for (int i = 0; i < 16; ++i) r.m[i] = 0;
        r.m[0] = f / aspect; r.m[5] = f;
        r.m[10] = (zfar + znear) / (znear - zfar);
        r.m[11] = -1.0f;
        r.m[14] = (2 * zfar * znear) / (znear - zfar);
        return r;
    }

    static Mat4 lookAt(const Vec3& eye, const Vec3& center, const Vec3& up) {
        Vec3 f = (center - eye).norm();
        Vec3 s = f.cross(up).norm();
        Vec3 u = s.cross(f);
        Mat4 r;
        r.m[0] = s.x; r.m[4] = s.y; r.m[8]  = s.z;
        r.m[1] = u.x; r.m[5] = u.y; r.m[9]  = u.z;
        r.m[2] = -f.x; r.m[6] = -f.y; r.m[10] = -f.z;
        r.m[12] = -s.dot(eye); r.m[13] = -u.dot(eye); r.m[14] = f.dot(eye);
        return r;
    }
};

}  // namespace nv

#endif  // NATIVEENGINE_VIEWER_GLMATH_HPP
