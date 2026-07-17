#pragma once
// Minimal sgc::Vec3 stub for STANDALONE TESTING of the MitiruEngine binding.
// In a real MitiruEngine build this file is shadowed by the actual
// external/sgc/include/sgc/math/Vec3.hpp (which has the full API). The binding
// only touches x/y/z, so this suffices to compile and test the conversions.
namespace sgc {
template <class T>
struct Vec3 {
    T x{}, y{}, z{};
    Vec3() = default;
    Vec3(T a, T b, T c) : x(a), y(b), z(c) {}
};
using Vec3f = Vec3<float>;
using Vec3d = Vec3<double>;
}  // namespace sgc
