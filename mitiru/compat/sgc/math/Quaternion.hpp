#pragma once
// Minimal sgc::Quaternion stub for STANDALONE TESTING of the MitiruEngine
// binding. Real builds use external/sgc/include/sgc/math/Quaternion.hpp.
// CRITICAL: sgc's layout is (x, y, z, w) with w LAST and w=1 the identity --
// the binding's conversion depends on this order (NativeEngine's own quaternion
// is (w, x, y, z)).
namespace sgc {
template <class T>
struct Quaternion {
    T x{}, y{}, z{}, w{1};
    Quaternion() = default;
    Quaternion(T x_, T y_, T z_, T w_) : x(x_), y(y_), z(z_), w(w_) {}
};
using Quaternionf = Quaternion<float>;
}  // namespace sgc
