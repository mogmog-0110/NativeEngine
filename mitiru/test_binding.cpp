// Standalone test of the MitiruEngine binding, using the minimal sgc shim on the
// include path. Proves the sgc<->NativeEngine conversions round-trip and that the
// backend drives a real scene through the sgc-typed API. A real MitiruEngine
// build uses the actual sgc headers instead of the shim.
#define MITIRU_HAS_NATIVEPHYS 1

#include <cmath>
#include <cstdio>

#include "native_physics_world.hpp"

using namespace mitiru::nativephys;

static int g_pass = 0, g_fail = 0;
static void check(bool ok, const char* name) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (ok) ++g_pass; else ++g_fail;
}

int main() {
    std::printf("=== NativeEngine <-> MitiruEngine binding test ===\n");

    // Quaternion order round-trip: sgc (x,y,z,w) <-> ne (w,x,y,z).
    sgc::Quaternionf q{0.1f, 0.2f, 0.3f, 0.9f};
    sgc::Quaternionf q2 = toSgc(toNe(q));
    check(std::fabs(q2.x - 0.1f) < 1e-6 && std::fabs(q2.y - 0.2f) < 1e-6 &&
              std::fabs(q2.z - 0.3f) < 1e-6 && std::fabs(q2.w - 0.9f) < 1e-6,
          "binding: quaternion (x,y,z,w) order round-trips");

    sgc::Vec3f v{1.5f, -2.5f, 3.5f};
    sgc::Vec3f v2 = toSgc(toNe(v));
    check(v2.x == 1.5f && v2.y == -2.5f && v2.z == 3.5f, "binding: vec3 round-trips");

    // Drive a real scene through the sgc-typed API: a sphere falls onto a floor.
    NativePhysicsWorld w;
    w.init();
    w.setGravity(sgc::Vec3f{0, -10, 0});
    w.setTimestep(1.0f / 240.0f);
    w.setRestitution(0.0f);
    w.setSleepEnabled(true);
    BodyId floor = w.makeStatic(
        w.createBox(sgc::Vec3f{0, -1, 0}, sgc::Quaternionf{0, 0, 0, 1},
                    sgc::Vec3f{20, 1, 20}, 1.0f));
    BodyId ball = w.createSphere(sgc::Vec3f{0, 5, 0}, 1.0f, 1.0f);
    (void)floor;

    for (int s = 0; s < 3000; ++s) w.update(1.0f / 240.0f);
    sgc::Vec3f p = w.getPosition(ball);
    check(p.y > 0.9f && p.y < 1.2f, "binding: sphere falls and rests on the floor (sgc API)");
    check(w.isSleeping(ball), "binding: settled body sleeps");

    // The distinctive capability: a native periodic box.
    NativePhysicsWorld pw;
    pw.setPeriodicBox(5.0f, true);
    pw.setTimestep(1.0f / 60.0f);
    BodyId a = pw.createSphere(sgc::Vec3f{4.5f, 0, 0}, 0.5f, 1.0f);
    pw.setVelocity(a, sgc::Vec3f{2, 0, 0});
    for (int s = 0; s < 100; ++s) pw.update(1.0f / 60.0f);
    sgc::Vec3f pa = pw.getPosition(a);
    check(std::fabs(pa.x) <= 5.0f + 1e-4f, "binding: native periodic box wraps a body (sgc API)");

    std::printf("\n=== Summary ===\n  Passed: %d\n  Failed: %d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
