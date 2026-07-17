// Proves the original PhysX project's science layer compiles and runs UNCHANGED
// against NativeEngine, through the PhysX-compat shim -- no PhysX SDK involved.
//
// It #includes the real bond_graph.cpp / metrics.cpp from ../../PhysxRender via
// the build command's include path (compat/ shadows PxPhysicsAPI.h), then
// exercises the graph on a hand-built target macrocycle.

#include <cstdio>

#include "bond_graph.hpp"
#include "metrics.hpp"

using namespace physx;

static int g_pass = 0, g_fail = 0;
static void check(bool ok, const char* name) {
    std::printf("  [%s] %s\n", ok ? "PASS" : "FAIL", name);
    if (ok) ++g_pass; else ++g_fail;
}

int main() {
    std::printf("=== science-layer reuse against NativeEngine (via shim) ===\n");

    // Build the target: 4 spheres + 4 linkers in an alternating 8-cycle.
    // Sphere i bonds linker i (top) and linker (i+3)%4 (bottom): a square ring.
    BondGraph g;
    for (int i = 0; i < 4; ++i) {
        g.addBond(i, i, /*isTop=*/true, nullptr);            // S_i - L_i  (top end)
        g.addBond(i, (i + 3) % 4, /*isTop=*/false, nullptr); // S_i - L_{i-1} (bottom)
    }

    auto comps = g.getComponents();
    check(comps.size() == 1, "target: single connected component");
    check(!comps.empty() && comps[0].numSpheres() == 4 && comps[0].numCylinders() == 4,
          "target: 4 spheres + 4 linkers");
    check(!comps.empty() && g.isTargetCycle(comps[0]),
          "target: recognised as the 4S-4L alternating target cycle");
    check(!comps.empty() && g.isViable(comps[0]),
          "target: viable (V1-V4)");

    // An open chain S0-L0-S1 must be viable but NOT a target.
    BondGraph chain;
    chain.addBond(0, 0, true, nullptr);
    chain.addBond(1, 0, false, nullptr);
    auto cc = chain.getComponents();
    check(cc.size() == 1 && !chain.isTargetCycle(cc[0]) && chain.isViable(cc[0]),
          "open chain: viable but not target");

    // Metrics on a scene = 1 target ring (4 S + 4 L bonded) + 8 free monomers of
    // each species. totalSpheres = 12, unreactedSpheres = 8, likewise linkers.
    PureMetrics m = computePureMetrics(g, /*totalSpheres=*/12, /*totalCylinders=*/12,
                                       /*unreactedSpheres=*/8, /*unreactedCylinders=*/8);
    check(m.N_target == 1, "metrics: N_target = 1");
    check(m.N_cycle == 1 && m.Nm[4] == 1, "metrics: one completed 4-sphere cycle");
    check(m.S4() == 1.0, "metrics: S4 = 1.0 (the only cycle is a square)");
    check(m.M_mono() == 16.0 / 24.0, "metrics: M_mono = 16/24 free monomers");

    std::printf("\n=== Summary ===\n  Passed: %d\n  Failed: %d\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
