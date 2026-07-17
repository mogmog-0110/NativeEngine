# NativeEngine — a native-PBC rigid-body engine for macrocycle assembly

A from-scratch rigid-body engine built to do the one thing no off-the-shelf
engine offers together: **hard, momentum-conserving rigid-body contact AND a
native periodic boundary**. Game engines (PhysX/Bullet/MuJoCo/Jolt/ODE) have hard
contact but cannot represent a periodic world (one transform, one AABB, no
lattice slot on constraints); MD/DEM codes have periodic boundaries but soften
the contact into a penalty potential. This engine keeps the rigid-body method and
puts the toroidal topology into the distance function itself.

It also makes the project's two hard-won properties trivial rather than delicate:
- **Determinism** is bit-for-bit by construction (fixed evaluation order, we own
  every reduction), not something reconstructed via single-threading + stable
  sorts as in the PhysX build.
- The **science layer** (bond graph, V1–V4 viability, S4, MSER, its unit tests)
  is engine-agnostic and is reused unchanged (increment 2).

## Status

Increment 1 (done): deterministic rigid-body core + native PBC + sphere contact.
- `vmath.hpp`  double-precision Vec3 / quaternion, orientation integration
- `body.hpp`   rigid body + sphere/cylinder inertia factories
- `box.hpp`    periodic/reflective box: minimum image + wrap
- `world.hpp/.cpp`  symplectic integrator, sphere–sphere impulse contact, walls
- `selftest.cpp`    18 physics tests (conservation, PBC transit, determinism)

Increments 2–4 (planned): reuse the science layer via a PhysX-compat shim; add
the cylinder body + narrow phase + bonds; validate against PhysX reflective, then
run true-PBC experiments to test whether the wall-templating artifact disappears.

## Build & test

    NativeEngine\build.bat
    NativeEngine\build\NativeEngine.exe --selftest

Requires only MSVC (VS 2022). No external SDK.
