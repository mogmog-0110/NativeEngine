# NativeEngine — a deterministic rigid-body engine with native periodic boundaries

A from-scratch, dependency-free C++ rigid-body physics engine built to do the one
thing no off-the-shelf engine offers together: **hard, momentum-conserving
rigid-body contact AND a native periodic boundary**. Game engines
(PhysX/Bullet/MuJoCo/Jolt/ODE) have hard contact but cannot represent a periodic
world — one transform, one AABB, no lattice slot on constraints. MD/DEM codes
have periodic boundaries but soften the contact into a penalty potential. This
engine keeps the game-engine rigid-body method and puts the toroidal topology
into the distance function itself (minimum image), with no ghost bodies.

It is designed to embed in **MitiruEngine** as a physics backend (see
`mitiru/INTEGRATION.md`), and doubles as the simulator for a coarse-grained
macrocycle self-assembly study (the science layer is reused unchanged).

## Features

- **Shapes:** sphere, box, flat-cap cylinder — one convex path (analytic for
  sphere-sphere / sphere-cylinder, GJK + EPA otherwise; box support = one line).
- **Contact:** hard impulse with restitution, Coulomb **friction**, an iterative
  sequential-impulse solver with restitution slop; **box-box clipping manifolds**
  so boxes rest flat and **stack**.
- **Forces:** gravity; **distance joints** — rigid rod and spring-damper.
- **Boundaries:** **native periodic** (minimum image, cross-boundary collisions,
  no ghosts) and reflective walls.
- **Performance:** **sleeping** (settled bodies stop integrating, wake on
  interaction) and a PBC-aware **uniform-grid broadphase** (bit-identical to the
  O(N^2) scan).
- **Determinism:** bit-for-bit for a given build, by construction (fixed
  evaluation order) — not reconstructed via single-threading + stable sorts.
- **API:** a handle-based `PhysicsWorld` (create/remove, pose get/set, impulse,
  contact events) and a MitiruEngine-facing `NativePhysicsWorld` speaking sgc
  types.

## Layout

```
src/
  vmath.hpp        double Vec3 / quaternion, orientation integration
  body.hpp         rigid body + sphere/box/cylinder factories, bounding radius
  box.hpp          periodic/reflective box: minimum image + wrap
  gjk.hpp          support functions + GJK overlap
  epa.hpp          EPA penetration depth + contact point
  narrowphase.hpp  analytic sphere-cylinder
  manifold.hpp     box-box clipping manifold
  detect.hpp       unified narrow-phase dispatch
  contact.hpp      impulse resolver (normal + friction) + effective mass
  world.hpp/.cpp   integrator, solver, joints, broadphase, sleeping, boundaries
  physics_world.hpp/.cpp   handle-based facade + contact events
  selftest.cpp     86 physics unit tests
compat/
  PxPhysicsAPI.h   PhysX shim so the original science layer compiles unchanged
mitiru/
  native_physics_world.hpp    sgc-typed MitiruEngine backend + Null stub
  native_physics_system.hpp   scene::ISystem adapter (drop-in)
  INTEGRATION.md              CMake + usage
```

## Build & test

    NativeEngine\build.bat
    NativeEngine\build\NativeEngine.exe --selftest        # 86/86 physics tests

MitiruEngine binding (via the sgc shim) and science-layer reuse:

    NativeEngine\mitiru\build_binding_test.bat  &&  build\test_binding.exe        # 5/5
    NativeEngine\build_science_test.bat         &&  build\test_science.exe        # 9/9

Requires only MSVC (VS 2022). No external SDK.

## Status

Game-physics core + MitiruEngine binding: complete and tested. Remaining: the
macrocycle science demo (reproduce the PhysX reflective results as a correctness
proof, then run the native-PBC experiment on the wall-templating artifact),
which reuses the bonding joints and the science layer already present.
