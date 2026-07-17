# NativeEngine — a deterministic rigid-body engine with native periodic boundaries

A from-scratch, dependency-free C++ rigid-body physics engine built to do the one
thing no off-the-shelf engine offers together: **hard, momentum-conserving
rigid-body contact AND a native periodic boundary**. Game engines
(PhysX/Bullet/MuJoCo/Jolt/ODE) have hard contact but cannot represent a periodic
world — one transform, one AABB, no lattice slot on constraints. MD/DEM codes
have periodic boundaries but soften the contact into a penalty potential. This
engine keeps the game-engine rigid-body method and puts the toroidal topology
into the distance function itself (minimum image), with no ghost bodies.

It is a standalone project (no dependency on any other repo). It is designed to
embed in **MitiruEngine** as a physics backend (see `mitiru/INTEGRATION.md`), and
its science-facing use (a coarse-grained macrocycle self-assembly study) lives in
the consuming project as a bridge — the engine itself knows nothing about it.

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
viewer/
  viewer.cpp       LIVE interactive demo viewer (steps the engine in real time)
  gl_core.hpp      OpenGL 3.3 core loader; glmath.hpp   float math for rendering
  build_viewer.bat freeglut + modern GL (instanced Blinn-Phong, MSAA)
mitiru/
  native_physics_world.hpp    sgc-typed MitiruEngine backend + Null stub
  native_physics_system.hpp   scene::ISystem adapter (drop-in)
  INTEGRATION.md              CMake + usage
```

## Build & test

    build.bat
    build\NativeEngine.exe --selftest        # 86/86 physics tests

MitiruEngine binding (via the sgc shim):

    mitiru\build_binding_test.bat  &&  build\test_binding.exe        # 5/5

Core + selftest + binding require only MSVC (VS 2022) — no external SDK.

## Live demos

The viewer steps the engine in real time (not a recording) and is interactive —
classic PhysX-snippet scenes on an open ground, plus the periodic-boundary
feature demos. Requires freeglut (`FREEGLUT_ROOT`, default `E:\dev\freeglut`).

    viewer\build_viewer.bat
    build\NativeViewer.exe 1        # 1 pyramid  2 brick-wall+ball  3 dominoes
                                    # 4 sphere pile  5 cylinder pile  6 rope
                                    # 7 elastic gas (box)  8 PERIODIC gas  9 PERIODIC pair

Keys: `1`-`9` scene · `n` next · `r` reset · `d` drop a body · `space` pause ·
`,`/`.` sim speed · `h` recenter camera · drag to orbit · wheel to zoom.

## Status

Game-physics core + MitiruEngine binding + live interactive viewer: complete and
tested. The macrocycle science use (reproduce the PhysX reflective results, then
run the native-PBC experiment) lives in the consuming project as a bridge, and
reuses the bonding joints and science layer.
