# Embedding NativeEngine as a MitiruEngine physics backend

NativeEngine slots into MitiruEngine the same way the Jolt/Box2D bridges do: a
concrete `World` class speaking `sgc` types, a `scene::ISystem` adapter, a Null
stub selected by `#ifdef`, and a CMake submodule. It adds one capability no other
MitiruEngine backend has: **native periodic boundaries**.

## What's here

| File | Role | Standalone-testable |
|---|---|---|
| `native_physics_world.hpp` | The backend: sgc-typed `NativePhysicsWorld` mirroring `JoltPhysicsWorld`, + a Null stub under `#else` | Yes (via the sgc shim) |
| `native_physics_system.hpp` | Illustrative `scene::ISystem` adapter. The one actually in use is MitiruEngine's `include/mitiru/physics/NativePhysicsBridge.hpp` — edit that one | No (needs MitiruEngine's scene) |
| `compat/sgc/math/*.hpp` | Minimal `sgc::Vec3f`/`Quaternionf` stubs for STANDALONE TESTING only | — |
| `test_binding.cpp` | Proves the conversions round-trip and the sgc API drives a real scene | Yes |

Standalone check (does not touch MitiruEngine):

    NativeEngine\mitiru\build_binding_test.bat
    NativeEngine\build\test_binding.exe        # 12/12 PASS

In a real MitiruEngine build the `compat/` shim is NOT on the include path -- the
real `external/sgc` headers are used instead.

**Verified against the real sgc headers (2026-07-18):** the binding compiles and
runs 12/12 against MitiruEngine's actual `external/sgc/include` (`Vec3<float>` /
`Quaternion<float>`, both `(x,y,z[,w])` with matching constructors) under
`/std:c++20`. The real sgc math uses C++20 concepts, so the translation units that
include `native_physics_world.hpp` must be C++20 (MitiruEngine already is); the
`NativeEngine` static lib itself stays C++17.

## Feature surface (all sgc-typed)

The backend exposes NativeEngine's full feature set through `sgc` types:

- **Shapes:** sphere, box, cylinder, **capsule**, plane, **static triangle mesh**,
  **heightfield** (`createSphere/Box/Cylinder/Capsule/Plane/StaticMesh/Heightfield`).
- **Body types:** dynamic / static / **kinematic** (`makeStatic/makeKinematic`,
  `setKinematicTarget` for moving platforms).
- **Per-body properties:** `setMaterial` (friction/restitution), `setDamping`,
  `setGravityScale`, `setLayerMask`, `setUserData`, `setSensor`, `setCcd`.
- **Events (sgc callback with `ContactEventS`):** `onContactBegin/Stay/End`,
  `onTriggerEnter/Stay/Exit`.
- **Queries:** `raycast`, `sphereCast`, `overlapSphere` (PBC-aware, layer-filtered).
- **Joints:** `addBall/Hinge/Fixed/SliderJoint` + `setJointMotor/Limit/Breakable`.
- **Character:** `NativeCharacter` (a capsule collide-and-slide controller).
- **Integration:** `interpolatedPosition/Rotation(alpha)` for smooth rendering,
  `debugDraw(lineFn, flags)` into MitiruEngine's own renderer, `saveState/loadState`.

### The component adapter -- turning a `PhysicsTrait` into a body in one call

`BodyDesc` mirrors `PhysicsTrait` / `RigidBodyComponent`; fill it from the
component and call `createBody`. `mass > 0` overrides density to hit that mass;
`friction`/`restitution < 0` inherit the world default. This makes the schema
fields (bodyType, collider, mass, material, damping, isTrigger, collisionLayer/mask)
actually take effect at runtime -- which the un-wired Jolt path never did.

```cpp
BodyDesc d;
d.type = trait.bodyType == "kinematic" ? BodyDesc::Type::Kinematic
       : trait.bodyType == "static"    ? BodyDesc::Type::Static
                                       : BodyDesc::Type::Dynamic;
d.shape = /* map trait.colliderType */ BodyDesc::Shape::Capsule;
d.position = transform.position; d.rotation = transform.rotationQuat();
d.radius = trait.colliderSize[0]; d.halfHeight = trait.colliderSize[1];
d.mass = trait.mass; d.friction = trait.friction; d.restitution = trait.restitution;
d.isTrigger = trait.isTrigger; d.layer = trait.collisionLayer; d.mask = trait.collisionMask;
d.userData = entityId;
BodyId body = phys.createBody(d);
```

## Conventions matched (from the MitiruEngine survey)

- **Types at the boundary:** `sgc::Vec3f`, `sgc::Quaternionf`. Right-handed,
  Y-up, SI units. sgc quaternion layout is **(x, y, z, w)** with w last;
  NativeEngine's own quaternion is (w, x, y, z) -- the conversion handles this.
- **Handles:** integer `BodyId` (0 = invalid).
- **Stepping:** `update(dt)` (matches `JoltPhysicsWorld::update`).
- **TransformComponent stores Euler angles** (not a quaternion), so the ISystem
  adapter converts quaternion -> Euler on write-back (`quatToEulerYXZ`; verify
  the engine's Euler order).
- **The POD game-DLL boundary does NOT apply here:** physics is an engine-side
  C++ subsystem, so full C++/templates/sgc types are fine.

## CMake wiring (as shipped in MitiruEngine's root `CMakeLists.txt`)

Two-stage opt-in: `-DMITIRU_USE_NATIVEPHYS=ON` decides whether NativeEngine is
built at all; linking `mitiru_nativephys` decides which targets can see it.
Nothing is added to the `mitiru` target itself, so turning the option on does
not force a rebuild of the engine.

```cmake
option(MITIRU_USE_NATIVEPHYS "..." OFF)
set(MITIRU_NATIVEPHYS_ROOT "" CACHE PATH "Path to the NativeEngine repo root")
if(MITIRU_USE_NATIVEPHYS)
  # root = MITIRU_NATIVEPHYS_ROOT, else external/NativeEngine, else $ENV{NATIVEENGINE_ROOT}
  add_library(NativeEngine STATIC "${root}/src/world.cpp" "${root}/src/physics_world.cpp")
  target_include_directories(NativeEngine PUBLIC "${root}/src")
  target_compile_features(NativeEngine PUBLIC cxx_std_17)

  add_library(mitiru_nativephys INTERFACE)
  target_link_libraries(mitiru_nativephys INTERFACE NativeEngine sgc::sgc)
  target_include_directories(mitiru_nativephys INTERFACE "${root}/mitiru")
  target_compile_definitions(mitiru_nativephys INTERFACE MITIRU_HAS_NATIVEPHYS=1)
  target_compile_features(mitiru_nativephys INTERFACE cxx_std_20)
endif()
```

`mitiru_nativephys` links `sgc::sgc` so `native_physics_world.hpp` finds the real
sgc headers, and requires C++20 because that sgc math uses concepts.

## Using it from game code

```cpp
#include <native_physics_world.hpp>
using namespace mitiru::nativephys;

NativePhysicsWorld phys;
phys.setGravity({0, -9.81f, 0});
phys.setTimestep(1.f/60.f);
auto floor = phys.makeStatic(phys.createBox({0,-1,0}, {0,0,0,1}, {50,1,50}, 1));
auto ball  = phys.createSphere({0,10,0}, 0.5f, 1);
// ... each frame:
phys.update(1.f/60.f);
sgc::Vec3f p = phys.getPosition(ball);   // feed to Pipeline3D::submitMesh

// The distinctive capability -- a periodic world:
phys.setPeriodicBox(25.f, true);         // seamless wrap-around with correct
                                         // cross-boundary collisions

// A player: a capsule character controller (collide-and-slide).
NativeCharacter player(phys, 0.4f, 0.8f);
player.setPosition({0, 2, 0});
// each frame: move by input*speed*dt + gravity; player.grounded() for jump logic
player.move(inputDir * (speed * dt));

// Gameplay queries + smooth rendering:
RayHitS hit = phys.raycast(muzzle, aimDir, 100.f);       // weapons / line of sight
sgc::Vec3f smooth = phys.interpolatedPosition(ball, alpha); // judder-free at fixed step
```

## Caveats

- **Native PBC has no hook in the rest of MitiruEngine.** Its space is unbounded
  Euclidean; raycasts, the scene graph, and `TransformComponent` consumers do not
  understand wrap-around. Use `setPeriodicBox(half, true)` only when the game
  genuinely wants a toroidal world.
- **Determinism** is bit-for-bit for a given build. If the game needs
  cross-platform determinism, pin the compiler/flags.
- **Contact manifolds** cover box-box; sphere/capsule/cylinder/mesh use
  single-point contact, kept stable across frames by warm starting (fine for
  dynamics and stacks, less rigid than a full manifold for large flat rests).
- **Queries skip mesh/convex/plane** for now (`raycast` against those returns no
  hit); ray-vs-triangle-mesh is future work.
- **Compound & convex inertia** are approximations (parallel-axis diagonal / vertex
  AABB) -- fine for gameplay, off for very asymmetric bodies.
- The `native_physics_system.hpp` component field names are illustrative --
  adjust to the game's actual component; prefer the `BodyDesc` adapter above.

## Status: wired into MitiruEngine (in its test suite, 2026-07-24)

The concrete ECS adapter, the CMake wiring and the tests live IN the MitiruEngine
repo, on branch `feat/nativephys-backend`:

- `include/mitiru/physics/NativePhysicsBridge.hpp` — `physics3d::NativePhysicsSystem`
  (a `scene::ISystem`): maps entities with `RigidBodyComponent3D + TransformComponent`
  to NativeEngine bodies via `BodyDesc`, fixed-steps, writes poses back
  (quaternion -> Euler), pushes per-frame material/damping changes, follows
  kinematic/static transforms, and removes the body when the entity goes away.
- Root `CMakeLists.txt` — the two-stage opt-in above (default OFF).
- `tests/mitiru/TestNativePhysicsBridge.cpp` — one source, two targets:
  `mitiru_tests_core` (backend absent: asserts the system is a true no-op) and
  `mitiru_tests_nativephys` (backend linked: landing, determinism, periodic
  wrap-around, body removal). Run: `ctest --test-dir build -C Debug -L nativephys`.
- `docs/PHYSICS_NATIVE_BACKEND.md` — when to reach for this backend and its limits.

**Verified against the real engine, not a mock:** a real `scene::GameWorld` drops
a ball from y=5.0 to rest at y=0.489 on an AABB floor; two identical runs agree
bit-for-bit; with `periodicHalf = 5` a body at +7 m/s wraps to x = -3.03.

**Build note:** MitiruEngine sources are UTF-8 with Japanese comments and require
`/utf-8`. MitiruEngine's CMake already sets `/utf-8` on the `mitiru` target and,
being INTERFACE in header-only mode, propagates it to consumers, so the binding
TUs get it automatically. A standalone compile-check of the binding must pass
`/utf-8` (else the multibyte comments corrupt brace parsing).
