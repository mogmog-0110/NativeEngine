# Embedding NativeEngine as a MitiruEngine physics backend

NativeEngine slots into MitiruEngine the same way the Jolt/Box2D bridges do: a
concrete `World` class speaking `sgc` types, a `scene::ISystem` adapter, a Null
stub selected by `#ifdef`, and a CMake submodule. It adds one capability no other
MitiruEngine backend has: **native periodic boundaries**.

## What's here

| File | Role | Standalone-testable |
|---|---|---|
| `native_physics_world.hpp` | The backend: sgc-typed `NativePhysicsWorld` mirroring `JoltPhysicsWorld`, + a Null stub under `#else` | Yes (via the sgc shim) |
| `native_physics_system.hpp` | `scene::ISystem` adapter syncing `TransformComponent` <-> bodies | No (needs MitiruEngine's scene) |
| `compat/sgc/math/*.hpp` | Minimal `sgc::Vec3f`/`Quaternionf` stubs for STANDALONE TESTING only | — |
| `test_binding.cpp` | Proves the conversions round-trip and the sgc API drives a real scene | Yes |

Standalone check (does not touch MitiruEngine):

    NativeEngine\mitiru\build_binding_test.bat
    NativeEngine\build\test_binding.exe        # 5/5 PASS

In a real MitiruEngine build the `compat/` shim is NOT on the include path -- the
real `external/sgc` headers are used instead.

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

## CMake wiring (mirror the Jolt/Box2D submodule pattern, ~lines 643-671)

Add NativeEngine as a submodule under `external/NativeEngine`, then in the root
`CMakeLists.txt`:

```cmake
# --- NativeEngine physics backend (native periodic boundaries) ---
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/NativeEngine/NativeEngine/src/world.cpp")
  add_library(NativeEngine STATIC
    external/NativeEngine/NativeEngine/src/world.cpp
    external/NativeEngine/NativeEngine/src/physics_world.cpp)
  target_include_directories(NativeEngine PUBLIC
    external/NativeEngine/NativeEngine/src
    external/NativeEngine/NativeEngine/mitiru)
  target_compile_features(NativeEngine PUBLIC cxx_std_17)
  target_link_libraries(mitiru ${MITIRU_TARGET_SCOPE} NativeEngine sgc)
  target_include_directories(mitiru ${MITIRU_TARGET_SCOPE}
    external/NativeEngine/NativeEngine/mitiru)
  target_compile_definitions(mitiru ${MITIRU_TARGET_SCOPE} MITIRU_HAS_NATIVEPHYS=1)
endif()
```

`NativeEngine` links `sgc` so `native_physics_world.hpp` finds the real sgc
headers. Always use `${MITIRU_TARGET_SCOPE}` (INTERFACE in header-only mode,
PUBLIC in static mode), never a hardcoded scope.

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
```

## Caveats

- **Native PBC has no hook in the rest of MitiruEngine.** Its space is unbounded
  Euclidean; raycasts, the scene graph, and `TransformComponent` consumers do not
  understand wrap-around. Use `setPeriodicBox(half, true)` only when the game
  genuinely wants a toroidal world.
- **Determinism** is bit-for-bit for a given build. If the game needs
  cross-platform determinism, pin the compiler/flags.
- **Contact manifolds** cover box-box; sphere/cylinder use single-point contact
  (fine for dynamics, and boxes stack).
- The `native_physics_system.hpp` component field names are illustrative --
  adjust to the game's actual `NativeRigidBodyComponent`.
