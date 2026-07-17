# NativeEngine roadmap — from physics core to game-engine backend

This lists what the engine still needs to be a drop-in game-physics backend (the
intent: replace MitiruEngine's physics). It is grounded in what MitiruEngine
actually asks of a backend — its `JoltPhysicsWorld` API shape, the `PhysicsTrait`
serialization schema, and the `RigidBodyComponent` fields — not a generic wishlist.

## What exists today (done)

- **Shapes:** sphere, box, cylinder (one convex path: analytic where possible,
  GJK + EPA otherwise; box-box clipping manifolds so boxes stack).
- **Dynamics:** sequential-impulse solver, Coulomb friction, restitution,
  Baumgarte position correction; angular-momentum integration (asymmetric bodies
  correct).
- **Boundaries:** native periodic (minimum image, cross-boundary collisions, no
  ghosts) + reflective walls + open ground. PBC verified (translation invariance,
  cross-boundary == free space, momentum/energy conservation — all machine-precision).
- **Joints:** distance (rigid rod / spring-damper).
- **Perf:** sleeping; PBC-aware uniform-grid broadphase (bit-identical to O(N²)).
- **Determinism:** bit-for-bit by construction (fixed evaluation order).
- **API:** handle-based `PhysicsWorld` (create/remove, pose get/set, applyForce /
  impulse / impulse-at-point, contact callback); MitiruEngine `sgc`-typed binding.
- **Tooling:** live ImGui **visual debugger** (scene hierarchy, inspector, stats,
  contact/velocity/AABB overlays, picking).

Two things here are genuinely ahead of MitiruEngine's own backends: **native PBC**
(unique) and **determinism + sleeping** (the Jolt wrapper barely exposes sleep).
Keep both as invariants when adding the features below.

## Gaps, by tier — ALL DONE ✅ (164 selftests, 5 binding tests)

### Tier 1 — game table stakes ✅
- ✅ Per-body material friction/restitution + combine modes (#20)
- ✅ Linear/angular damping (#21)
- ✅ Collision layers + masks (#23)
- ✅ Body userData → ECS entity (#24)
- ✅ Kinematic bodies (script-driven, imparts velocity) (#22)
- ✅ Triggers/sensors + enter/stay/exit (#25)
- ✅ Contact events begin/stay/end (#26)

### Tier 2 — characters & queries ✅
- ✅ Capsule shape (#27)
- ✅ Scene queries: raycast / overlap / sphere-cast (#28)

### Tier 3 — world geometry & articulation ✅
- ✅ Static triangle mesh + heightfield + BVH midphase (#29)
- ✅ Convex hull + plane (#30); compound (#36)
- ✅ Joints: ball / hinge / fixed / slider (+ limits, motors, breakable) (#31)

### Tier 4 — advanced motion ✅
- ✅ Kinematic capsule character controller (collide-and-slide, grounded, slopes) (#32)
- ✅ Continuous collision detection (#33)

### Tier 5 — MitiruEngine integration ✅
- ✅ Debug-draw hook (DebugDrawFlags), render interpolation, per-body gravity,
  full force API (torque, force-at-point, angular impulse) (#34)

### Tier 6 — solver robustness & scale ✅
- ✅ Warm-started contacts, island sleeping, snapshot save/load (#35)

## Remaining / future (small, non-blocking)

- ray-vs-mesh / ray-vs-convex / ray-vs-plane (queries currently skip these shapes)
- exact convex & compound inertia (currently AABB / parallel-axis-diagonal
  approximations — fine for gameplay, off for very asymmetric bodies)
- box-box-style multi-point manifolds for capsule/convex/mesh resting contacts
  (single-point today; stable via warm starting but not as rigid as a manifold)
- optional deterministic-order multithreaded solver (single-thread is the
  deterministic default)

The engine now covers the full game-physics feature surface — shapes, queries,
joints, character controller, CCD, events, filtering, materials, debug draw,
interpolation, and snapshotting — while keeping its two differentiators (native
PBC, determinism) intact.
