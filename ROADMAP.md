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

## Gaps, by tier (→ task IDs)

### Tier 1 — game table stakes (MitiruEngine's schema/components require these)
The `PhysicsTrait`/`RigidBodyComponent` fields that today have no runtime effect:
- ~~**Per-body material** friction/restitution + combine modes~~ ✅ done (#20)
- ~~**Linear/angular damping** per body~~ ✅ done (#21)
- ~~**Collision layers + masks**~~ ✅ done (#23)
- ~~**Body userData** (body → ECS entity)~~ ✅ done (#24)
- **Kinematic bodies** (script-driven, infinite mass, imparts velocity) — only
  static/dynamic exist → **#22**
- **Triggers/sensors** + enter/stay/exit — `isTrigger` has no real semantics → **#25**
- **Contact events begin/stay/end** — callback is per-frame, no transitions → **#26**

### Tier 2 — characters & queries
- **Capsule shape** — the character/prop shape; in schema + both components → **#27**
- **Scene queries** raycast / overlap / shapecast — gameplay staple; only a
  viewer-side ray-vs-sphere exists → **#28**

### Tier 3 — world geometry & articulation
- **Static triangle mesh + heightfield** (level/terrain collision) + BVH midphase
  — the biggest hole for real levels → **#29**
- **Convex hull + compound + plane** shapes → **#30**
- **Joints** hinge / ball-socket / fixed / slider (+ limits, motors, breakable) —
  only distance exists → **#31**

### Tier 4 — advanced motion
- **Kinematic character controller** (capsule, slopes, step-up, collide-and-slide)
  — absent everywhere; needs capsule + shapecast → **#32**
- **Continuous collision detection** for fast bodies (no tunneling) → **#33**

### Tier 5 — MitiruEngine integration polish
- **Debug-draw hook** feeding MitiruEngine's renderer (DebugDrawFlags) instead of
  only the standalone viewer; **render interpolation** (prev/curr + factor);
  **per-body gravity factor**; round out the force API (torque, force-at-point) → **#34**

### Tier 6 — solver robustness & scale (keep determinism as default)
- **Warm-started persistent manifolds** (stable stacks, fewer iterations),
  **simulation islands**, **world snapshot serialization**; optional
  deterministic-order parallel solver later → **#35**

## Suggested order

Tier 1 is a cohesive low-risk batch (all Body + facade changes, no new
narrow-phase) that immediately makes the `PhysicsTrait` schema real — do it first.
Then capsule (#27) + queries (#28) unlock the character controller (#32). Mesh
colliders (#29) unlock real levels. Joints (#31) unlock ragdolls/vehicles. Tiers
5–6 are integration and scale, valuable once the feature surface is complete.

Not on the critical path but worth stating: MitiruEngine currently wires **no**
physics into a shipping game (usage is unit tests + the un-wired Jolt submodule),
so this backend can define the API rather than chase call sites.
