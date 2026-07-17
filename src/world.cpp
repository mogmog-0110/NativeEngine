#include "world.hpp"

#include "contact.hpp"
#include "narrowphase.hpp"
#include "epa.hpp"
#include "detect.hpp"
#include "manifold.hpp"
#include "queries.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>

namespace ne {

// A one-point contact constraint for the iterative solver. Single point per pair
// for now (box-box face manifolds via clipping come next); the accumulated
// impulses + restitution bias make it a proper sequential-impulse solver that
// handles several simultaneous contacts on one body.
namespace {
struct Constraint {
    std::size_t i, j;
    V3 n;              // b -> a
    V3 rA, rB;         // lever arms (minimum-imaged)
    double overlap;
    double restBias;   // target separating speed from restitution
    double mu;         // combined Coulomb friction for this contact
    double nImp = 0.0, tImp = 0.0;
    std::uint64_t pairKey = 0;   // warm-start key (body-index pair)
    int slot = 0;                // contact index within the pair
};

inline std::uint64_t pairKeyOf(std::size_t i, std::size_t j) {
    return ((std::uint64_t)i << 32) | (std::uint32_t)j;
}

// Below this approach speed, restitution is suppressed so resting contacts do
// not jitter or gain energy.
constexpr double kRestitutionSlop = 0.5;

double combineMat(World::Combine m, double a, double b) {
    switch (m) {
        case World::Combine::Average: return 0.5 * (a + b);
        case World::Combine::Min: return a < b ? a : b;
        case World::Combine::Max: return a > b ? a : b;
        case World::Combine::Multiply: return a * b;
        case World::Combine::GeometricMean: return std::sqrt((a > 0 ? a : 0) * (b > 0 ? b : 0));
    }
    return 0.5 * (a + b);
}

// Contacts between `other` (any non-mesh body) and a static triangle mesh, one
// per overlapping triangle. Each triangle is treated as a thin prism (a valid 3D
// convex) so the existing GJK/EPA path is robust; normals point mesh -> other.
void meshContacts(const Body& other, const Body& meshBody, const Box& box,
                  std::vector<Contact>& out) {
    if (!meshBody.mesh) return;
    const MeshData& md = *meshBody.mesh;
    V3 rel = box.minImage(other.x - meshBody.x);
    V3 cLocal = meshBody.q.inverseRotate(rel);
    double rr = other.boundingRadius() + 0.05;
    std::vector<int> tris;
    md.query(cLocal - V3{rr, rr, rr}, cLocal + V3{rr, rr, rr}, tris);
    const double thick = 1e-3;
    Body tri;                                  // reused temp convex (thin prism)
    tri.shape = Shape::Convex; tri.dynamic = false;
    tri.invMass = 0.0; tri.invInertiaBody = {0, 0, 0};
    tri.vertices.resize(6);
    for (int ti : tris) {
        V3 w0 = meshBody.x + meshBody.q.rotate(md.v(ti, 0));
        V3 w1 = meshBody.x + meshBody.q.rotate(md.v(ti, 1));
        V3 w2 = meshBody.x + meshBody.q.rotate(md.v(ti, 2));
        V3 n = (w1 - w0).cross(w2 - w0);
        double nl = n.norm();
        if (nl < 1e-12) continue;              // degenerate triangle
        n = n / nl;
        tri.vertices[0] = w0 + n * thick; tri.vertices[1] = w1 + n * thick; tri.vertices[2] = w2 + n * thick;
        tri.vertices[3] = w0 - n * thick; tri.vertices[4] = w1 - n * thick; tri.vertices[5] = w2 - n * thick;
        Contact c = detectContact(other, tri, box);   // normal tri(=mesh) -> other
        if (c.hit) out.push_back(c);
    }
}

// All contacts for a pair, normal always pointing b -> a. Box-box uses a clipping
// manifold; a mesh expands into per-triangle contacts; a compound expands into
// per-child contacts; everything else is one narrow-phase point.
void gatherContacts(const Body& a, const Body& b, const Box& box, std::vector<Contact>& out) {
    if (b.shape == Shape::Mesh) { meshContacts(a, b, box, out); return; }        // mesh(b) -> a == b -> a
    if (a.shape == Shape::Mesh) {
        meshContacts(b, a, box, out);
        for (Contact& c : out) c.normal = -c.normal;                            // -> b -> a
        return;
    }
    // Compound: test each child (as a temp body) against the other; normals are
    // reported b -> a as usual. Compound-vs-compound expands both.
    if (a.shape == Shape::Compound || b.shape == Shape::Compound) {
        if (a.shape == Shape::Compound && b.shape == Shape::Compound) {
            for (const ChildShape& ca : a.children) {
                Body ba = compoundChildWorld(a, ca);
                for (const ChildShape& cb : b.children) {
                    Body bb = compoundChildWorld(b, cb);
                    Contact c = detectContact(ba, bb, box);                     // bb -> ba == b -> a
                    if (c.hit) out.push_back(c);
                }
            }
        } else if (a.shape == Shape::Compound) {
            for (const ChildShape& ca : a.children) {
                Body ba = compoundChildWorld(a, ca);
                Contact c = detectContact(ba, b, box);                          // b -> ba == b -> a
                if (c.hit) out.push_back(c);
            }
        } else {
            for (const ChildShape& cb : b.children) {
                Body bb = compoundChildWorld(b, cb);
                Contact c = detectContact(a, bb, box);                          // bb(=b) -> a == b -> a
                if (c.hit) out.push_back(c);
            }
        }
        return;
    }
    Contact c = detectContact(a, b, box);
    if (!c.hit) return;
    if (a.shape == Shape::Box && b.shape == Shape::Box) {
        Manifold mf = boxBoxManifold(a, b, c.normal, box);
        for (int p = 0; p < mf.count; ++p) out.push_back({true, c.normal, mf.point[p], mf.depth[p]});
        if (mf.count > 0) return;
    }
    out.push_back(c);
}
}  // namespace

void World::step() {
    applyJointForces();     // spring joints add forces before integration
    for (Body& b : bodies) { b.prevX = b.x; b.prevQ = b.q; }   // snapshot: CCD + interpolation
    integrate();
    ccdPass();              // stop fast CCD bodies before they tunnel
    collide();
    solveRigidJoints(8);    // bilateral distance constraints
    solveJoints(8);         // ball / hinge / fixed / slider
    if (box.periodic)
        wrapPositions();
    else if (!openBoundary)
        applyWalls();
    if (sleepEnabled) updateSleep();
}

// Spring-damper distance joints: F = -k (dist - rest) - c (relative speed) along
// the anchor line, applied at each anchor (so it also torques the bodies).
void World::applyJointForces() {
    for (const DistanceJoint& j : distanceJoints) {
        if (j.stiffness <= 0.0) continue;                 // rigid ones are impulse-solved
        Body& a = bodies[j.a];
        Body& b = bodies[j.b];
        V3 wa = a.x + a.q.rotate(j.localA);
        V3 wb = b.x + b.q.rotate(j.localB);
        V3 d = box.minImage(wa - wb);
        double dist = d.norm();
        if (dist < 1e-9) continue;
        V3 nrm = d / dist;
        V3 rA = wa - a.x, rB = box.minImage(wb - b.x);
        V3 vrel = (a.v + a.w.cross(rA)) - (b.v + b.w.cross(rB));
        double vn = vrel.dot(nrm);
        double f = -j.stiffness * (dist - j.rest) - j.damping * vn;
        V3 F = nrm * f;
        a.force += F;  a.torque += rA.cross(F);
        b.force -= F;  b.torque -= rB.cross(F);
    }
}

// Rigid distance joints: a bilateral velocity constraint (impulse may be either
// sign) with a Baumgarte bias that removes the length error over time.
void World::solveRigidJoints(int iterations) {
    for (int it = 0; it < iterations; ++it) {
        for (const DistanceJoint& j : distanceJoints) {
            if (j.stiffness > 0.0) continue;
            Body& a = bodies[j.a];
            Body& b = bodies[j.b];
            V3 wa = a.x + a.q.rotate(j.localA);
            V3 wb = b.x + b.q.rotate(j.localB);
            V3 d = box.minImage(wa - wb);
            double dist = d.norm();
            if (dist < 1e-9) continue;
            V3 nrm = d / dist;
            V3 rA = wa - a.x, rB = box.minImage(wb - b.x);
            double K = effMass(a, b, rA, rB, nrm);
            if (K <= 1e-18) continue;
            V3 vrel = (a.v + a.w.cross(rA)) - (b.v + b.w.cross(rB));
            double vn = vrel.dot(nrm);
            double bias = 0.2 * (dist - j.rest) / dt;     // Baumgarte length correction
            double dj = -(vn + bias) / K;                 // bilateral: no clamp
            V3 J = nrm * dj;
            a.v += J * a.invMass; a.w += a.applyInvInertiaWorld(rA.cross(J));
            b.v -= J * b.invMass; b.w -= b.applyInvInertiaWorld(rB.cross(J));
        }
    }
}

namespace {
V3 perpendicular(const V3& v) {
    V3 a = std::fabs(v.x) < 0.9 ? V3{1, 0, 0} : V3{0, 1, 0};
    return v.cross(a).normalized();
}
double clampd(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }
}  // namespace

// General joints: sequential impulse. Every type first solves the point (ball)
// constraint pinning the anchors; Hinge/Fixed then add angular constraints and,
// for Hinge, an optional motor and angle limit; Slider frees one axis.
void World::solveJoints(int iterations) {
    const double beta = 0.2;
    for (Joint& j : joints) { j.appliedImpulse = 0.0; j.motorAcc = 0.0; }

    for (int it = 0; it < iterations; ++it) {
        for (Joint& j : joints) {
            if (j.broken) continue;
            Body& A = bodies[j.a];
            Body& B = bodies[j.b];
            const Mat3 IA = A.invInertiaWorld(), IB = B.invInertiaWorld();

            // Slider frees the axis translation; the others pin the point fully.
            V3 rA = A.q.rotate(j.localA), rB = B.q.rotate(j.localB);
            V3 posErr = box.minImage((A.x + rA) - (B.x + rB));
            V3 vRel = (A.v + A.w.cross(rA)) - (B.v + B.w.cross(rB));

            if (j.type == JointType::Slider) {
                // Constrain the two translation axes perpendicular to the slider
                // axis (allow sliding along it); axis defined on A.
                V3 ax = A.q.rotate(j.axisA).normalized();
                V3 t1 = perpendicular(ax), t2 = ax.cross(t1).normalized();
                Mat3 K = Mat3::scaled(A.invMass + B.invMass)
                         - Mat3::skew(rA) * IA * Mat3::skew(rA)
                         - Mat3::skew(rB) * IB * Mat3::skew(rB);
                for (const V3& t : {t1, t2}) {
                    double jv = vRel.dot(t) + posErr.dot(t) * (beta / dt);
                    double k = t.dot(K * t);
                    if (k > 1e-18) {
                        double l = -jv / k;
                        V3 P = t * l;
                        A.v += P * A.invMass; A.w += IA * (rA.cross(P));
                        B.v -= P * B.invMass; B.w -= IB * (rB.cross(P));
                        j.appliedImpulse += std::fabs(l);
                        vRel = (A.v + A.w.cross(rA)) - (B.v + B.w.cross(rB));
                    }
                }
            } else {
                Mat3 K = Mat3::scaled(A.invMass + B.invMass)
                         - Mat3::skew(rA) * IA * Mat3::skew(rA)
                         - Mat3::skew(rB) * IB * Mat3::skew(rB);
                V3 P = K.solve(-(vRel + posErr * (beta / dt)));
                A.v += P * A.invMass; A.w += IA * (rA.cross(P));
                B.v -= P * B.invMass; B.w -= IB * (rB.cross(P));
                j.appliedImpulse += P.norm();
            }

            // Angular constraints.
            if (j.type == JointType::Fixed || j.type == JointType::Slider) {
                // Lock the full relative orientation to the reference.
                V3 wRel = A.w - B.w;
                Q dq = (A.q * B.q.conjugate()) * j.refRel.conjugate();
                V3 angErr{dq.x, dq.y, dq.z};
                if (dq.w < 0.0) angErr = -angErr;
                angErr = angErr * 2.0;
                Mat3 Kang = IA + IB;
                V3 L = Kang.solve(-(wRel + angErr * (beta / dt)));
                A.w += IA * L; B.w -= IB * L;
            } else if (j.type == JointType::Hinge) {
                V3 aA = A.q.rotate(j.axisA).normalized();
                V3 aB = B.q.rotate(j.axisB).normalized();
                V3 t1 = perpendicular(aA), t2 = aA.cross(t1).normalized();
                V3 alignErr = aA.cross(aB);              // zero when the axes align
                Mat3 Kang = IA + IB;
                for (const V3& t : {t1, t2}) {
                    V3 wRel = A.w - B.w;
                    double jv = wRel.dot(t) - alignErr.dot(t) * (beta / dt);   // Baumgarte: reduce misalignment
                    double k = t.dot(Kang * t);
                    if (k > 1e-18) {
                        V3 L = t * (-jv / k);
                        A.w += IA * L; B.w -= IB * L;
                    }
                }
                // Motor about the hinge axis (clamped total impulse per step).
                if (j.useMotor) {
                    double k = aA.dot(Kang * aA);
                    if (k > 1e-18) {
                        double jv = (A.w - B.w).dot(aA) - j.motorSpeed;
                        double l = -jv / k;
                        double old = j.motorAcc;
                        double nn = clampd(old + l, -j.maxMotorImpulse, j.maxMotorImpulse);
                        l = nn - old; j.motorAcc = nn;
                        V3 L = aA * l; A.w += IA * L; B.w -= IB * L;
                    }
                }
                // Angle limit about the hinge axis.
                if (j.useLimit) {
                    Q dq = (A.q * B.q.conjugate()) * j.refRel.conjugate();
                    V3 dv{dq.x, dq.y, dq.z};
                    double angle = 2.0 * std::atan2(dv.dot(aA), dq.w);
                    double c = 0.0; bool active = false;
                    if (angle <= j.lower) { c = angle - j.lower; active = true; }
                    else if (angle >= j.upper) { c = angle - j.upper; active = true; }
                    if (active) {
                        double k = aA.dot(Kang * aA);
                        if (k > 1e-18) {
                            double jv = (A.w - B.w).dot(aA) + c * (beta / dt);
                            V3 L = aA * (-jv / k); A.w += IA * L; B.w -= IB * L;
                        }
                    }
                }
            }
        }
    }
    for (Joint& j : joints)
        if (j.breakable && !j.broken && j.appliedImpulse > j.breakImpulse) j.broken = true;
}

void World::wake(Body& b) {
    if (!b.sleeping) return;
    b.invMass = b.invMassStore;
    b.invInertiaBody = b.invInertiaStore;
    b.sleeping = false;
    b.sleepTimer = 0.0;
}

// Island-based sleeping: a connected group of awake dynamic bodies (linked by
// contacts or joints) sleeps only when EVERY member has been slow long enough, so
// a settling stack doesn't sleep one box at a time and re-wake its neighbours.
void World::updateSleep() {
    const std::size_t n = bodies.size();
    auto awakeDyn = [&](std::size_t i) { return bodies[i].dynamic && !bodies[i].sleeping && bodies[i].invMass > 0.0; };

    std::vector<std::size_t> parent(n);
    for (std::size_t i = 0; i < n; ++i) parent[i] = i;
    auto find = [&](std::size_t x) {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };
    auto join = [&](std::size_t a, std::size_t b) {
        if (awakeDyn(a) && awakeDyn(b)) parent[find(a)] = find(b);
    };
    for (const auto& pr : islandPairs_) join(pr.first, pr.second);
    for (const DistanceJoint& j : distanceJoints) join(j.a, j.b);
    for (const Joint& j : joints) if (!j.broken) join(j.a, j.b);

    // Advance each awake body's slow-timer.
    for (std::size_t i = 0; i < n; ++i) {
        if (!awakeDyn(i)) continue;
        Body& b = bodies[i];
        if (b.v.norm() < sleepLinVel && b.w.norm() < sleepAngVel) b.sleepTimer += dt;
        else b.sleepTimer = 0.0;
    }
    // An island is ready iff its slowest member's timer has reached sleepTime.
    std::unordered_map<std::size_t, double> islandMin;
    for (std::size_t i = 0; i < n; ++i) {
        if (!awakeDyn(i)) continue;
        std::size_t r = find(i);
        auto it = islandMin.find(r);
        if (it == islandMin.end() || bodies[i].sleepTimer < it->second) islandMin[r] = bodies[i].sleepTimer;
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!awakeDyn(i)) continue;
        if (islandMin[find(i)] < sleepTime) continue;
        Body& b = bodies[i];
        b.invMassStore = b.invMass;
        b.invInertiaStore = b.invInertiaBody;
        b.invMass = 0.0;
        b.invInertiaBody = {0, 0, 0};
        b.v = {}; b.w = {};
        b.sleeping = true;
    }
}

void World::integrate() {
    for (Body& b : bodies) {
        if (b.kinematic) {               // infinite mass, but its pose follows its
            b.x += b.v * dt;             // own velocity (user-driven); no forces,
            b.q = integrateOrientation(b.q, b.w, dt);   // no gravity, no impulses
            b.force = {};
            b.torque = {};
            continue;
        }
        if (b.invMass <= 0.0) {          // static / sleeping: never moves
            b.force = {};
            b.torque = {};
            continue;
        }
        // Linear: symplectic Euler. Gravity is an acceleration (mass-independent),
        // scaled per body; applied forces divide by mass.
        b.v += (gravity * b.gravityScale + b.force * b.invMass) * dt;
        b.x += b.v * dt;

        // Rotational: integrate ANGULAR MOMENTUM, not angular velocity. For an
        // asymmetric body (the cylinder, Ixx != Iyy) a torque-free w is NOT
        // constant -- only L is -- so advancing w directly would leak angular
        // momentum. Advance L by the torque, rotate by the current w, then
        // recompute w = I_world(q_new)^-1 L so L is conserved exactly at zero
        // torque. (For the isotropic sphere this reduces to constant w.)
        V3 L = b.angularMomentum();      // I_world(q) w, at the current pose
        L += b.torque * dt;
        b.q = integrateOrientation(b.q, b.w, dt);
        b.w = b.applyInvInertiaWorld(L); // w consistent with L at the new pose

        // Per-body velocity damping (drag), applied implicitly so it is stable
        // for any dt. 0 leaves v/w untouched, so undamped bodies are unchanged.
        if (b.linearDamping > 0.0) b.v = b.v * (1.0 / (1.0 + b.linearDamping * dt));
        if (b.angularDamping > 0.0) b.w = b.w * (1.0 / (1.0 + b.angularDamping * dt));

        b.force = {};
        b.torque = {};
    }
}

// All-pairs hard contact. O(N^2) with minimum-image separation; a cell list
// replaces this once correctness is established. Iterating i<j in index order
// keeps the impulse application deterministic. Cylinder-cylinder narrow phase
// (GJK/EPA) is added in the next increment; here we handle sphere-sphere and
// sphere-cylinder.
// PBC-aware uniform-grid broadphase. Cell size = 2*max bounding radius, so any
// overlapping pair lies in the same or an adjacent cell. Candidate pairs are
// normalised to (i<j), de-duplicated, and SORTED, so the downstream constraint
// list is identical (and bit-for-bit deterministic) to the brute-force scan --
// the grid only prunes pairs that could not touch. Falls back to O(N^2) for
// small N or a degenerate grid.
std::vector<std::pair<std::size_t, std::size_t>> World::broadphasePairs() const {
    const size_t n = bodies.size();
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    double maxR = 0.0;
    bool anyPlane = false;
    for (const Body& b : bodies) {
        if (b.gridExcluded()) { anyPlane = true; continue; }        // infinite/large: not gridded
        maxR = std::max(maxR, b.boundingRadius());
    }
    const double cell = 2.0 * maxR;

    if (forceBruteForce || n < 64 || cell < 1e-9) {
        for (size_t i = 0; i < n; ++i)
            for (size_t j = i + 1; j < n; ++j) pairs.push_back({i, j});
        return pairs;
    }

    // Planes are grid-excluded (infinite extent); pair each with every other body.
    auto appendPlanePairs = [&]() {
        if (!anyPlane) return;
        for (size_t p = 0; p < n; ++p) {
            if (!bodies[p].gridExcluded()) continue;
            for (size_t j = 0; j < n; ++j)
                if (j != p) pairs.push_back({std::min(p, j), std::max(p, j)});
        }
    };

    // Per-axis cell counts for the periodic wrap (>=1). Reflective uses unbounded
    // integer cells (no wrap).
    int nc[3] = {1, 1, 1};
    if (box.periodic) {
        double L = box.edge();
        for (int a = 0; a < 3; ++a) nc[a] = std::max(1, (int)std::floor(L / cell));
    }
    auto cellOf = [&](const V3& p, int axis) -> long {
        double c = axis == 0 ? p.x : (axis == 1 ? p.y : p.z);
        if (box.periodic) {
            long k = (long)std::floor((c + box.half) / cell);
            long m = nc[axis];
            return ((k % m) + m) % m;
        }
        return (long)std::floor(c / cell);
    };
    auto key = [](long x, long y, long z) {
        return (std::uint64_t)(x & 0x1FFFFF) | ((std::uint64_t)(y & 0x1FFFFF) << 21) |
               ((std::uint64_t)(z & 0x1FFFFF) << 42);
    };

    std::unordered_map<std::uint64_t, std::vector<size_t>> grid;
    grid.reserve(n * 2);
    for (size_t i = 0; i < n; ++i) {
        if (bodies[i].gridExcluded()) continue;
        long cx = cellOf(bodies[i].x, 0), cy = cellOf(bodies[i].x, 1), cz = cellOf(bodies[i].x, 2);
        grid[key(cx, cy, cz)].push_back(i);
    }
    for (size_t i = 0; i < n; ++i) {
        if (bodies[i].gridExcluded()) continue;
        long cx = cellOf(bodies[i].x, 0), cy = cellOf(bodies[i].x, 1), cz = cellOf(bodies[i].x, 2);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    long nx = cx + dx, ny = cy + dy, nz = cz + dz;
                    if (box.periodic) {
                        nx = ((nx % nc[0]) + nc[0]) % nc[0];
                        ny = ((ny % nc[1]) + nc[1]) % nc[1];
                        nz = ((nz % nc[2]) + nc[2]) % nc[2];
                    }
                    auto it = grid.find(key(nx, ny, nz));
                    if (it == grid.end()) continue;
                    for (size_t j : it->second)
                        if (j != i) pairs.push_back({std::min(i, j), std::max(i, j)});
                }
    }
    appendPlanePairs();
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

void World::collide() {
    // 1. Detect: build a constraint per overlapping pair, in sorted (i<j) order.
    std::vector<Constraint> cons;
    if (captureContacts) debugContacts.clear();
    islandPairs_.clear();
    for (const auto& pr : broadphasePairs()) {
        {
            Body& a = bodies[pr.first];
            Body& b = bodies[pr.second];
            const size_t i = pr.first, j = pr.second;
            if (a.invMass + b.invMass <= 0.0) continue;
            if (!layersCollide(a, b)) continue;         // collision filtering
            if (a.sensor || b.sensor) continue;         // triggers: no physical response

            std::vector<Contact> contacts;
            gatherContacts(a, b, box, contacts);        // each normal points b -> a
            if (contacts.empty()) continue;

            const double e = combineMat(restitutionCombine, effRestitution(a), effRestitution(b));
            const double mu = combineMat(frictionCombine, effFriction(a), effFriction(b));

            // Wake a sleeping body if a genuinely moving awake body contacts it;
            // if both are asleep the pair is stable and needs no constraint.
            if (sleepEnabled) {
                if (a.sleeping && !b.sleeping && (b.dynamic || b.kinematic) && b.v.norm() > wakeVel) wake(a);
                if (b.sleeping && !a.sleeping && (a.dynamic || a.kinematic) && a.v.norm() > wakeVel) wake(b);
                if (a.sleeping && b.sleeping) continue;
            }

            if (a.invMass > 0.0 && b.invMass > 0.0) islandPairs_.push_back({i, j});   // for island sleeping

            const std::uint64_t key = pairKeyOf(i, j);
            int slot = 0;
            for (const Contact& c : contacts) {
                Constraint k;
                k.i = i; k.j = j; k.n = c.normal; k.overlap = c.overlap; k.mu = mu;
                k.rA = box.minImage(c.point - a.x);
                k.rB = box.minImage(c.point - b.x);
                V3 vc = (a.v + a.w.cross(k.rA)) - (b.v + b.w.cross(k.rB));
                double vn = vc.dot(k.n);
                k.restBias = (vn < -kRestitutionSlop) ? -e * vn : 0.0;
                k.pairKey = key; k.slot = slot++;
                if (warmStart) {                          // seed from last step's impulse
                    auto it = warmCache_.find(key);
                    if (it != warmCache_.end() && k.slot < (int)it->second.size())
                        k.nImp = it->second[k.slot];
                }
                cons.push_back(k);
                if (captureContacts) debugContacts.push_back({c.point, c.normal, c.overlap, i, j});
            }
        }
    }

    // Warm start: apply each contact's seeded normal impulse up front, so the
    // iterations begin near last step's solution (stable, fast-settling stacks).
    if (warmStart) {
        for (Constraint& k : cons) {
            if (k.nImp == 0.0) continue;
            Body& a = bodies[k.i];
            Body& b = bodies[k.j];
            V3 J = k.n * k.nImp;
            a.v += J * a.invMass; a.w += a.applyInvInertiaWorld(k.rA.cross(J));
            b.v -= J * b.invMass; b.w -= b.applyInvInertiaWorld(k.rB.cross(J));
        }
    }

    // 2. Velocity solve: sequential impulse with accumulated clamping. Several
    // iterations so simultaneous contacts on one body settle consistently.
    const int kIters = 8;
    for (int it = 0; it < kIters; ++it) {
        for (Constraint& k : cons) {
            Body& a = bodies[k.i];
            Body& b = bodies[k.j];
            // Normal impulse toward the target separating speed (restBias).
            V3 vc = (a.v + a.w.cross(k.rA)) - (b.v + b.w.cross(k.rB));
            double vn = vc.dot(k.n);
            double Kn = effMass(a, b, k.rA, k.rB, k.n);
            if (Kn > 1e-18) {
                double dj = (k.restBias - vn) / Kn;
                double old = k.nImp;
                k.nImp = (old + dj > 0.0) ? old + dj : 0.0;   // clamp >= 0
                dj = k.nImp - old;
                V3 J = k.n * dj;
                a.v += J * a.invMass; a.w += a.applyInvInertiaWorld(k.rA.cross(J));
                b.v -= J * b.invMass; b.w -= b.applyInvInertiaWorld(k.rB.cross(J));
            }
            // Friction: oppose tangential velocity, clamped to the Coulomb cone.
            if (k.mu > 0.0 && k.nImp > 0.0) {
                V3 vc2 = (a.v + a.w.cross(k.rA)) - (b.v + b.w.cross(k.rB));
                V3 vt = vc2 - k.n * vc2.dot(k.n);
                double vtl = vt.norm();
                if (vtl > 1e-12) {
                    V3 t = vt / vtl;
                    double Kt = effMass(a, b, k.rA, k.rB, t);
                    if (Kt > 1e-18) {
                        double djt = -vtl / Kt;
                        double maxF = k.mu * k.nImp;
                        double old = k.tImp;
                        double sum = old + djt;
                        k.tImp = (sum < -maxF) ? -maxF : (sum > maxF ? maxF : sum);
                        djt = k.tImp - old;
                        V3 Jt = t * djt;
                        a.v += Jt * a.invMass; a.w += a.applyInvInertiaWorld(k.rA.cross(Jt));
                        b.v -= Jt * b.invMass; b.w -= b.applyInvInertiaWorld(k.rB.cross(Jt));
                    }
                }
            }
        }
    }

    // 3. Positional correction (geometric, energy-neutral).
    for (const Constraint& k : cons)
        correctPenetration(bodies[k.i], bodies[k.j], k.n, k.overlap, contactBeta);

    // Store this step's converged normal impulses for next step's warm start.
    if (warmStart) {
        warmCache_.clear();
        for (const Constraint& k : cons) {
            std::vector<double>& slots = warmCache_[k.pairKey];
            if ((int)slots.size() <= k.slot) slots.resize(k.slot + 1, 0.0);
            slots[k.slot] = k.nImp;
        }
    }
}

// Reflective walls at +-half, shape-aware for ALL shapes (a sphere-only version
// let cylinders and boxes fall straight through). The body's furthest point
// toward each wall is its support point; if that point crosses the wall the
// centre is projected back so the support just touches, and the outward velocity
// component is reflected (scaled by restitution). Center-based reflection: no
// wall torque, but no escape and correct extent per shape.
void World::applyWalls() {
    const double h = box.half;
    for (Body& b : bodies) {
        if (b.invMass <= 0.0) continue;
        double* px[3] = {&b.x.x, &b.x.y, &b.x.z};
        double* pv[3] = {&b.v.x, &b.v.y, &b.v.z};
        for (int axis = 0; axis < 3; ++axis) {
            for (int sgn = -1; sgn <= 1; sgn += 2) {
                V3 outward{axis == 0 ? double(sgn) : 0.0,
                           axis == 1 ? double(sgn) : 0.0,
                           axis == 2 ? double(sgn) : 0.0};
                V3 sp = supportWorld(b, outward);
                double coord = (axis == 0 ? sp.x : axis == 1 ? sp.y : sp.z) * sgn;
                double pen = coord - h;
                if (pen > 0.0) {
                    *px[axis] -= double(sgn) * pen;                 // support -> on the wall
                    if (*pv[axis] * sgn > 0.0)                      // moving into the wall
                        *pv[axis] = -restitution * *pv[axis];
                }
            }
        }
    }
}

// Continuous collision: for each fast CCD body, sweep a sphere of its bounding
// radius from its pre-step pose along this step's displacement; if it would reach
// another body's surface before the move completes, clamp it to that impact and
// cancel the velocity going into the surface. The next discrete step then resolves
// the contact normally. Conservative (bounding-sphere sweep) -> never tunnels.
void World::ccdPass() {
    const size_t n = bodies.size();
    for (size_t i = 0; i < n; ++i) {
        Body& bi = bodies[i];
        if (!bi.ccd || bi.invMass <= 0.0 || bi.sleeping) continue;
        V3 disp = bi.x - bi.prevX;
        double dist = disp.norm();
        double R = bi.boundingRadius();
        if (dist < R * 0.5 || dist < 1e-9) continue;   // too slow to tunnel
        V3 dir = disp / dist;
        double bestT = dist; V3 bestN; bool hit = false;
        for (size_t j = 0; j < n; ++j) {
            if (j == i) continue;
            const Body& bj = bodies[j];
            if (bj.sensor) continue;
            if (!layersCollide(bi, bj)) continue;
            Body infl = bj;                              // inflate bj by the sweep radius
            infl.radius += R;
            if (bj.shape == Shape::Box) infl.halfExtents = bj.halfExtents + V3{R, R, R};
            else if (bj.shape == Shape::Cylinder) infl.halfHeight = bj.halfHeight + R;
            V3 cImg = bi.prevX + box.minImage(bj.x - bi.prevX);
            double t; V3 p, nn;
            if (rayVsBody(infl, cImg, bi.prevX, dir, bestT, t, p, nn) && t < bestT) {
                bestT = t; bestN = nn; hit = true;
            }
        }
        if (hit) {
            bi.x = bi.prevX + dir * bestT;               // stop at the impact
            double vn = bi.v.dot(bestN);
            if (vn < 0.0) bi.v -= bestN * vn;            // cancel motion into the surface
        }
    }
}

void World::wrapPositions() {
    for (Body& b : bodies)
        if (b.invMass > 0.0)
            b.x = box.wrap(b.x);
}

V3 World::totalLinearMomentum() const {
    V3 p;
    for (const Body& b : bodies)
        if (b.invMass > 0.0) p += b.v * (1.0 / b.invMass);
    return p;
}

V3 World::totalAngularMomentum() const {
    // L_origin = sum( r x m v + I_world w ). Reflective only -- under PBC the
    // origin-referenced L is not conserved (wrapping shifts r by a lattice vector).
    V3 L;
    for (const Body& b : bodies) {
        if (b.invMass <= 0.0) continue;
        double m = 1.0 / b.invMass;
        L += b.x.cross(b.v * m);
        L += b.angularMomentum();
    }
    return L;
}

double World::totalKinetic() const {
    double e = 0.0;
    for (const Body& b : bodies)
        if (b.invMass > 0.0) e += b.kinetic();
    return e;
}

}  // namespace ne
