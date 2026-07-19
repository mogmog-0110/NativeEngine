#include "physics_world.hpp"

#include <algorithm>
#include <cstring>

#include "detect.hpp"
#include "queries.hpp"

namespace ne {

namespace {
template <class T> void put(std::vector<char>& b, const T& v) {
    const char* p = reinterpret_cast<const char*>(&v);
    b.insert(b.end(), p, p + sizeof(T));
}
template <class T> bool get(const std::vector<char>& b, std::size_t& off, T& v) {
    if (off + sizeof(T) > b.size()) return false;
    std::memcpy(&v, b.data() + off, sizeof(T));
    off += sizeof(T);
    return true;
}
constexpr std::uint32_t kSnapMagic = 0x4E45'5331;   // "NES1"
}  // namespace

std::vector<char> PhysicsWorld::saveState() const {
    std::vector<char> b;
    put(b, kSnapMagic);
    put(b, (std::uint32_t)w_.bodies.size());
    for (const Body& bd : w_.bodies) {
        put(b, bd.x); put(b, bd.q); put(b, bd.v); put(b, bd.w);
        put(b, bd.invMass); put(b, bd.invInertiaBody);
        put(b, bd.invMassStore); put(b, bd.invInertiaStore);
        put(b, bd.sleepTimer);
        put(b, (std::uint8_t)(bd.sleeping ? 1 : 0));
    }
    put(b, (std::uint32_t)gjoints_.size());
    for (const Joint& j : gjoints_) put(b, (std::uint8_t)(j.broken ? 1 : 0));
    return b;
}

bool PhysicsWorld::loadState(const std::vector<char>& data) {
    std::size_t off = 0;
    std::uint32_t magic = 0, n = 0;
    if (!get(data, off, magic) || magic != kSnapMagic) return false;
    if (!get(data, off, n) || n != w_.bodies.size()) return false;
    for (Body& bd : w_.bodies) {
        std::uint8_t sleeping = 0;
        if (!get(data, off, bd.x) || !get(data, off, bd.q) || !get(data, off, bd.v) ||
            !get(data, off, bd.w) || !get(data, off, bd.invMass) || !get(data, off, bd.invInertiaBody) ||
            !get(data, off, bd.invMassStore) || !get(data, off, bd.invInertiaStore) ||
            !get(data, off, bd.sleepTimer) || !get(data, off, sleeping))
            return false;
        bd.sleeping = (sleeping != 0);
    }
    std::uint32_t jn = 0;
    if (!get(data, off, jn) || jn != gjoints_.size()) return false;
    for (Joint& j : gjoints_) {
        std::uint8_t broken = 0;
        if (!get(data, off, broken)) return false;
        j.broken = (broken != 0);
    }
    w_.resetWarmCache();   // impulse cache is transient; rebuilds next step
    return true;
}

namespace {
bool layerAllowed(std::uint32_t mask, const Body& b) { return (mask & (1u << (b.layer & 31u))) != 0; }

void drawCircle(const PhysicsWorld::LineFn& line, const V3& c, const V3& u, const V3& v,
                double r, const V3& col, int segs = 20) {
    V3 prev = c + u * r;
    for (int i = 1; i <= segs; ++i) {
        double a = (double)i / segs * 2.0 * PI;
        V3 p = c + (u * std::cos(a) + v * std::sin(a)) * r;
        line(prev, p, col); prev = p;
    }
}
void drawAabb(const PhysicsWorld::LineFn& line, const V3& c, double h, const V3& col) {
    V3 vtx[8];
    for (int i = 0; i < 8; ++i)
        vtx[i] = {c.x + ((i & 1) ? h : -h), c.y + ((i & 2) ? h : -h), c.z + ((i & 4) ? h : -h)};
    int e[12][2] = {{0,1},{1,3},{3,2},{2,0},{4,5},{5,7},{7,6},{6,4},{0,4},{1,5},{2,6},{3,7}};
    for (auto& ed : e) line(vtx[ed[0]], vtx[ed[1]], col);
}
}  // namespace

void PhysicsWorld::debugDraw(const LineFn& line, std::uint32_t flags) const {
    const V3 cCol{0.4, 0.8, 0.5}, aCol{0.3, 0.5, 0.35}, vCol{0.4, 0.9, 0.5};
    const V3 ctCol{1.0, 0.9, 0.2}, nCol{0.2, 0.85, 1.0}, jCol{0.9, 0.6, 0.2}, comCol{1, 1, 1};

    for (const Body& b : w_.bodies) {
        if (flags & DbgColliders) {
            Mat3 R = b.rotationMatrix();
            V3 ex{R.m[0], R.m[3], R.m[6]}, ey{R.m[1], R.m[4], R.m[7]}, ez{R.m[2], R.m[5], R.m[8]};
            if (b.shape == Shape::Sphere) {
                drawCircle(line, b.x, ex, ey, b.radius, cCol);
                drawCircle(line, b.x, ey, ez, b.radius, cCol);
                drawCircle(line, b.x, ex, ez, b.radius, cCol);
            } else if (b.shape == Shape::Box) {
                V3 he = b.halfExtents; V3 c[8];
                for (int i = 0; i < 8; ++i)
                    c[i] = b.x + ex * ((i & 1) ? he.x : -he.x) + ey * ((i & 2) ? he.y : -he.y) + ez * ((i & 4) ? he.z : -he.z);
                int e[12][2] = {{0,1},{1,3},{3,2},{2,0},{4,5},{5,7},{7,6},{6,4},{0,4},{1,5},{2,6},{3,7}};
                for (auto& ed : e) line(c[ed[0]], c[ed[1]], cCol);
            } else if (b.shape == Shape::Capsule || b.shape == Shape::Cylinder) {
                V3 top = b.x + ey * b.halfHeight, bot = b.x - ey * b.halfHeight;
                drawCircle(line, top, ex, ez, b.radius, cCol);
                drawCircle(line, bot, ex, ez, b.radius, cCol);
                for (int k = 0; k < 4; ++k) {
                    double a = k * PI / 2.0; V3 off = (ex * std::cos(a) + ez * std::sin(a)) * b.radius;
                    line(bot + off, top + off, cCol);
                }
            } else {
                drawAabb(line, b.x, b.boundingRadius(), cCol);   // convex / plane
            }
        }
        if (flags & DbgAABB) drawAabb(line, b.x, b.boundingRadius(), aCol);
        if ((flags & DbgVelocities) && b.invMass > 0.0) line(b.x, b.x + b.v * 0.1, vCol);
        if (flags & DbgCOM) {
            double s = 0.1;
            line(b.x - V3{s, 0, 0}, b.x + V3{s, 0, 0}, comCol);
            line(b.x - V3{0, s, 0}, b.x + V3{0, s, 0}, comCol);
            line(b.x - V3{0, 0, s}, b.x + V3{0, 0, s}, comCol);
        }
    }
    if (flags & DbgContacts) {
        for (const ContactViz& c : w_.debugContacts) {
            double s = 0.12;
            line(c.point - V3{s, 0, 0}, c.point + V3{s, 0, 0}, ctCol);
            line(c.point - V3{0, s, 0}, c.point + V3{0, s, 0}, ctCol);
            line(c.point - V3{0, 0, s}, c.point + V3{0, 0, s}, ctCol);
            line(c.point, c.point + c.normal * 0.4, nCol);
        }
    }
    if (flags & DbgJoints) {
        auto anchorWorld = [&](BodyId h, const V3& local) -> V3 {
            const Body* b = body(h); return b ? b->x + b->q.rotate(local) : V3{};
        };
        for (const FacadeJoint& j : joints_)
            line(anchorWorld(j.a, j.localA), anchorWorld(j.b, j.localB), jCol);
        for (const Joint& j : gjoints_)
            line(anchorWorld((BodyId)j.a, j.localA), anchorWorld((BodyId)j.b, j.localB), jCol);
    }
}

RayHit PhysicsWorld::raycast(const V3& origin, const V3& dir, double maxDist,
                             std::uint32_t mask) const {
    RayHit hit;
    V3 d = dir.normalized();
    double best = maxDist;
    for (const Body& b : w_.bodies) {
        if (!layerAllowed(mask, b)) continue;
        V3 cImg = origin + w_.box.minImage(b.x - origin);   // nearest periodic image
        double t; V3 p, n;
        if (rayVsBody(b, cImg, origin, d, best, t, p, n) && t < best) {
            best = t;
            hit.hit = true; hit.distance = t; hit.point = p; hit.normal = n;
            // resolve the body index back to its handle
            size_t idx = (size_t)(&b - w_.bodies.data());
            hit.body = handle_[idx]; hit.userData = b.userData;
        }
    }
    return hit;
}

std::vector<BodyId> PhysicsWorld::overlapSphere(const V3& center, double radius,
                                                std::uint32_t mask) const {
    std::vector<BodyId> out;
    Body q = makeSphere(-1, center, radius, 1.0);
    for (size_t i = 0; i < w_.bodies.size(); ++i) {
        const Body& b = w_.bodies[i];
        if (!layerAllowed(mask, b)) continue;
        if (detectContact(q, b, w_.box).hit) out.push_back(handle_[i]);
    }
    return out;
}

std::vector<BodyId> PhysicsWorld::overlapBox(const V3& center, const Q& rot, const V3& half,
                                             std::uint32_t mask) const {
    std::vector<BodyId> out;
    Body q = makeBox(-1, center, rot, half, 1.0);
    for (size_t i = 0; i < w_.bodies.size(); ++i) {
        const Body& b = w_.bodies[i];
        if (!layerAllowed(mask, b)) continue;
        if (detectContact(q, b, w_.box).hit) out.push_back(handle_[i]);
    }
    return out;
}

std::vector<PhysicsWorld::OverlapHit> PhysicsWorld::overlapContacts(const Body& query,
                                                                    std::uint32_t mask) const {
    std::vector<OverlapHit> out;
    for (size_t i = 0; i < w_.bodies.size(); ++i) {
        const Body& b = w_.bodies[i];
        if (b.sensor) continue;                       // triggers don't block
        if (!layerAllowed(mask, b)) continue;
        Contact c = detectContact(query, b, w_.box);   // normal b -> query (push-out)
        if (!c.hit) continue;
        out.push_back({handle_[i], b.userData, c.normal, c.overlap});
    }
    return out;
}

RayHit PhysicsWorld::sphereCast(const V3& origin, double radius, const V3& dir,
                                double maxDist, std::uint32_t mask) const {
    // Sweeping a sphere of radius R == raycasting against each body inflated by R
    // (its Minkowski sum with the sphere). Exact for spheres/capsules; a
    // conservative over-approximation at box/cylinder edges (safe: hits early).
    RayHit hit;
    V3 d = dir.normalized();
    double best = maxDist;
    for (const Body& b : w_.bodies) {
        if (!layerAllowed(mask, b)) continue;
        Body infl = b;
        infl.radius += radius;
        if (b.shape == Shape::Box) infl.halfExtents = b.halfExtents + V3{radius, radius, radius};
        else if (b.shape == Shape::Cylinder) infl.halfHeight = b.halfHeight + radius;
        V3 cImg = origin + w_.box.minImage(b.x - origin);
        double t; V3 p, n;
        if (rayVsBody(infl, cImg, origin, d, best, t, p, n) && t < best) {
            best = t;
            hit.hit = true; hit.distance = t; hit.normal = n;
            hit.point = p - n * radius;             // contact on the real surface
            size_t idx = (size_t)(&b - w_.bodies.data());
            hit.body = handle_[idx]; hit.userData = b.userData;
        }
    }
    return hit;
}

void PhysicsWorld::processContacts() {
    const std::vector<Body>& B = w_.bodies;
    const size_t n = B.size();

    // 1. Gather this step's overlapping pairs, split solid vs trigger.
    std::map<PairKey, ContactInfo> curContacts, curTriggers;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const bool anySensor = B[i].sensor || B[j].sensor;
            // Skip only if there is nothing to report: two immovable solids.
            if (!anySensor && B[i].invMass + B[j].invMass <= 0.0) continue;
            if (!World::layersCollide(B[i], B[j])) continue;
            if (World::pairExcluded(B[i], (std::uint32_t)j)) continue;   // bonded/jointed pair
            Contact c = detectContact(B[i], B[j], w_.box);
            if (!c.hit) continue;

            ContactInfo ci;
            ci.a = handle_[i]; ci.b = handle_[j];
            ci.aUser = B[i].userData; ci.bUser = B[j].userData;
            ci.point = c.point; ci.normal = c.normal; ci.depth = c.overlap;
            PairKey key{std::min(ci.a, ci.b), std::max(ci.a, ci.b)};
            (anySensor ? curTriggers : curContacts)[key] = ci;
        }
    }

    // 2. Solid pairs: legacy per-step callback + begin/stay transitions.
    for (const auto& kv : curContacts) {
        if (onContact_) onContact_(kv.second);
        const bool wasThere = prevContacts_.count(kv.first) != 0;
        if (!wasThere) { if (onBegin_) onBegin_(kv.second); }
        else { if (onStay_) onStay_(kv.second); }
    }
    // Solid pairs that ended (present last step, gone now).
    if (onEnd_)
        for (const auto& kv : prevContacts_)
            if (curContacts.count(kv.first) == 0) onEnd_(kv.second);

    // 3. Trigger pairs: enter/stay/exit.
    for (const auto& kv : curTriggers) {
        const bool wasThere = prevTriggers_.count(kv.first) != 0;
        if (!wasThere) { if (onTrigEnter_) onTrigEnter_(kv.second); }
        else { if (onTrigStay_) onTrigStay_(kv.second); }
    }
    if (onTrigExit_)
        for (const auto& kv : prevTriggers_)
            if (curTriggers.count(kv.first) == 0) onTrigExit_(kv.second);

    prevContacts_.swap(curContacts);
    prevTriggers_.swap(curTriggers);
}

}  // namespace ne
