#include "physics_world.hpp"

#include <algorithm>

#include "detect.hpp"
#include "queries.hpp"

namespace ne {

namespace {
bool layerAllowed(std::uint32_t mask, const Body& b) { return (mask & (1u << (b.layer & 31u))) != 0; }
}  // namespace

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
