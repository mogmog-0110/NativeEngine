#include "physics_world.hpp"

#include <algorithm>

#include "detect.hpp"

namespace ne {

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
