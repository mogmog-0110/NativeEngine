#include "physics_world.hpp"

#include "detect.hpp"

namespace ne {

void PhysicsWorld::emitContacts() {
    const std::vector<Body>& B = w_.bodies;
    const size_t n = B.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (B[i].invMass + B[j].invMass <= 0.0) continue;
            if (!World::layersCollide(B[i], B[j])) continue;   // same filter as the solver
            Contact c = detectContact(B[i], B[j], w_.box);
            if (!c.hit) continue;
            ContactInfo ci;
            ci.a = handle_[i];
            ci.b = handle_[j];
            ci.aUser = B[i].userData;
            ci.bUser = B[j].userData;
            ci.point = c.point;
            ci.normal = c.normal;
            ci.depth = c.overlap;
            onContact_(ci);
        }
    }
}

}  // namespace ne
