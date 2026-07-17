#include "physics_world.hpp"

#include "narrowphase.hpp"
#include "epa.hpp"

namespace ne {

// Contact detection for EVENT REPORTING only (read-only; the solver in World has
// already resolved the step). Normal points from b toward a, matching ContactInfo.
static Contact detectForReport(const Body& a, const Body& b, const Box& box) {
    Contact c;
    if (a.isSphere() && b.isSphere()) {
        V3 d = box.minImage(a.x - b.x);
        double sumR = a.radius + b.radius, dist2 = d.norm2();
        if (dist2 >= sumR * sumR || dist2 < 1e-18) return c;
        double dist = std::sqrt(dist2);
        c.hit = true;
        c.normal = d / dist;                       // b -> a
        c.overlap = sumR - dist;
        c.point = a.x - c.normal * a.radius;        // on a's surface
        return c;
    }
    const bool aCyl = a.shape == Shape::Cylinder;
    const bool bCyl = b.shape == Shape::Cylinder;
    if (a.isSphere() && bCyl) return sphereVsCylinder(a, b, box);   // normal b->a
    if (b.isSphere() && aCyl) {                                     // normal was a->b
        Contact t = sphereVsCylinder(b, a, box);
        if (t.hit) t.normal = -t.normal;                            // -> b->a
        return t;
    }
    return convexContact(a, b, box);               // normal b->a
}

void PhysicsWorld::emitContacts() {
    const std::vector<Body>& B = w_.bodies;
    const size_t n = B.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (B[i].invMass + B[j].invMass <= 0.0) continue;
            Contact c = detectForReport(B[i], B[j], w_.box);
            if (!c.hit) continue;
            ContactInfo ci;
            ci.a = handle_[i];
            ci.b = handle_[j];
            ci.point = c.point;
            ci.normal = c.normal;
            ci.depth = c.overlap;
            onContact_(ci);
        }
    }
}

}  // namespace ne
