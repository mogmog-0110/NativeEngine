#pragma once
#ifndef NATIVEENGINE_CHARACTER_HPP
#define NATIVEENGINE_CHARACTER_HPP

#include "body.hpp"
#include "physics_world.hpp"

namespace ne {

// A kinematic capsule character controller: not a rigid body, but moved by
// "collide and slide" against the world. The capsule stands upright (axis +Y),
// with cylinder half-height `halfHeight` and radius `radius` (total height
// 2*halfHeight + 2*radius). move() applies the desired displacement then recovers
// from penetration, which both stops the character at walls/floors and lets it
// slide along them; a short downward probe sets the grounded flag.
class CharacterController {
public:
    CharacterController(PhysicsWorld* world, double radius, double halfHeight,
                        std::uint32_t mask = 0xFFFFFFFFu)
        : world_(world), radius_(radius), halfHeight_(halfHeight), mask_(mask) {}

    void setPosition(const V3& p) { pos_ = p; }
    V3 position() const { return pos_; }
    bool grounded() const { return grounded_; }
    // A floor is walkable if its normal.y exceeds this cosine (default ~60 deg).
    void setSlopeLimit(double cosAngle) { slopeCos_ = cosAngle; }
    void setSkinWidth(double s) { skin_ = s; }

    // Move by the desired displacement, then resolve penetration (a few passes,
    // deepest first for stability). Displacement is caller-driven (apply gravity /
    // input outside), so this stays a pure kinematic sweep.
    void move(const V3& disp) {
        pos_ += disp;
        Body q = capsuleAt(pos_);
        for (int iter = 0; iter < 5; ++iter) {
            auto hits = world_->overlapContacts(q, mask_);
            if (hits.empty()) break;
            const PhysicsWorld::OverlapHit* deepest = &hits[0];
            for (const auto& h : hits) if (h.depth > deepest->depth) deepest = &h;
            if (deepest->depth <= skin_) break;                  // within the skin: settled
            pos_ += deepest->normal * (deepest->depth - skin_);  // push out, keep a skin
            q.x = pos_;
        }
        updateGrounded();
    }

private:
    Body capsuleAt(const V3& p) const {
        return makeCapsule(-1, p, Q{1, 0, 0, 0}, radius_, halfHeight_, 1.0);
    }
    void updateGrounded() {
        RayHit g = world_->sphereCast(pos_, radius_, V3{0, -1, 0}, halfHeight_ + 0.12, mask_);
        grounded_ = g.hit && g.normal.y > slopeCos_;
    }

    PhysicsWorld* world_;
    double radius_, halfHeight_;
    std::uint32_t mask_;
    V3 pos_;
    bool grounded_ = false;
    double slopeCos_ = 0.5;
    double skin_ = 0.01;
};

}  // namespace ne

#endif  // NATIVEENGINE_CHARACTER_HPP
