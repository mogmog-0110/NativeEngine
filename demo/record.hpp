#pragma once
#ifndef NATIVEENGINE_DEMO_RECORD_HPP
#define NATIVEENGINE_DEMO_RECORD_HPP

// Record a NativeEngine World to the project's engine-agnostic .pxrf format, so
// the existing OpenGL playback viewer (x64\ReleaseRender\PhysxRender.exe) can
// replay a NativeEngine run for visual inspection. The viewer draws spheres and
// (Y-axis) cylinders, so demos use those shapes; boxes are not rendered there.

#include <string>
#include <vector>

#include "world.hpp"
#include "recording.hpp"   // reused unchanged from PhysxRender/ (PhysX-free)

namespace ne {

class Recorder {
public:
    // interval = record every `interval` steps; roomSize is the box half-edge
    // (the viewer uses it for the camera / wall box).
    void init(const World& w, int interval, double roomSize) {
        rec_.version = 2;
        rec_.recordInterval = (uint32_t)interval;
        rec_.roomSize = (float)roomSize;
        rec_.sphereRadius = 1.0f;
        rec_.cylinderRadius = 1.0f;
        rec_.cylinderHeight = 4.0f;
        rec_.actorMeta.clear();
        for (size_t i = 0; i < w.bodies.size(); ++i) {
            const Body& b = w.bodies[i];
            ActorMeta m;
            m.id = (int32_t)i;
            m.isSphere = b.isSphere() ? 1 : 0;      // boxes fall here too; avoid in demos
            m.radius = (float)b.radius;
            m.halfHeight = (float)b.halfHeight;
            rec_.actorMeta.push_back(m);
        }
        interval_ = interval;
    }

    void maybeCapture(const World& w, int step) {
        if (interval_ <= 0 || step % interval_ != 0) return;
        capture(w, (uint32_t)step);
    }

    void capture(const World& w, uint32_t step) {
        std::vector<ActorSnapshot> snaps;
        snaps.reserve(w.bodies.size());
        for (const Body& b : w.bodies) {
            ActorSnapshot s;
            s.px = (float)b.x.x; s.py = (float)b.x.y; s.pz = (float)b.x.z;
            s.qx = (float)b.q.x; s.qy = (float)b.q.y; s.qz = (float)b.q.z; s.qw = (float)b.q.w;
            snaps.push_back(s);
        }
        rec_.addFrame(step, snaps);
    }

    bool save(const std::string& path) const { return rec_.save(path); }
    uint32_t frames() const { return rec_.numFrames(); }

private:
    Recording rec_;
    int interval_ = 1;
};

}  // namespace ne

#endif  // NATIVEENGINE_DEMO_RECORD_HPP
