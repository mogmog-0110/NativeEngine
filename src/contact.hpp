#pragma once
#ifndef NATIVEENGINE_CONTACT_HPP
#define NATIVEENGINE_CONTACT_HPP

#include "body.hpp"

namespace ne {

// General rigid-body contact impulse at an off-centre point.
//
// This is the piece sphere-sphere didn't need (its contact normal passes through
// both centres, so the lever arms are parallel to n and produce no torque) but
// every flat-cap cylinder contact does: a hit on a cap or a side edge is
// off-centre and must impart spin. Standard formulation (Baraff / Catto):
//
//   v_rel = (v_a + w_a x r_a) - (v_b + w_b x r_b)          at the contact point
//   K     = 1/m_a + 1/m_b + n . ( Iinv_a (r_a x n) x r_a
//                               + Iinv_b (r_b x n) x r_b )  effective mass
//   j     = -(1+e) (v_rel . n) / K
//   apply +j n at r_a on a, -j n at r_b on b.
//
// Applying +J and -J at the SAME world contact point conserves linear momentum
// (equal/opposite) and angular momentum about any origin (the orbital terms
// cancel; the spin terms are the bodies' response). r_a, r_b are lever arms from
// each body centre to the contact point, already minimum-imaged by the caller,
// so this function is PBC-agnostic.
inline void resolveContact(Body& a, Body& b,
                           const V3& ra, const V3& rb,
                           const V3& n, double e) {
    V3 vContact = (a.v + a.w.cross(ra)) - (b.v + b.w.cross(rb));
    double vn = vContact.dot(n);
    if (vn >= 0.0) return;   // separating (or resting): no impulse

    V3 raxn = ra.cross(n);
    V3 rbxn = rb.cross(n);
    V3 angA = a.applyInvInertiaWorld(raxn).cross(ra);
    V3 angB = b.applyInvInertiaWorld(rbxn).cross(rb);
    double K = a.invMass + b.invMass + (angA + angB).dot(n);
    if (K <= 1e-18) return;

    double j = -(1.0 + e) * vn / K;
    V3 J = n * j;
    a.v += J * a.invMass;
    a.w += a.applyInvInertiaWorld(ra.cross(J));
    b.v -= J * b.invMass;
    b.w -= b.applyInvInertiaWorld(rb.cross(J));
}

// Positional (overlap) correction along the contact normal, split by inverse
// mass. Moving centres does not change kinetic energy or momentum. `beta` is the
// fraction of the overlap removed per step.
inline void correctPenetration(Body& a, Body& b, const V3& n,
                               double overlap, double beta) {
    double invSum = a.invMass + b.invMass;
    if (invSum <= 0.0 || overlap <= 0.0) return;
    V3 corr = n * (beta * overlap / invSum);
    a.x += corr * a.invMass;
    b.x -= corr * b.invMass;
}

}  // namespace ne

#endif  // NATIVEENGINE_CONTACT_HPP
