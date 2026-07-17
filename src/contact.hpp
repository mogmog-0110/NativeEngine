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
// Effective mass along a unit direction u for a contact at lever arms ra, rb.
inline double effMass(const Body& a, const Body& b,
                      const V3& ra, const V3& rb, const V3& u) {
    V3 angA = a.applyInvInertiaWorld(ra.cross(u)).cross(ra);
    V3 angB = b.applyInvInertiaWorld(rb.cross(u)).cross(rb);
    return a.invMass + b.invMass + (angA + angB).dot(u);
}

// mu = Coulomb friction coefficient (0 = frictionless). Friction is applied
// after the normal impulse, opposing the tangential relative velocity and
// clamped to mu * (normal impulse) -- the standard Coulomb cone. It imparts spin
// (a glancing hit makes a body roll) and dissipates tangential energy.
inline void resolveContact(Body& a, Body& b,
                           const V3& ra, const V3& rb,
                           const V3& n, double e, double mu = 0.0) {
    V3 vContact = (a.v + a.w.cross(ra)) - (b.v + b.w.cross(rb));
    double vn = vContact.dot(n);
    if (vn >= 0.0) return;   // separating (or resting): no impulse

    double K = effMass(a, b, ra, rb, n);
    if (K <= 1e-18) return;

    double jn = -(1.0 + e) * vn / K;
    V3 J = n * jn;
    a.v += J * a.invMass;
    a.w += a.applyInvInertiaWorld(ra.cross(J));
    b.v -= J * b.invMass;
    b.w -= b.applyInvInertiaWorld(rb.cross(J));

    if (mu <= 0.0) return;
    // Tangential relative velocity after the normal impulse.
    V3 vc2 = (a.v + a.w.cross(ra)) - (b.v + b.w.cross(rb));
    V3 vt = vc2 - n * vc2.dot(n);
    double vtl = vt.norm();
    if (vtl < 1e-12) return;
    V3 t = vt / vtl;
    double Kt = effMass(a, b, ra, rb, t);
    if (Kt <= 1e-18) return;
    double jt = -vtl / Kt;
    double maxF = mu * jn;                       // jn > 0 here
    if (jt < -maxF) jt = -maxF; else if (jt > maxF) jt = maxF;
    V3 Jt = t * jt;
    a.v += Jt * a.invMass;
    a.w += a.applyInvInertiaWorld(ra.cross(Jt));
    b.v -= Jt * b.invMass;
    b.w -= b.applyInvInertiaWorld(rb.cross(Jt));
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
