/******************************************************************************
 * Cosserat Shell — Element Stiffness and Residual Assembler                 *
 *                                                                            *
 * Computes the element tangent stiffness K_e (24×24) and internal force     *
 * vector F_int_e (24×1) for one Q4 Cosserat shell element.                  *
 *                                                                            *
 * Strategy:                                                                  *
 *   - Strains evaluated at centroid (shear-locking free, constant per elem) *
 *   - K = K_M (material) + K_G (geometric)                                  *
 *   - Integration: 2×2 Gauss, but strains from centroid                     *
 *                                                                            *
 * Reference: cosserat_shell_v0.pdf, Sections 4-5                           *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                   *
 ******************************************************************************/
#pragma once

#include "ShellElement.h"
#include "SE3ShellKinematics.h"
#include "ConstitutiveLaw.h"
#include <Eigen/Dense>

namespace sofa::component::cosserat::shell {

/**
 * @brief Compute element internal force vector and tangent stiffness.
 *
 * Given the current configuration of the 4 nodes and the pre-computed
 * reference data, this function assembles:
 *
 *   F_int_e ∈ ℝ^{24}   — internal force vector (residual contribution)
 *   K_e     ∈ ℝ^{24×24} — element tangent stiffness (material + geometric)
 *
 * Both are returned by output parameters.
 *
 * Algorithm:
 * 1. Evaluate strain twists X_t at centroid → constant over element
 * 2. Compute differential strain E = X_t - X_0
 * 3. Compute stress resultants S^α via D^{αβ}
 * 4. Loop over 2×2 Gauss points for integration:
 *    a. Compute B_alpha at Gauss point (uses centroid strains, varies due to N^i)
 *    b. F_int += Σ_α B_α^T S^α j_0^g w_g
 *    c. K_M   += Σ_{α,β} B_α^T D^{αβ} B_β j_0^g w_g
 *    d. K_G   += Σ_α B_α^T Ãd_{S^α}^T B_α j_0^g w_g
 *
 * @param nodes      Current node configurations (position + orientation)
 * @param elemRef    Pre-computed reference element data
 * @param tensors    Constitutive material tensors D^{αβ}
 * @param[out] elem  Output: current element state (Fint, K, strains, stresses)
 */
inline void computeElementForceAndStiffness(
    const std::array<NodeConfig, NODES_PER_ELEM>& nodes,
    const ShellElementRef& elemRef,
    const ConstitutiveTensors& tensors,
    ShellElementCurrent& elem)
{
    // ── 1. Strain twists at centroid ──────────────────────────────────────
    elem.Xt = computeStrainTwistsAtCentroid(nodes);

    // ── 2. Differential strain E = X_t - X_0 ─────────────────────────────
    elem.E = elem.Xt - elemRef.X0;

    // ── 3. Stress resultants S^α ──────────────────────────────────────────
    elem.S = computeAllStressResultants(tensors, elem.E);

    // ── 4. Zero output ────────────────────────────────────────────────────
    elem.Fint.setZero();
    elem.K.setZero();

    // ── 5. Gauss integration ──────────────────────────────────────────────
    const auto gaussPts = gaussPoints2x2();

    for (int g = 0; g < N_GAUSS; ++g) {
        const double xi  = gaussPts[g].xi;
        const double eta = gaussPts[g].eta;
        const double w   = gaussPts[g].w;
        const double j0g = elemRef.j0Gauss[g];
        const double dA  = j0g * w;

        // B_alpha matrices at Gauss point g
        // (strains are constant from centroid, but B depends on N^i(xi,eta))
        Eigen::Matrix<double,6,DOF_PER_ELEM> B[2];
        for (int alpha = 0; alpha < 2; ++alpha) {
            B[alpha] = computeBMatrix(elem.Xt.col(alpha), xi, eta, alpha);
        }

        // ── F_int += Σ_α B_α^T · S^α · j0 · w ───────────────────────────
        for (int alpha = 0; alpha < 2; ++alpha) {
            elem.Fint += B[alpha].transpose() * elem.S.col(alpha) * dA;
        }

        // ── K_M += Σ_{α,β} B_α^T · D^{αβ} · B_β · j0 · w ───────────────
        for (int alpha = 0; alpha < 2; ++alpha) {
            for (int beta = 0; beta < 2; ++beta) {
                const Eigen::Matrix<double,6,6>& Dab = tensors.D[2*alpha + beta];
                elem.K += B[alpha].transpose() * Dab * B[beta] * dA;
            }
        }

        // ── K_G += Σ_α B_α^T · (-ad_{ζ_α}^T) · S^α_cross · B_α · j0 · w
        // From the proof: geometric stiffness comes from ⟨S^α, ad_{K_α η} κ⟩
        // = ⟨-ad*_{S^α} κ, K_α η⟩
        // The co-adjoint term: for a given S^α ∈ ℝ⁶, build a 6×6 matrix
        // such that the geometric contribution is:
        //   K_G = Σ_α B_α^T · ad_S^α_matrix^T · B_α
        // where ad_S^α_matrix η = ad_{η} S^α = -ad_{S^α} η (antisymmetry of ad)
        // i.e., the matrix [−ad_{S^α}] applied to η gives [η, S^α].
        //
        // More precisely: ⟨S^α, [K_α η, κ]⟩ = -⟨ad_{S^α} κ, K_α η⟩
        //                                    = κ^T (ad_{S^α})^T B_α η
        // → contribution to K: B_α^T (ad_{S^α})^T B_α (from the κ side)
        // But since we want K symmetric, we take the symmetric part:
        // K_G += ½ B_α^T [adS^T + adS] B_α — at equilibrium K_G is symmetric.
        for (int alpha = 0; alpha < 2; ++alpha) {
            const Eigen::Matrix<double,6,6> adS  = adMatrix(elem.S.col(alpha));
            // Geometric stiffness: -adS^T is the co-adjoint matrix
            const Eigen::Matrix<double,DOF_PER_ELEM,DOF_PER_ELEM> KG_elem =
                B[alpha].transpose() * (-adS.transpose()) * B[alpha] * dA;
            elem.K += KG_elem;
        }
    }
}

// ─── Reference element initializer ────────────────────────────────────────

/**
 * @brief Initialize reference element data from reference node configurations.
 *
 * Computes:
 *   - X_0 = strain twists at reference centroid
 *   - X_0* = pseudo-inverse of X_0
 *   - j_0 = area Jacobian at centroid
 *   - j_0^g = area Jacobian at each Gauss point
 *
 * @param refNodes   Reference node configurations
 * @param[out] elem  Initialized reference element data
 */
inline void initShellElementRef(
    const std::array<NodeConfig, NODES_PER_ELEM>& refNodes,
    ShellElementRef& elem)
{
    // ── Reference strains at centroid ─────────────────────────────────────
    elem.X0 = computeStrainTwistsAtCentroid(refNodes);

    // ── Pseudo-inverse X_0* = X_0^T (X_0 X_0^T)^{-1} ────────────────────
    const auto& X0 = elem.X0;
    const Eigen::Matrix<double,2,2> XtX = X0.transpose() * X0;
    // Use LDLT for the 2×2 symmetric system
    elem.X0star = XtX.ldlt().solve(X0.transpose());  // ∈ ℝ^{2×6}

    // ── Area Jacobian at centroid ─────────────────────────────────────────
    std::array<Eigen::Vector3d, NODES_PER_ELEM> refPositions;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        refPositions[i] = refNodes[i].position;
    }
    elem.j0 = computeAreaJacobian(refPositions, 0.0, 0.0);

    // ── Area Jacobian at Gauss points ─────────────────────────────────────
    const auto gaussPts = gaussPoints2x2();
    for (int g = 0; g < N_GAUSS; ++g) {
        elem.j0Gauss[g] = computeAreaJacobian(refPositions,
                                               gaussPts[g].xi,
                                               gaussPts[g].eta);
    }
}

// ─── SE(3) Newton update ──────────────────────────────────────────────────

/**
 * @brief Apply a Newton increment to a node configuration using SE(3) update.
 *
 *   φ_{new} = φ_{old} + δv   (translational increment)
 *   R_{new}  = R_{old} * exp(δω)  (SO(3) Lie update)
 *
 * The 6D increment is laid out as [δω (3) ; δv (3)].
 *
 * @param node      Current node configuration (modified in place)
 * @param increment 6D increment vector [δω; δv] ∈ ℝ⁶
 */
inline void applyNodeUpdate(NodeConfig& node,
                            const Eigen::Matrix<double,6,1>& increment)
{
    const Eigen::Vector3d delta_omega = increment.head<3>();
    const Eigen::Vector3d delta_v     = increment.tail<3>();

    // Rotation update: R_new = R_old * exp(δω)
    node.orientation = node.orientation * SO3d::exp(delta_omega);
    // Position update
    node.position += delta_v;
}

} // namespace sofa::component::cosserat::shell
