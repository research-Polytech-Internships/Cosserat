/******************************************************************************
 * Cosserat Shell -- Full element stiffness and residual assembler           *
 *                                                                            *
 * Extends ShellElementComputer.h with:                                       *
 *   1. Correct geometric stiffness K_G = B^T * ad_{S^alpha} * B             *
 *      (was -ad_{S^alpha}^T in the simple version)                           *
 *   2. Optional symmetrization of K_G                                        *
 *   3. Accepts ConstitutiveTensorsFull (orthotropic + cylindrical)           *
 *   4. Richer diagnostic output: K_M, K_G, W_int stored separately          *
 *                                                                            *
 * Comparison:                                                                *
 *   ShellElementComputer.h     -- simple, K_G approximate, isotropic only   *
 *   ShellElementComputerFull.h -- correct K_G, full material support         *
 *                                                                            *
 * Why K_G matters:                                                           *
 *   The geometric stiffness encodes the change in the B-matrix with the     *
 *   configuration. Without it (K_G=0) equilibrium is still correct but     *
 *   Newton may require many more iterations. With the wrong formula          *
 *   (-ad^T instead of ad), Newton can diverge for large deformations.       *
 *                                                                            *
 * Derivation (see cosserat_shell_formulation.md S11):                       *
 *   d/de <S^alpha, K_alpha delta_eta>|_{e=0}                                *
 *   = <S^alpha, ad_{K_alpha Delta_eta} delta_eta>                           *
 *   = delta_eta^T * ad_{S^alpha}^T * B_alpha * Delta_eta_e                 *
 *   -> K_G = sum_alpha B_alpha^T * ad_{S^alpha} * B_alpha                   *
 *                                                                            *
 * Reference: cosserat_shell_formulation.md S11, formulation_review.md S5    *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                   *
 ******************************************************************************/
#pragma once

#include "ShellElement.h"
#include "SE3ShellKinematics.h"
#include "ConstitutiveLawFull.h"
#include <Eigen/Dense>

namespace sofa::component::cosserat::shell {

// --- Extended element state -------------------------------------------------

/**
 * @brief Element state with separated K_M / K_G and strain energy.
 *
 * Provides richer diagnostic output than ShellElementCurrent, enabling:
 *   - Comparison of material vs. geometric stiffness norms
 *   - Verification that K_G -> 0 at the reference state
 *   - Energy tracking during Newton iterations
 */
struct ShellElementCurrentFull {
    // Kinematics
    Eigen::Matrix<Scalar, 6, 2> Xt  = Eigen::Matrix<Scalar,6,2>::Zero();  ///< Current strain twists
    Eigen::Matrix<Scalar, 6, 2> E   = Eigen::Matrix<Scalar,6,2>::Zero();  ///< Differential strain E = X_t - X_0
    Eigen::Matrix<Scalar, 6, 2> S   = Eigen::Matrix<Scalar,6,2>::Zero();  ///< Stress resultants S^alpha

    // Force and stiffness
    Eigen::Matrix<Scalar, DOF_PER_ELEM, 1>             Fint
        = Eigen::Matrix<Scalar,DOF_PER_ELEM,1>::Zero();        ///< Internal force vector (24)

    Eigen::Matrix<Scalar, DOF_PER_ELEM, DOF_PER_ELEM>  K
        = Eigen::Matrix<Scalar,DOF_PER_ELEM,DOF_PER_ELEM>::Zero(); ///< K_M + K_G (24x24)

    Eigen::Matrix<Scalar, DOF_PER_ELEM, DOF_PER_ELEM>  KM
        = Eigen::Matrix<Scalar,DOF_PER_ELEM,DOF_PER_ELEM>::Zero(); ///< Material stiffness only

    Eigen::Matrix<Scalar, DOF_PER_ELEM, DOF_PER_ELEM>  KG
        = Eigen::Matrix<Scalar,DOF_PER_ELEM,DOF_PER_ELEM>::Zero(); ///< Geometric stiffness only

    // Energy
    Scalar strainEnergy{0.0};  ///< W = 0.5 * sum_alpha <S^alpha, E_{t alpha}>

    void setZero() {
        Xt.setZero(); E.setZero(); S.setZero();
        Fint.setZero(); K.setZero(); KM.setZero(); KG.setZero();
        strainEnergy = 0.0;
    }
};

// --- Main assembler ---------------------------------------------------------

/**
 * @brief Compute element internal force and full tangent stiffness.
 *
 * Assembles F_int_e in R^{24} and K_e = K_M + K_G in R^{24x24} for one Q4
 * Cosserat shell element.
 *
 * Algorithm (same structure as ShellElementComputer, with corrections):
 *
 *   1. Strain twists X_t at centroid (with J_R^{-1}, from SE3ShellKinematics)
 *   2. Differential strain E = X_t - X_0
 *   3. Stress resultants S^alpha = sum_beta D^{alpha beta} E_{t beta}
 *      (uses ConstitutiveTensorsFull: isotropic/orthotropic + cylindrical)
 *   4. Strain energy W = 0.5 * sum_alpha <S^alpha, E_{t alpha}>
 *   5. Gauss integration loop (2x2):
 *      a. B_alpha at Gauss point
 *      b. F_int += sum_alpha B_alpha^T S^alpha j0 w
 *      c. K_M   += sum_{alpha,beta} B_alpha^T D^{alpha beta} B_beta j0 w
 *      d. K_G   += sum_alpha B_alpha^T ad_{S^alpha} B_alpha j0 w  [CORRECTED]
 *         optionally symmetrized: K_G <- 0.5*(K_G + K_G^T)
 *
 * Correction vs. ShellElementComputer:
 *   Simple version: K_G uses -ad_{S^alpha}^T (co-adjoint)
 *   This version:   K_G uses  ad_{S^alpha}   (correct adjoint, Eq. 80)
 *
 * @param nodes        Current node configurations (position + orientation)
 * @param elemRef      Pre-computed reference element data
 * @param tensors      Full constitutive tensors (isotropic/orthotropic + cyl)
 * @param[out] elem    Extended element state (Fint, K, KM, KG, strainEnergy)
 * @param symmetrize   If true (default), symmetrize K_G to ensure K_e is SPD
 *                     at reasonable configurations. Recommended for Newton solver.
 */
inline void computeElementForceAndStiffnessFull(
    const std::array<NodeConfig, NODES_PER_ELEM>& nodes,
    const ShellElementRef& elemRef,
    const ConstitutiveTensorsFull& tensors,
    ShellElementCurrentFull& elem,
    bool symmetrize = true)
{
    elem.setZero();

    // -- 1. Strain twists at centroid (J_R^{-1} already applied) -------------
    elem.Xt = computeStrainTwistsAtCentroid(nodes);

    // -- 2. Differential strain E = X_t - X_0 --------------------------------
    elem.E = elem.Xt - elemRef.X0;

    // -- 3. Stress resultants S^alpha ----------------------------------------
    elem.S = tensors.allStressResultants(elem.E);

    // -- 4. Strain energy W = 0.5 * sum_alpha <S^alpha, E_{t alpha}> ---------
    for (int alpha = 0; alpha < 2; ++alpha)
        elem.strainEnergy += 0.5 * elem.S.col(alpha).dot(elem.E.col(alpha));

    // -- 5. Gauss integration loop -------------------------------------------
    const auto gaussPts = gaussPoints2x2();

    for (int g = 0; g < N_GAUSS; ++g) {
        const double xi  = gaussPts[g].xi;
        const double eta = gaussPts[g].eta;
        const double w   = gaussPts[g].w;
        const double dA  = elemRef.j0Gauss[g] * w;

        // B_alpha at Gauss point g
        // Strains are constant from centroid (anti-locking), shape functions vary
        Eigen::Matrix<double, 6, DOF_PER_ELEM> B[2];
        for (int alpha = 0; alpha < 2; ++alpha)
            B[alpha] = computeBMatrix(elem.Xt.col(alpha), xi, eta, alpha);

        // -- F_int += sum_alpha B_alpha^T * S^alpha * dA ---------------------
        for (int alpha = 0; alpha < 2; ++alpha)
            elem.Fint += B[alpha].transpose() * elem.S.col(alpha) * dA;

        // -- K_M += sum_{alpha,beta} B_alpha^T * D^{alpha beta} * B_beta * dA
        for (int alpha = 0; alpha < 2; ++alpha) {
            for (int beta = 0; beta < 2; ++beta) {
                const Eigen::Matrix<double,6,6>& Dab = tensors.D[2*alpha + beta];
                elem.KM += B[alpha].transpose() * Dab * B[beta] * dA;
            }
        }

        // -- K_G += sum_alpha B_alpha^T * ad_{S^alpha} * B_alpha * dA --------
        //
        // CORRECTION vs. ShellElementComputer:
        //   Wrong (simple): -adMatrix(S)^T = adStarMatrix(S)   co-adjoint
        //   Correct (full):  adMatrix(S)                        adjoint
        //
        // Derivation: the geometric stiffness comes from
        //   d/de <S^alpha, K_alpha delta_eta>|_{e=0} * Delta_eta
        //   = <S^alpha, ad_{K_alpha Delta_eta} delta_eta>
        //   = (B_alpha Delta_eta)^T * ad_{S^alpha}^T * delta_eta
        //   (by Lie algebra duality: <a, [b,c]> = <ad_b^T a, c> = <-ad_a b, c>)
        //   -> contribution to tangent: B_alpha^T * ad_{S^alpha} * B_alpha
        //
        // Symmetrization: ad_{S} is generally not symmetric. We symmetrize
        // 0.5*(K_G + K_G^T) to enforce symmetry at all configurations
        // (the exact tangent is symmetric only at equilibrium).
        for (int alpha = 0; alpha < 2; ++alpha) {
            const Eigen::Matrix<double,6,6> adS = adMatrix(elem.S.col(alpha));
            Eigen::Matrix<double, DOF_PER_ELEM, DOF_PER_ELEM> KG_alpha =
                B[alpha].transpose() * adS * B[alpha] * dA;

            if (symmetrize)
                KG_alpha = 0.5 * (KG_alpha + KG_alpha.transpose());

            elem.KG += KG_alpha;
        }
    }

    // -- 6. Total tangent: K = K_M + K_G -------------------------------------
    elem.K = elem.KM + elem.KG;
}

// --- Reference element initializer (shared with simple version) -------------

/**
 * @brief Initialize reference element data.
 *
 * Identical to initShellElementRef in ShellElementComputer.h.
 * Provided here so ShellElementComputerFull.h is self-contained.
 *
 * Computes:
 *   X_0      = strain twists at reference centroid (with J_R^{-1})
 *   X_0*     = (X_0^T X_0)^{-1} X_0^T   (left pseudo-inverse, 2x6)
 *   j_0      = area Jacobian at centroid
 *   j_0^g    = area Jacobian at each Gauss point
 *
 * @param refNodes   Reference node configurations
 * @param[out] elem  Initialized reference element data (ShellElementRef)
 */
inline void initShellElementRefFull(
    const std::array<NodeConfig, NODES_PER_ELEM>& refNodes,
    ShellElementRef& elem)
{
    // Reference strains at centroid
    elem.X0 = computeStrainTwistsAtCentroid(refNodes);

    // Pseudo-inverse X_0* = (X_0^T X_0)^{-1} X_0^T  in R^{2x6}
    // X_0 in R^{6x2} -> X_0^T X_0 is 2x2 (always invertible for valid element)
    const Eigen::Matrix<double,2,2> XtX = elem.X0.transpose() * elem.X0;
    elem.X0star = XtX.ldlt().solve(elem.X0.transpose());

    // Area Jacobians
    std::array<Eigen::Vector3d, NODES_PER_ELEM> refPos;
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        refPos[i] = refNodes[i].position;

    elem.j0 = computeAreaJacobian(refPos, 0.0, 0.0);

    const auto gaussPts = gaussPoints2x2();
    for (int g = 0; g < N_GAUSS; ++g)
        elem.j0Gauss[g] = computeAreaJacobian(refPos, gaussPts[g].xi, gaussPts[g].eta);
}

// --- SE(3) Newton update (identical to simple version) ----------------------

/**
 * @brief Apply a Newton increment to a node using SE(3) multiplicative update.
 *
 *   R_new = R_old * exp(delta_omega)   -- right multiplication on SO(3)
 *   phi_new = phi_old + delta_v
 *
 * The increment is [delta_omega (3); delta_v (3)] in R^6.
 *
 * @param node      Node configuration (modified in place)
 * @param increment 6D increment [delta_omega; delta_v]
 */
inline void applyNodeUpdateFull(NodeConfig& node,
                                const Eigen::Matrix<double,6,1>& increment)
{
    node.orientation = node.orientation * SO3d::exp(increment.head<3>());
    node.position   += increment.tail<3>();
}

// --- Convenience: comparison helper -----------------------------------------

/**
 * @brief Compute ratio ||K_G|| / ||K_M|| for one element.
 *
 * Useful to assess when the geometric stiffness becomes significant.
 * Rule of thumb:
 *   < 0.01  -> geometric effects negligible (small-deformation regime)
 *   > 0.10  -> geometric stiffness important (large-deformation regime)
 *   > 1.00  -> possible buckling / snap-through
 *
 * @param elem  Computed element state (after computeElementForceAndStiffnessFull)
 * @return      ||K_G||_F / (||K_M||_F + epsilon)
 */
inline double geometricStiffnessRatio(const ShellElementCurrentFull& elem)
{
    const double normKM = elem.KM.norm();
    const double normKG = elem.KG.norm();
    return normKG / (normKM + 1e-300);
}

} // namespace sofa::component::cosserat::shell
