/******************************************************************************
 * Cosserat Shell — SE(3) Kinematics for one Q4 element                      *
 *                                                                            *
 * Computes:                                                                  *
 *   - Strain twists ζ_{tα} at element centroid (shear-locking free)         *
 *   - Differential strains E = X_t - X_0                                    *
 *   - Local tangent operators K_α (for weak form discretization)            *
 *   - Lie adjoint operators ad_{ζ_{tα}} (for geometric stiffness)           *
 *                                                                            *
 * Reference: cosserat_shell_v0.pdf, Sections 3–5                           *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                   *
 ******************************************************************************/
#pragma once

#include "ShellElement.h"
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>
#include <Eigen/Dense>

namespace sofa::component::cosserat::shell {

using SE3d = sofa::component::cosserat::liegroups::SE3<double>;
using SO3d = sofa::component::cosserat::liegroups::SO3<double>;

// ─── Helper: skew-symmetric matrix (hat operator) ─────────────────────────

/**
 * @brief Returns the 3×3 skew-symmetric matrix [v]_×.
 */
inline Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d S;
    S <<  0.0, -v(2),  v(1),
         v(2),   0.0, -v(0),
        -v(1),  v(0),   0.0;
    return S;
}

// ─── Helper: 6×6 ad_xi matrix (Lie bracket on se(3)) ─────────────────────

/**
 * @brief Returns the 6×6 adjoint (small) matrix ad_ξ for ξ = [ω; v] ∈ ℝ⁶.
 *
 * ad_ξ = | [ω]×   0    |    so that ad_ξ η = [ω×ω_η ; ω×v_η + v×ω_η]
 *         | [v]×  [ω]×  |
 *
 * Note: this follows the convention [angular; translational] for se(3) elements.
 */
inline Eigen::Matrix<double,6,6> adMatrix(const Eigen::Matrix<double,6,1>& xi) {
    const Eigen::Vector3d omega = xi.head<3>();
    const Eigen::Vector3d v     = xi.tail<3>();
    Eigen::Matrix<double,6,6> ad = Eigen::Matrix<double,6,6>::Zero();
    ad.topLeftCorner<3,3>()     = skew(omega);
    ad.bottomRightCorner<3,3>() = skew(omega);
    ad.bottomLeftCorner<3,3>()  = skew(v);
    return ad;
}

/**
 * @brief Returns the transpose of ad_ξ (co-adjoint ad*_ξ).
 *
 * ad*_ξ · f = - ad_ξ^T · f
 * Used in the geometric stiffness contribution.
 */
inline Eigen::Matrix<double,6,6> adStarMatrix(const Eigen::Matrix<double,6,1>& xi) {
    return -adMatrix(xi).transpose();
}

// ─── Node configuration ────────────────────────────────────────────────────

/**
 * @brief Configuration of one shell node (position + orientation).
 */
struct NodeConfig {
    Eigen::Vector3d position    = Eigen::Vector3d::Zero();
    SO3d            orientation = SO3d::identity();
};

// ─── Centroid strain twist computation ────────────────────────────────────

/**
 * @brief Compute strain twists ζ_{tα} at the element centroid (x=0, y=0).
 *
 * ζ_{tα} = (g_t^{-1} ∂g_t/∂ξ^α)^∨  ∈ ℝ⁶    (Eq. 17)
 *
 * Discretized using bilinear shape functions:
 *   φ_t(0,0)         = Σ_i N^i(0,0) φ^i = ¼ Σ_i φ^i
 *   ∂φ_t/∂ξ^α(0,0)  = Σ_i ∂N^i/∂ξ^α(0,0) · φ^i
 *   ∂R_t/∂ξ^α(0,0)  = R_t(0,0) · [Σ_i ∂N^i/∂ξ^α(0,0) · log(R_t(0,0)^{-1} R^i)]
 *                      (linearized tangent interpolation)
 *
 * Returns X_t = [ζ_{t1} | ζ_{t2}] ∈ ℝ^{6×2}.
 *
 * @param nodes   Array of 4 node configurations
 * @return        Strain matrix X_t ∈ ℝ^{6×2}
 */
inline Eigen::Matrix<double,6,2> computeStrainTwistsAtCentroid(
    const std::array<NodeConfig, NODES_PER_ELEM>& nodes)
{
    // ── Interpolate position and orientation at centroid (x=0, y=0) ──────
    Eigen::Vector3d phi_c = Eigen::Vector3d::Zero();
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const double Ni = shapeFunction(i, 0.0, 0.0);  // = 0.25 for all i
        phi_c += Ni * nodes[i].position;
    }

    // For orientation at centroid, use a linearized average on SO(3):
    //   R_c = R_0 * exp( Σ_i N^i * log(R_0^{-1} * R_i) )
    // where R_0 is the reference rotation (first node, or any fixed frame).
    // Simpler for Q4: use node 0 as base and interpolate in its tangent space.
    const SO3d& R_base = nodes[0].orientation;
    Eigen::Vector3d avg_log = Eigen::Vector3d::Zero();
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const double Ni = shapeFunction(i, 0.0, 0.0);
        // log( R_base^{-1} * R_i ) ∈ ℝ³
        avg_log += Ni * (R_base.inverse() * nodes[i].orientation).log();
    }
    const SO3d R_c = R_base * SO3d::exp(avg_log);

    // ── Compute ∂φ/∂ξ^α and ∂R/∂ξ^α at centroid ─────────────────────────
    Eigen::Matrix<double,6,2> Xt = Eigen::Matrix<double,6,2>::Zero();

    for (int alpha = 0; alpha < 2; ++alpha) {
        Eigen::Vector3d dphi_dxi = Eigen::Vector3d::Zero();
        Eigen::Vector3d domega_dxi = Eigen::Vector3d::Zero();

        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            const Eigen::Vector2d dNi = shapeFunctionGrad(i, 0.0, 0.0);
            const double dNi_alpha = dNi(alpha);

            // Translational part: ∂φ/∂ξ^α
            dphi_dxi += dNi_alpha * nodes[i].position;

            // Angular part: interpolate angular velocity via log on so(3)
            // ∂R/∂ξ^α ≈ R_c * [dNi_alpha * log(R_c^{-1} R_i)]_×
            // → angular strain: R_c^T ∂R/∂ξ^α = [Σ_i dNi_alpha * log(R_c^{-1} R_i)]_×
            domega_dxi += dNi_alpha * (R_base.inverse() * nodes[i].orientation).log();
        }

        // ζ_{tα} = [R_c^T ∂R/∂ξ^α ; R_c^T ∂φ/∂ξ^α]
        // Angular part: the log gives us ω ∈ ℝ³ such that R_c^T ∂R/∂ξ^α = [ω]×
        Xt.col(alpha).head<3>() = domega_dxi;
        Xt.col(alpha).tail<3>() = R_c.inverse().act(dphi_dxi);
    }

    return Xt;
}

// ─── Reference area Jacobian ───────────────────────────────────────────────

/**
 * @brief Compute the reference area Jacobian j_0 at a given point (x,y).
 *
 * j_0 = || ∂φ_0/∂ξ¹ × ∂φ_0/∂ξ² ||
 *
 * @param refNodes  Reference node positions ∈ ℝ³ × 4
 * @param x, y      Reference coordinates ∈ [-1,1]
 * @return          Area Jacobian (scalar)
 */
inline double computeAreaJacobian(
    const std::array<Eigen::Vector3d, NODES_PER_ELEM>& refPositions,
    double x, double y)
{
    Eigen::Vector3d dphi1 = Eigen::Vector3d::Zero();
    Eigen::Vector3d dphi2 = Eigen::Vector3d::Zero();

    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const Eigen::Vector2d dNi = shapeFunctionGrad(i, x, y);
        dphi1 += dNi(0) * refPositions[i];
        dphi2 += dNi(1) * refPositions[i];
    }

    return dphi1.cross(dphi2).norm();
}

// ─── B-matrix (discrete tangent operator K_α) ─────────────────────────────

/**
 * @brief Compute the discrete tangent operator B_alpha ∈ ℝ^{6 × DOF_PER_ELEM}.
 *
 * B_alpha discretizes the operator K_α η = ∂η/∂ξ^α + ad_{ζ_{tα}} η (Lemma 1).
 * It maps the element DOF vector η_e ∈ ℝ^{24} to K_α η ∈ ℝ⁶.
 *
 * For a Q4 element with shape functions N^i:
 *   K_α η ≈ Σ_i ∂N^i/∂ξ^α · η^i + ad_{ζ_{tα}} Σ_i N^i · η^i
 *          = Σ_i [∂N^i/∂ξ^α · I₆ + N^i · ad_{ζ_{tα}}] · η^i
 *
 * Therefore: B_alpha[:, 6i:6i+6] = ∂N^i/∂ξ^α · I₆ + N^i(x,y) · ad_{ζ_{tα}}
 *
 * @param zeta_alpha  Current strain twist ζ_{tα} ∈ ℝ⁶ at evaluation point
 * @param x, y        Reference coordinates of evaluation point
 * @param alpha       Direction index α ∈ {0,1}
 * @return            B matrix ∈ ℝ^{6 × 24}
 */
inline Eigen::Matrix<double, 6, DOF_PER_ELEM> computeBMatrix(
    const Eigen::Matrix<double,6,1>& zeta_alpha,
    double x, double y,
    int alpha)
{
    Eigen::Matrix<double, 6, DOF_PER_ELEM> B = Eigen::Matrix<double,6,DOF_PER_ELEM>::Zero();

    const Eigen::Matrix<double,6,6> ad_zeta = adMatrix(zeta_alpha);

    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const double Ni = shapeFunction(i, x, y);
        const double dNi_alpha = shapeFunctionGrad(i, x, y)(alpha);

        // Block for node i: [∂N^i/∂ξ^α · I₆ + N^i · ad_{ζ_{tα}}]
        B.template block<6,6>(0, 6*i) =
            dNi_alpha * Eigen::Matrix<double,6,6>::Identity()
            + Ni * ad_zeta;
    }

    return B;
}

} // namespace sofa::component::cosserat::shell
