/******************************************************************************
 * Cosserat Shell -- SE(3) Kinematics for one Q4 element                     *
 *                                                                            *
 * Computes:                                                                  *
 *   - Right Jacobian inverse J_R^{-1}(omega) on SO(3)                       *
 *   - Strain twists xi_{t alpha} at element centroid (shear-locking free)   *
 *     with J_R^{-1} for accuracy at large inter-node rotations              *
 *   - Differential strains E = X_t - X_0                                    *
 *   - Local tangent operators K_alpha (for weak form discretization)        *
 *   - Lie adjoint operators ad_{xi_{t alpha}} (for geometric stiffness)     *
 *                                                                            *
 * Key addition vs. simple version:                                           *
 *   Without J_R^{-1}: phi_{t alpha} ~ sum_i dN^i/dX_alpha * log(R_c^-1 Ri) *
 *   With    J_R^{-1}: phi_{t alpha} = J_R^{-1}(omega_bar) * above sum       *
 *   Error without J_R^{-1} is O(|omega_bar|^2) -- critical in cylindrical   *
 *   coords where inter-node rotations ~= 2pi/N_theta can reach 45 deg for   *
 *   N_theta = 8 elements.                                                    *
 *                                                                            *
 * Reference: cosserat_shell_formulation.md S9, formulation_review.md S4     *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                   *
 ******************************************************************************/
#pragma once

#include "ShellElement.h"
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>
#include <Eigen/Dense>
#include <cmath>

namespace sofa::component::cosserat::shell {

using SE3d = sofa::component::cosserat::liegroups::SE3<double>;
using SO3d = sofa::component::cosserat::liegroups::SO3<double>;

// --- Helper: skew-symmetric matrix (hat operator) ---------------------------

/**
 * @brief Returns the 3x3 skew-symmetric matrix [v]_x.
 */
inline Eigen::Matrix3d skew(const Eigen::Vector3d& v) {
    Eigen::Matrix3d S;
    S <<  0.0, -v(2),  v(1),
         v(2),   0.0, -v(0),
        -v(1),  v(0),   0.0;
    return S;
}

// --- Right Jacobian inverse J_R^{-1}(omega) on SO(3) ------------------------

/**
 * @brief Compute the right Jacobian inverse J_R^{-1}(omega) in R^{3x3}.
 *
 * This maps the derivative of the log map back to the Lie algebra:
 *
 *   J_R^{-1}(omega) = I + (1/2)[omega]_x + c(theta)[omega]_x^2
 *
 * where theta = ||omega|| and:
 *
 *   c(theta) = 1/theta^2 - (1 + cos(theta)) / (2 * theta * sin(theta))
 *
 * Near-zero Taylor expansion (avoids division by zero for theta -> 0):
 *   c(theta) = 1/12 + theta^2/720 + O(theta^4)
 *
 * Role in kinematics:
 *   Without J_R^{-1}, the angular strain interpolation
 *       phi_{t alpha} ~= sum_i dN^i/dX_alpha * log(R_base^{-1} R_i)
 *   carries an O(|omega_bar|^2) error where omega_bar is the mean log
 *   at the centroid.  The corrected formula is:
 *       phi_{t alpha} = J_R^{-1}(omega_bar) * sum_i dN^i/dX_alpha * log(R_base^{-1} R_i)
 *
 * This correction is critical for cylindrical shells (inter-node rotations
 * ~= 2pi/N_theta can be 45 deg or more) and negligible for flat plates with
 * small rotations (< ~15 deg).
 *
 * @param omega   Mean rotation vector omega_bar in R^3
 * @return        J_R^{-1}(omega_bar) in R^{3x3}
 */
inline Eigen::Matrix3d rightJacobianInverse(const Eigen::Vector3d& omega)
{
    const double theta2 = omega.squaredNorm();

    // Regular branch: theta > epsilon
    constexpr double kEps2 = 1e-6;  // threshold on theta^2 (~= theta > 1e-3 rad)
    if (theta2 > kEps2) {
        const double theta = std::sqrt(theta2);
        const double c = 1.0 / theta2
                       - (1.0 + std::cos(theta)) / (2.0 * theta * std::sin(theta));
        const Eigen::Matrix3d ox = skew(omega);
        return Eigen::Matrix3d::Identity()
             + 0.5 * ox
             + c   * ox * ox;
    }

    // Near-zero branch: Taylor expansion c(theta) -> 1/12 + theta^2/720
    {
        const double c = 1.0/12.0 + theta2/720.0;
        const Eigen::Matrix3d ox = skew(omega);
        return Eigen::Matrix3d::Identity()
             + 0.5 * ox
             + c   * ox * ox;
    }
}

// --- Helper: 6x6 ad_xi matrix (Lie bracket on se(3)) ------------------------

/**
 * @brief Returns the 6x6 adjoint (small) matrix ad_xi for xi = [omega; v] in R^6.
 *
 * Convention: [angular (head<3>); translational (tail<3>)] = [omega; v]
 *
 *   ad_xi = | [omega]_x    0      |
 *           | [v]_x     [omega]_x |
 *
 * so that ad_xi * eta = [omega x omega_eta ; omega x v_eta + v x omega_eta].
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
 * @brief Returns the transpose of ad_xi (co-adjoint ad*_xi).
 *
 * ad*_xi = -ad_xi^T
 * Used in the geometric stiffness contribution.
 */
inline Eigen::Matrix<double,6,6> adStarMatrix(const Eigen::Matrix<double,6,1>& xi) {
    return -adMatrix(xi).transpose();
}

// --- Node configuration -----------------------------------------------------

/**
 * @brief Configuration of one shell node (position + orientation).
 */
struct NodeConfig {
    Eigen::Vector3d position    = Eigen::Vector3d::Zero();
    SO3d            orientation = SO3d::identity();
};

// --- Centroid strain twist computation --------------------------------------

/**
 * @brief Compute strain twists xi_{t alpha} at element centroid (x=0, y=0).
 *
 * xi_{t alpha} = (g_t^{-1} dg_t/dX_alpha)^vee  in R^6    (Eq. 17)
 *
 * Discretization on the Q4 element:
 *
 *   Step 1 -- per-node logs in the tangent space of R_base (= R_0):
 *     log_i = log(R_base^{-1} * R_i)  in R^3
 *
 *   Step 2 -- mean log at centroid (interpolated orientation):
 *     omega_bar = sum_i N^i(0,0) * log_i   [= 0.25 * sum for Q4]
 *
 *   Step 3 -- centroid orientation:
 *     R_c = R_base * exp(omega_bar)
 *
 *   Step 4 -- right Jacobian inverse J_R^{-1}(omega_bar):
 *     Corrects the angular strain for large inter-node rotations.
 *     Identity-like for small omega_bar; non-trivial for cylindrical shells.
 *
 *   Step 5 -- strain twists for each direction alpha:
 *
 *     Translational part (indices 3-5):
 *       rho_{t alpha} = R_c^{-1} * sum_i dN^i/dX_alpha * phi_i
 *
 *     Angular part (indices 0-2) WITH J_R^{-1}:
 *       delta_omega_alpha = sum_i dN^i/dX_alpha * log_i
 *       phi_{t alpha}     = J_R^{-1}(omega_bar) * delta_omega_alpha
 *
 * Returns X_t = [xi_{t1} | xi_{t2}] in R^{6x2}.
 *
 * @param nodes   Array of 4 node configurations
 * @return        Strain matrix X_t in R^{6x2}
 */
inline Eigen::Matrix<double,6,2> computeStrainTwistsAtCentroid(
    const std::array<NodeConfig, NODES_PER_ELEM>& nodes)
{
    // Step 1: per-node logs in the tangent space of R_base
    const SO3d& R_base = nodes[0].orientation;

    std::array<Eigen::Vector3d, NODES_PER_ELEM> log_i;
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        log_i[i] = (R_base.inverse() * nodes[i].orientation).log();

    // Step 2: mean log at centroid -- omega_bar = sum_i N^i(0,0) * log_i
    // N^i(0,0) = 0.25 for all i in a Q4 element
    Eigen::Vector3d omega_bar = Eigen::Vector3d::Zero();
    for (int i = 0; i < NODES_PER_ELEM; ++i)
        omega_bar += shapeFunction(i, 0.0, 0.0) * log_i[i];

    // Step 3: centroid orientation R_c = R_base * exp(omega_bar)
    const SO3d R_c = R_base * SO3d::exp(omega_bar);

    // Step 4: right Jacobian inverse J_R^{-1}(omega_bar)
    // This corrects angular strains for large inter-node rotations.
    // For flat plates with small rotations: J_R^{-1} ~= I (no-op).
    // For cylindrical shells: significant correction (see header).
    const Eigen::Matrix3d JRinv = rightJacobianInverse(omega_bar);

    // Step 5: strain twists xi_{t alpha} for each surface direction
    Eigen::Matrix<double,6,2> Xt = Eigen::Matrix<double,6,2>::Zero();

    for (int alpha = 0; alpha < 2; ++alpha) {
        Eigen::Vector3d dphi_dalpha = Eigen::Vector3d::Zero();   // sum_i dN^i/dX_alpha * phi_i
        Eigen::Vector3d dlog_dalpha = Eigen::Vector3d::Zero();   // sum_i dN^i/dX_alpha * log_i

        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            const double dNi_alpha = shapeFunctionGrad(i, 0.0, 0.0)(alpha);
            dphi_dalpha += dNi_alpha * nodes[i].position;
            dlog_dalpha += dNi_alpha * log_i[i];
        }

        // Angular part: phi_{t alpha} = J_R^{-1}(omega_bar) * delta_omega_alpha
        // Without J_R^{-1} this would be just dlog_dalpha (first-order approximation)
        Xt.col(alpha).head<3>() = JRinv * dlog_dalpha;

        // Translational part: rho_{t alpha} = R_c^{-1} * d(phi)/dX_alpha
        Xt.col(alpha).tail<3>() = R_c.inverse().act(dphi_dalpha);
    }

    return Xt;
}

// --- Reference area Jacobian ------------------------------------------------

/**
 * @brief Compute the reference area Jacobian j_0 at a given point (x,y).
 *
 * j_0 = || d(phi_0)/dX_1  x  d(phi_0)/dX_2 ||
 *
 * @param refPositions  Reference node positions in R^3 x 4
 * @param x, y          Reference coordinates in [-1,1]
 * @return              Area Jacobian (scalar, > 0 for a valid element)
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

// --- B-matrix (discrete tangent operator K_alpha) ---------------------------

/**
 * @brief Compute the discrete tangent operator B_alpha in R^{6 x DOF_PER_ELEM}.
 *
 * B_alpha discretizes the operator K_alpha eta = d(eta)/dX_alpha + ad_{xi_{t alpha}} eta
 * (Lemma 1 of the paper). It maps the element DOF vector eta_e in R^{24} to
 * K_alpha eta in R^6.
 *
 * For a Q4 element with shape functions N^i:
 *   K_alpha eta ~= sum_i dN^i/dX_alpha * eta^i + ad_{xi_{t alpha}} sum_i N^i * eta^i
 *               = sum_i [dN^i/dX_alpha * I_6 + N^i * ad_{xi_{t alpha}}] * eta^i
 *
 * Therefore:
 *   B_alpha[:, 6i:6i+6] = dN^i/dX_alpha * I_6 + N^i(x,y) * ad_{xi_{t alpha}}
 *
 * @param xi_alpha    Current strain twist xi_{t alpha} in R^6 at evaluation point
 * @param x, y        Reference coordinates of the evaluation point
 * @param alpha       Direction index: 0 for X_1, 1 for X_2
 * @return            B matrix in R^{6 x 24}
 */
inline Eigen::Matrix<double, 6, DOF_PER_ELEM> computeBMatrix(
    const Eigen::Matrix<double,6,1>& xi_alpha,
    double x, double y,
    int alpha)
{
    Eigen::Matrix<double, 6, DOF_PER_ELEM> B =
        Eigen::Matrix<double, 6, DOF_PER_ELEM>::Zero();

    const Eigen::Matrix<double,6,6> ad_xi = adMatrix(xi_alpha);

    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const double Ni        = shapeFunction(i, x, y);
        const double dNi_alpha = shapeFunctionGrad(i, x, y)(alpha);

        // Block for node i: [dN^i/dX_alpha * I_6 + N^i * ad_{xi_{t alpha}}]
        B.template block<6,6>(0, 6*i) =
            dNi_alpha * Eigen::Matrix<double,6,6>::Identity()
            + Ni * ad_xi;
    }

    return B;
}

} // namespace sofa::component::cosserat::shell
