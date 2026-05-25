/******************************************************************************
 * Cosserat Shell — Constitutive Law (isotropic linear elasticity)            *
 *                                                                            *
 * Assembles the 6×6 material stiffness matrices D^{αβ} for an isotropic    *
 * linear-elastic Cosserat shell (Reissner–Mindlin-type with Cosserat        *
 * structure).                                                                *
 *                                                                            *
 * Reference: cosserat_shell_v0.pdf, Eqs. 40-41                             *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                   *
 ******************************************************************************/
#pragma once

#include <Eigen/Dense>
#include <array>

namespace sofa::component::cosserat::shell {

using Scalar = double;

/**
 * @brief Material parameters for an isotropic linear-elastic Cosserat shell.
 */
struct ShellMaterialParams {
    Scalar E{1.0};       ///< Young's modulus [Pa]
    Scalar nu{0.3};      ///< Poisson's ratio [-]
    Scalar h{0.001};     ///< Thickness [m]
    Scalar kappa{5.0/6.0}; ///< Shear correction factor [-]
};

/**
 * @brief Pre-computed constitutive matrices for the shell.
 *
 * Stores D^{αβ} ∈ ℝ^{6×6} for α,β ∈ {1,2}.
 * Convention: D[0] = D^{11}, D[1] = D^{12}, D[2] = D^{21}, D[3] = D^{22}
 *
 * Each D^{αβ} is block-diagonal:
 *   D^{αβ} = block-diag(D_bend^{αβ}, D_memb^{αβ})
 *
 * The 6-vector strain/stress layout is [ω_x, ω_y, ω_z, v_x, v_y, v_z],
 * i.e. angular (curvature/torsion) DOF first, then translational (extension/shear).
 */
struct ConstitutiveTensors {
    /// D^{αβ} matrices: index = 2*(α-1) + (β-1), so D11=D[0], D12=D[1], D21=D[2], D22=D[3]
    std::array<Eigen::Matrix<Scalar, 6, 6>, 4> D;

    ConstitutiveTensors() {
        for (auto& m : D) m.setZero();
    }
};

/**
 * @brief Compute the 4 constitutive matrices D^{αβ} for an isotropic shell.
 *
 * For an isotropic material:
 *   - D^{11} = D^{22} = block-diag(C_bend, C_memb)
 *   - D^{12} = D^{21} = block-diag(C_bend_12, C_memb_12)
 *
 * Membrane stiffness (3×3 in-plane):
 *   C_memb = Eh/(1-ν²) * | 1   ν   0          |
 *                         | ν   1   0          |
 *                         | 0   0   (1-ν)/2    |
 *
 * Bending stiffness (3×3 curvature):
 *   C_bend = Eh³/(12(1-ν²)) * same structure
 *
 * Shear transverse (added to translational diagonal, DOF [3] and [4]):
 *   k_s = G * κ * h   (appended to the diagonal blocks)
 *
 * @param params  Material parameters
 * @return        ConstitutiveTensors with all D^{αβ} filled
 */
inline ConstitutiveTensors computeConstitutiveTensors(const ShellMaterialParams& params)
{
    const Scalar E  = params.E;
    const Scalar nu = params.nu;
    const Scalar h  = params.h;
    const Scalar ks = params.kappa;

    const Scalar G      = E / (2.0 * (1.0 + nu));
    const Scalar factor = E * h / (1.0 - nu * nu);
    const Scalar factorB = factor * h * h / 12.0;  // = Eh³/[12(1-ν²)]
    const Scalar ksGh   = ks * G * h;              // shear transverse stiffness

    // ── Membrane stiffness 3×3 ──
    Eigen::Matrix<Scalar,3,3> Cm = Eigen::Matrix<Scalar,3,3>::Zero();
    Cm(0,0) = factor;          Cm(0,1) = nu * factor;
    Cm(1,0) = nu * factor;     Cm(1,1) = factor;
    Cm(2,2) = factor * (1.0 - nu) / 2.0;

    // ── Bending stiffness 3×3 ──
    Eigen::Matrix<Scalar,3,3> Cb = Eigen::Matrix<Scalar,3,3>::Zero();
    Cb(0,0) = factorB;         Cb(0,1) = nu * factorB;
    Cb(1,0) = nu * factorB;    Cb(1,1) = factorB;
    Cb(2,2) = factorB * (1.0 - nu) / 2.0;

    // ── Build D^{αβ} ─────────────────────────────────────────────────────
    // Layout: [angular (bend/torsion): 0-2 | translational (memb+transv-shear): 3-5]
    // For an isotropic shell: only diagonal blocks D^{11} = D^{22} are non-zero
    // (D^{12} = D^{21} = 0 for isotropic in the principal material axes).
    //
    // D^{αα} = | Cb   0  |   with transverse shear on diag [3,4]
    //           |  0   Cm |
    //
    // D^{αβ} (α≠β) = 0 (isotropic, no coupling between ξ¹ and ξ² directions)

    ConstitutiveTensors tensors;

    // D^{11}  (index 0)
    Eigen::Matrix<Scalar,6,6>& D11 = tensors.D[0];
    D11.setZero();
    D11.template block<3,3>(0,0) = Cb;   // bending
    D11.template block<3,3>(3,3) = Cm;   // membrane
    // Transverse shear correction on v_y, v_z (translational DOF 4 and 5)
    D11(4,4) += ksGh;
    D11(5,5) += ksGh;

    // D^{22}  (index 3) = same structure for isotropic
    tensors.D[3] = D11;

    // D^{12} = D^{21} = 0 (isotropic)
    tensors.D[1].setZero();
    tensors.D[2].setZero();

    return tensors;
}

/**
 * @brief Compute stress resultant S^α from strain E_{tβ} via D^{αβ}.
 *
 * S^α = Σ_{β=1}^{2} D^{αβ} · E_{tβ}   (Eq. 39 of the paper)
 *
 * @param tensors   Pre-computed constitutive tensors
 * @param E         Differential strain matrix E ∈ ℝ^{6×2}, column β = E_{tβ}
 * @param alpha     Index α ∈ {0,1} (for α=1 or α=2 in 1-based notation)
 * @return          Stress resultant S^α ∈ ℝ⁶
 */
inline Eigen::Matrix<Scalar,6,1> computeStressResultant(
    const ConstitutiveTensors& tensors,
    const Eigen::Matrix<Scalar,6,2>& E,
    int alpha)
{
    Eigen::Matrix<Scalar,6,1> S = Eigen::Matrix<Scalar,6,1>::Zero();
    for (int beta = 0; beta < 2; ++beta) {
        // D^{alpha+1, beta+1} is stored at index 2*alpha + beta
        S += tensors.D[2 * alpha + beta] * E.col(beta);
    }
    return S;
}

/**
 * @brief Compute all stress resultants S^1 and S^2.
 *
 * @param tensors   Constitutive tensors
 * @param E         Differential strain ∈ ℝ^{6×2}
 * @return          Matrix S ∈ ℝ^{6×2}, column α = S^{α+1}
 */
inline Eigen::Matrix<Scalar,6,2> computeAllStressResultants(
    const ConstitutiveTensors& tensors,
    const Eigen::Matrix<Scalar,6,2>& E)
{
    Eigen::Matrix<Scalar,6,2> S;
    S.col(0) = computeStressResultant(tensors, E, 0);
    S.col(1) = computeStressResultant(tensors, E, 1);
    return S;
}

} // namespace sofa::component::cosserat::shell
