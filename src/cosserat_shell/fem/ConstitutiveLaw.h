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
 * Strain/stress 6-vector layout: [φ_x, φ_y, φ_z, ρ_x, ρ_y, ρ_z]
 *   Angular (0-2): torsion, bending-Y, bending-Z (curvature components)
 *   Linear  (3-5): elongation-X, shear-Y, transverse-shear-Z
 *
 * Physical mapping for a flat plate (X_1 = first, X_2 = second direction):
 *   E_{t1}[0]=κ_11  E_{t1}[3]=ε_11  E_{t1}[4]=ε_12  E_{t1}[5]=γ_1 (transverse)
 *   E_{t2}[1]=κ_22  E_{t2}[3]=ε_21  E_{t2}[4]=ε_22  E_{t2}[5]=γ_2 (transverse)
 *
 * Non-zero blocks:
 *   D^{11}: [0,0]=D_b(κ_11²), [1,1]=C_b(twist), [2,2]=C_b, [3,3]=A_m, [4,4]=C_m, [5,5]=ksGh
 *   D^{22}: [0,0]=C_b(twist), [1,1]=D_b(κ_22²), [2,2]=C_b, [3,3]=C_m, [4,4]=A_m, [5,5]=ksGh
 *   D^{12}: [0,1]=νD_b (Poisson bending), [3,4]=νA_m (Poisson membrane)
 *   D^{21}: [1,0]=νD_b,                   [4,3]=νA_m
 *
 * NOTE — known approximation: the off-diagonal shear coupling
 *   D^{12}[4,3]=C_m and D^{12}[2,2]=C_b (symmetric shear split)
 *   is present in ConstitutiveLawFull but omitted here for simplicity.
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
    ConstitutiveTensors tensors;

    // ── D^{11} : strains in direction X_1 ────────────────────────────────
    // κ_11 → M_11 (bending), twist, ε_11 → N_11 (membrane), ε_12, γ_1
    {
        Eigen::Matrix<Scalar,6,6>& D = tensors.D[0];
        D.setZero();
        D(0, 0) = factorB;                          // κ_11² (bending)
        D(1, 1) = factorB * (1.0 - nu) / 2.0;      // twist from X_1 (C_b)
        D(2, 2) = factorB * (1.0 - nu) / 2.0;      // twist cross-term
        D(3, 3) = factor;                           // ε_11² (membrane)
        D(4, 4) = factor * (1.0 - nu) / 2.0;       // ε_12² (in-plane shear, C_m)
        D(5, 5) = ksGh;                             // γ_1² (transverse shear) ← FIX: was [4,4]+[5,5]
    }

    // ── D^{22} : strains in direction X_2 ────────────────────────────────
    // κ_22 → M_22 (bending), twist, ε_22 → N_22 (membrane), ε_21, γ_2
    {
        Eigen::Matrix<Scalar,6,6>& D = tensors.D[3];
        D.setZero();
        D(0, 0) = factorB * (1.0 - nu) / 2.0;      // twist cross-term (C_b)
        D(1, 1) = factorB;                          // κ_22² (bending)
        D(2, 2) = factorB * (1.0 - nu) / 2.0;      // twist from X_2
        D(3, 3) = factor * (1.0 - nu) / 2.0;       // ε_21² (in-plane shear, C_m)
        D(4, 4) = factor;                           // ε_22² (membrane)
        D(5, 5) = ksGh;                             // γ_2² (transverse shear)
    }

    // ── D^{12} : Poisson coupling X_2→S^1 ────────────────────────────────
    // FIX: was zero — adds missing Poisson effect between the two directions
    {
        Eigen::Matrix<Scalar,6,6>& D = tensors.D[1];
        D.setZero();
        D(0, 1) = nu * factorB;    // M_11 += ν·D_b·κ_22  (Poisson bending)
        D(3, 4) = nu * factor;     // N_11 += ν·A_m·ε_22  (Poisson membrane)
    }

    // ── D^{21} : Poisson coupling X_1→S^2  (transpose of D^{12}) ─────────
    {
        Eigen::Matrix<Scalar,6,6>& D = tensors.D[2];
        D.setZero();
        D(1, 0) = nu * factorB;    // M_22 += ν·D_b·κ_11
        D(4, 3) = nu * factor;     // N_22 += ν·A_m·ε_11
    }

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
