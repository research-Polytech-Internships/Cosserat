/******************************************************************************
 * Cosserat Shell — Full constitutive law                                     *
 *                                                                            *
 * Extends ConstitutiveLaw.h with:                                            *
 *   1. Physically complete D^{αβ} (Poisson + symmetric shear split)         *
 *   2. Orthotropic materials (E1, E2, G12, G13, G23, ν12)                   *
 *   3. Cylindrical coordinate scaling (radius R)                             *
 *                                                                            *
 * Comparison:                                                                *
 *   ConstitutiveLaw.h     — isotropic, flat, Poisson minimal fix             *
 *   ConstitutiveLawFull.h — isotropic + ortho + cyl, complete shear split    *
 *                                                                            *
 * Reference: cosserat_shell_formulation.md §7, formulation_cylindrique.md §7 *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                    *
 ******************************************************************************/
#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>

namespace sofa::component::cosserat::shell {

using Scalar = double;

// ─── Material parameter structs ───────────────────────────────────────────────

/**
 * @brief Isotropic material parameters.
 *
 * Identical to ShellMaterialParams in ConstitutiveLaw.h — kept separate
 * to avoid include conflicts.
 */
struct IsotropicParams {
    Scalar E{1.0};          ///< Young's modulus [Pa]
    Scalar nu{0.3};         ///< Poisson's ratio [-]
    Scalar h{0.001};        ///< Thickness [m]
    Scalar kappa{5.0/6.0};  ///< Shear correction factor (Timoshenko) [-]
};

/**
 * @brief Orthotropic material parameters (engineering constants).
 *
 * The shell plane is (X_1, X_2); X_3 is the normal.
 * Symmetry: ν12/E1 = ν21/E2.
 */
struct OrthotropicParams {
    Scalar E1{1.0};         ///< Young's modulus in X_1 direction [Pa]
    Scalar E2{1.0};         ///< Young's modulus in X_2 direction [Pa]
    Scalar G12{0.385};      ///< In-plane shear modulus [Pa]
    Scalar G13{0.385};      ///< Transverse shear modulus (X_1–normal plane) [Pa]
    Scalar G23{0.385};      ///< Transverse shear modulus (X_2–normal plane) [Pa]
    Scalar nu12{0.3};       ///< Poisson ratio: contraction in X_2 when loading in X_1
    Scalar h{0.001};        ///< Thickness [m]
    Scalar kappa1{5.0/6.0}; ///< Shear correction factor for G13
    Scalar kappa2{5.0/6.0}; ///< Shear correction factor for G23

    /// Derived: ν21 = ν12 * E2 / E1
    Scalar nu21() const { return nu12 * E2 / E1; }
    /// Plane-stress denominator: 1 - ν12*ν21
    Scalar denom() const { return 1.0 - nu12 * nu21(); }
};

/**
 * @brief Cylindrical scaling parameters.
 *
 * When cylinderRadius > 0, the constitutive matrices are scaled according to
 * formulation_cylindrique.md §7:
 *   D^{θθ} = D^{11} / R²
 *   D^{zz} = D^{22}
 *   D^{θz} = D^{12}_Poisson / R
 *
 * Set cylinderRadius = 0 to use the flat (Cartesian) formulation.
 */
struct CylindricalScaling {
    Scalar cylinderRadius{0.0};  ///< Cylinder radius R [m], 0 = flat plate
    bool   isCylindrical() const { return cylinderRadius > 0.0; }
};

/**
 * @brief Pre-computed constitutive matrices for the full law.
 *
 * Same storage as ConstitutiveTensors:
 *   D[0]=D^{11}, D[1]=D^{12}, D[2]=D^{21}, D[3]=D^{22}
 */
struct ConstitutiveTensorsFull {
    std::array<Eigen::Matrix<Scalar, 6, 6>, 4> D;

    ConstitutiveTensorsFull() {
        for (auto& m : D) m.setZero();
    }

    /// Compute stress resultant S^α = Σ_β D^{αβ} · E_{tβ}
    Eigen::Matrix<Scalar,6,1> stressResultant(
        const Eigen::Matrix<Scalar,6,2>& E, int alpha) const
    {
        Eigen::Matrix<Scalar,6,1> S = Eigen::Matrix<Scalar,6,1>::Zero();
        for (int beta = 0; beta < 2; ++beta)
            S += D[2 * alpha + beta] * E.col(beta);
        return S;
    }

    /// Compute both stress resultants at once → ℝ^{6×2}
    Eigen::Matrix<Scalar,6,2> allStressResultants(
        const Eigen::Matrix<Scalar,6,2>& E) const
    {
        Eigen::Matrix<Scalar,6,2> S;
        S.col(0) = stressResultant(E, 0);
        S.col(1) = stressResultant(E, 1);
        return S;
    }
};

// ─── Internal helpers ─────────────────────────────────────────────────────────

namespace detail {

/**
 * @brief Build the 4 D^{αβ} matrices for a general orthotropic shell.
 *
 * Strain layout: [φ_x, φ_y, φ_z, ρ_x, ρ_y, ρ_z]
 *
 * Physical mapping (flat plate):
 *   Direction 1 (X_1): κ_11 → [0], twist1 → [1], twist2 → [2],
 *                       ε_11 → [3], ε_12   → [4], γ_1    → [5]
 *   Direction 2 (X_2): twist1→ [0], κ_22  → [1], twist2 → [2],
 *                       ε_21  → [3], ε_22   → [4], γ_2    → [5]
 *
 * Complete shear split: the symmetric shear energy
 *   C_m(ε_12 + ε_21)²  and  C_b(κ_12 + κ_21)²
 * is distributed equally between D^{11}, D^{22}, D^{12}, D^{21}.
 *
 * @param Am1   Membrane stiffness in X_1: E1*h/(1-ν12*ν21)
 * @param Am2   Membrane stiffness in X_2: E2*h/(1-ν12*ν21)
 * @param Bm    Poisson membrane coupling: ν12*E2*h/(1-ν12*ν21) = ν21*E1*h/(1-ν12*ν21)
 * @param Cm    In-plane shear membrane: G12*h
 * @param Db1   Bending stiffness in X_1: E1*h³/[12(1-ν12*ν21)]
 * @param Db2   Bending stiffness in X_2: E2*h³/[12(1-ν12*ν21)]
 * @param Bb    Poisson bending coupling: ν12*E2*h³/[12(1-ν12*ν21)]
 * @param Cb    Torsional stiffness:      G12*h³/12
 * @param ks1Gh Transverse shear X_1: κ1*G13*h
 * @param ks2Gh Transverse shear X_2: κ2*G23*h
 */
inline std::array<Eigen::Matrix<Scalar,6,6>, 4> buildDMatrices(
    Scalar Am1, Scalar Am2, Scalar Bm, Scalar Cm,
    Scalar Db1, Scalar Db2, Scalar Bb, Scalar Cb,
    Scalar ks1Gh, Scalar ks2Gh)
{
    std::array<Eigen::Matrix<Scalar,6,6>, 4> D;
    for (auto& m : D) m.setZero();

    // ── D^{11} ── strains from direction X_1 ────────────────────────────
    // Diagonal: full stiffness for the primary strains in direction 1
    // Shear energy split: half of C_m*ε_12² goes here, half to D^{22}
    D[0](0, 0) = Db1;          // κ_11² (primary bending)
    D[0](1, 1) = Cb / 2.0;    // twist from X_1 — half energy
    D[0](2, 2) = Cb / 2.0;    // second twist component — half energy
    D[0](3, 3) = Am1;          // ε_11² (primary membrane)
    D[0](4, 4) = Cm / 2.0;    // ε_12² — half energy
    D[0](5, 5) = ks1Gh;        // γ_1² (transverse shear X_1)

    // ── D^{22} ── strains from direction X_2 ────────────────────────────
    D[3](0, 0) = Cb / 2.0;    // twist from X_2 — half energy
    D[3](1, 1) = Db2;          // κ_22² (primary bending)
    D[3](2, 2) = Cb / 2.0;    // second twist component
    D[3](3, 3) = Cm / 2.0;    // ε_21² — half energy
    D[3](4, 4) = Am2;          // ε_22² (primary membrane)
    D[3](5, 5) = ks2Gh;        // γ_2² (transverse shear X_2)

    // ── D^{12} ── E_{t2} → S^1 ──────────────────────────────────────────
    // Poisson coupling + symmetric shear cross-terms
    D[1](0, 1) = Bb;           // M_11 += ν·Db·κ_22     (Poisson bending)
    D[1](1, 0) = Cb / 2.0;    // twist coupling (symmetric split)
    D[1](2, 2) = Cb / 2.0;    // second twist coupling
    D[1](3, 4) = Bm;           // N_11 += ν·Am·ε_22     (Poisson membrane)
    D[1](4, 3) = Cm / 2.0;    // shear: S^1[ε_12] from E_{t2}[ε_21]

    // ── D^{21} ── E_{t1} → S^2  (transpose of D^{12} for energy symmetry) ─
    D[2](1, 0) = Bb;           // M_22 += ν·Db·κ_11
    D[2](0, 1) = Cb / 2.0;    // twist coupling
    D[2](2, 2) = Cb / 2.0;
    D[2](4, 3) = Bm;           // N_22 += ν·Am·ε_11
    D[2](3, 4) = Cm / 2.0;    // shear: S^2[ε_21] from E_{t1}[ε_12]

    return D;
}

} // namespace detail

// ─── Public API ───────────────────────────────────────────────────────────────

/**
 * @brief Build the full constitutive tensors for an isotropic shell.
 *
 * Includes:
 *   - Physically correct D^{αβ} with symmetric shear energy split
 *   - Complete Poisson coupling in D^{12}, D^{21}
 *   - Transverse shear only on index [5] (ρ_z)
 *   - Optional cylindrical scaling
 *
 * @param iso    Isotropic material parameters
 * @param cyl    Cylindrical scaling (default: flat, R=0)
 * @return       ConstitutiveTensorsFull
 */
inline ConstitutiveTensorsFull computeConstitutiveTensorsFull(
    const IsotropicParams& iso,
    const CylindricalScaling& cyl = CylindricalScaling{})
{
    const Scalar E  = iso.E;
    const Scalar nu = iso.nu;
    const Scalar h  = iso.h;
    const Scalar ks = iso.kappa;

    const Scalar G   = E / (2.0 * (1.0 + nu));
    const Scalar den = 1.0 - nu * nu;

    // Membrane scalars
    const Scalar Am = E * h / den;
    const Scalar Bm = nu * Am;
    const Scalar Cm = G * h;           // in-plane shear

    // Bending scalars
    const Scalar Db = E * h * h * h / (12.0 * den);
    const Scalar Bb = nu * Db;
    const Scalar Cb = G * h * h * h / 12.0;  // torsional stiffness

    // Transverse shear
    const Scalar ksGh = ks * G * h;

    ConstitutiveTensorsFull result;
    result.D = detail::buildDMatrices(Am, Am, Bm, Cm, Db, Db, Bb, Cb, ksGh, ksGh);

    // ── Cylindrical scaling ────────────────────────────────────────────────
    // From formulation_cylindrique.md §7:
    //   D^{θθ} → D^{11} / R²   (X_1 = θ direction)
    //   D^{zz} → D^{22}        (X_2 = z direction, unchanged)
    //   D^{θz} → D^{12} / R    (cross-term scaled by 1/R)
    //   D^{zθ} → D^{21} / R
    if (cyl.isCylindrical()) {
        const Scalar R  = cyl.cylinderRadius;
        const Scalar R2 = R * R;
        result.D[0] /= R2;   // D^{11} = D^{θθ}
        result.D[1] /= R;    // D^{12} = D^{θz}
        result.D[2] /= R;    // D^{21} = D^{zθ}
        // D^{22} = D^{zz}: unchanged
    }

    return result;
}

/**
 * @brief Build the full constitutive tensors for an orthotropic shell.
 *
 * @param orth   Orthotropic material parameters (E1, E2, G12, G13, G23, ν12)
 * @param cyl    Cylindrical scaling (default: flat)
 * @return       ConstitutiveTensorsFull
 */
inline ConstitutiveTensorsFull computeConstitutiveTensorsOrthotropic(
    const OrthotropicParams& orth,
    const CylindricalScaling& cyl = CylindricalScaling{})
{
    const Scalar h   = orth.h;
    const Scalar den = orth.denom();   // 1 - ν12*ν21

    // Membrane
    const Scalar Am1 = orth.E1 * h / den;
    const Scalar Am2 = orth.E2 * h / den;
    const Scalar Bm  = orth.nu12 * orth.E2 * h / den;  // = ν21*E1*h/den
    const Scalar Cm  = orth.G12 * h;

    // Bending (factor h³/12)
    const Scalar h3_12 = h * h * h / 12.0;
    const Scalar Db1 = orth.E1 * h3_12 / den;
    const Scalar Db2 = orth.E2 * h3_12 / den;
    const Scalar Bb  = orth.nu12 * orth.E2 * h3_12 / den;
    const Scalar Cb  = orth.G12 * h3_12;

    // Transverse shear
    const Scalar ks1Gh = orth.kappa1 * orth.G13 * h;
    const Scalar ks2Gh = orth.kappa2 * orth.G23 * h;

    ConstitutiveTensorsFull result;
    result.D = detail::buildDMatrices(Am1, Am2, Bm, Cm, Db1, Db2, Bb, Cb, ks1Gh, ks2Gh);

    if (cyl.isCylindrical()) {
        const Scalar R  = cyl.cylinderRadius;
        const Scalar R2 = R * R;
        result.D[0] /= R2;
        result.D[1] /= R;
        result.D[2] /= R;
    }

    return result;
}

// ─── Convenience wrappers (same signature as ConstitutiveLaw.h) ───────────────

/**
 * @brief Compute stress resultant S^α from differential strain E.
 *
 * S^α = Σ_β D^{αβ} · E_{tβ}
 *
 * @param tensors  Full constitutive tensors
 * @param E        Differential strain ∈ ℝ^{6×2}
 * @param alpha    Direction index: 0 → X_1, 1 → X_2
 */
inline Eigen::Matrix<Scalar,6,1> computeStressResultantFull(
    const ConstitutiveTensorsFull& tensors,
    const Eigen::Matrix<Scalar,6,2>& E,
    int alpha)
{
    return tensors.stressResultant(E, alpha);
}

/**
 * @brief Compute both stress resultants S^1 and S^2.
 *
 * @return Matrix S ∈ ℝ^{6×2}, column α = S^{α+1}
 */
inline Eigen::Matrix<Scalar,6,2> computeAllStressResultantsFull(
    const ConstitutiveTensorsFull& tensors,
    const Eigen::Matrix<Scalar,6,2>& E)
{
    return tensors.allStressResultants(E);
}

} // namespace sofa::component::cosserat::shell
