/******************************************************************************
 * Cosserat Shell — FEM Element Data Structures                               *
 *                                                                            *
 * Geometrically-exact Cosserat shell element (Q4 isoparametric, 4 nodes).   *
 * Reference: cosserat_shell_v0.pdf                                          *
 *                                                                            *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                   *
 ******************************************************************************/
#pragma once

#include <array>
#include <Eigen/Dense>

namespace sofa::component::cosserat::shell {

// ─── Scalar type ───────────────────────────────────────────────────────────
using Scalar   = double;

// ─── Compile-time constants ────────────────────────────────────────────────
static constexpr int NODES_PER_ELEM = 4;   ///< Q4 element
static constexpr int DOF_PER_NODE   = 6;   ///< 3 position + 3 rotation (se(3))
static constexpr int DOF_PER_ELEM   = NODES_PER_ELEM * DOF_PER_NODE; ///< 24

// ─── Local node positions in the reference element [-1,1]^2 ────────────────
/// Node coordinates in reference element: (-1,-1),(+1,-1),(+1,+1),(-1,+1)
static constexpr double LOCAL_COORDS[4][2] = {
    {-1.0, -1.0},
    {+1.0, -1.0},
    {+1.0, +1.0},
    {-1.0, +1.0}
};

// ─── Gauss quadrature (2x2 full integration) ──────────────────────────────
static constexpr int    N_GAUSS   = 4;
static constexpr double GAUSS_PT  = 1.0 / 1.7320508075688772;  // 1/√3
static constexpr double GAUSS_W   = 1.0;

/// 2×2 Gauss points (xi, eta) and weights
struct GaussPoint {
    double xi, eta, w;
};

inline std::array<GaussPoint, N_GAUSS> gaussPoints2x2() {
    return {{
        {-GAUSS_PT, -GAUSS_PT, GAUSS_W},
        {+GAUSS_PT, -GAUSS_PT, GAUSS_W},
        {+GAUSS_PT, +GAUSS_PT, GAUSS_W},
        {-GAUSS_PT, +GAUSS_PT, GAUSS_W}
    }};
}

// ─── Bilinear shape functions ──────────────────────────────────────────────

/**
 * @brief Bilinear shape function N^i(x,y) = ¼(1+x·xi)(1+y·yi)
 *
 * @param i   Node index in [0,3]
 * @param x   Reference coordinate ξ¹ ∈ [-1,1]
 * @param y   Reference coordinate ξ² ∈ [-1,1]
 * @return    Value of N^i at (x,y)
 */
inline double shapeFunction(int i, double x, double y) {
    return 0.25 * (1.0 + LOCAL_COORDS[i][0] * x)
                * (1.0 + LOCAL_COORDS[i][1] * y);
}

/**
 * @brief Gradient of shape function ∂N^i/∂ξ^α at (x,y).
 *
 * @param i     Node index in [0,3]
 * @param x     Reference coordinate ξ¹ ∈ [-1,1]
 * @param y     Reference coordinate ξ² ∈ [-1,1]
 * @return      [∂N^i/∂x, ∂N^i/∂y]
 */
inline Eigen::Vector2d shapeFunctionGrad(int i, double x, double y) {
    const double xi_i = LOCAL_COORDS[i][0];
    const double eta_i = LOCAL_COORDS[i][1];
    return {
        0.25 * xi_i  * (1.0 + eta_i * y),
        0.25 * eta_i * (1.0 + xi_i  * x)
    };
}

// ─── Pre-computed element data (reference configuration) ──────────────────

/**
 * @brief Pre-computed data for one shell element in the REFERENCE configuration.
 *
 * These quantities depend only on g_0 and are computed once in init().
 */
struct ShellElementRef {
    /// Global node indices (4 nodes, in counter-clockwise order)
    std::array<int, NODES_PER_ELEM> nodeIds{};

    /// Reference strain twists at CENTROID: X_0 = [ζ_{01} | ζ_{02}] ∈ ℝ^{6×2}
    Eigen::Matrix<Scalar, 6, 2> X0;

    /// Pseudo-inverse of X_0: X_0* = X_0^T (X_0 X_0^T)^{-1} ∈ ℝ^{2×6}
    Eigen::Matrix<Scalar, 2, 6> X0star;

    /// Area Jacobian at centroid: j_0 = ||∂φ_0/∂ξ¹ × ∂φ_0/∂ξ²|| ∈ ℝ
    Scalar j0{1.0};

    /// Gauss point Jacobians j_0^g (for integration — may differ from centroid)
    std::array<Scalar, N_GAUSS> j0Gauss{};

    ShellElementRef() : X0(Eigen::Matrix<Scalar,6,2>::Zero()),
                        X0star(Eigen::Matrix<Scalar,2,6>::Zero()) {}
};

// ─── Current element state ─────────────────────────────────────────────────

/**
 * @brief Current-configuration quantities computed each iteration.
 */
struct ShellElementCurrent {
    /// Current strain twists at CENTROID: X_t = [ζ_{t1} | ζ_{t2}] ∈ ℝ^{6×2}
    Eigen::Matrix<Scalar, 6, 2> Xt;

    /// Differential strain: E = X_t - X_0 ∈ ℝ^{6×2}
    Eigen::Matrix<Scalar, 6, 2> E;

    /// Stress resultants: S^α = Σ_β D^{αβ} E_β ∈ ℝ⁶, stored as 6×2 matrix
    Eigen::Matrix<Scalar, 6, 2> S;

    /// Element residual force vector (24 DOF)
    Eigen::Matrix<Scalar, DOF_PER_ELEM, 1> Fint;

    /// Element tangent stiffness matrix (24×24)
    Eigen::Matrix<Scalar, DOF_PER_ELEM, DOF_PER_ELEM> K;

    ShellElementCurrent()
        : Xt(Eigen::Matrix<Scalar,6,2>::Zero()),
          E(Eigen::Matrix<Scalar,6,2>::Zero()),
          S(Eigen::Matrix<Scalar,6,2>::Zero()),
          Fint(Eigen::Matrix<Scalar,DOF_PER_ELEM,1>::Zero()),
          K(Eigen::Matrix<Scalar,DOF_PER_ELEM,DOF_PER_ELEM>::Zero()) {}
};

} // namespace sofa::component::cosserat::shell
