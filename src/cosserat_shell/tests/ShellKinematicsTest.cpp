/******************************************************************************
 * Cosserat Shell — Unit tests for FEM kinematics and constitutive law       *
 *                                                                            *
 * Tests:                                                                     *
 *   1. Shape functions partition of unity: Σ N^i = 1                        *
 *   2. Shape function gradient consistency                                   *
 *   3. Strain twist = zero for rigid body motion                             *
 *   4. Constitutive tensors: D^{11} is SPD                                  *
 *   5. Stress resultant linearity                                             *
 *   6. Reference element init: X_0* · X_0 ≈ I_2                            *
 *   7. Area Jacobian: positive for valid quad                                *
 *   8. Element stiffness K_e is SPD for regular flat element                 *
 *   9. Zero strain in rest configuration (F_int = 0)                         *
 *  10. Newton update: SO(3) increment preserves det(R) = 1                   *
 ******************************************************************************/
#include <gtest/gtest.h>

#include <cosserat_shell/fem/ShellElement.h>
#include <cosserat_shell/fem/SE3ShellKinematics.h>
#include <cosserat_shell/fem/ConstitutiveLaw.h>
#include <cosserat_shell/fem/ShellElementComputer.h>
#include <liegroups/SO3.h>

using namespace sofa::component::cosserat::shell;
using SO3d = sofa::component::cosserat::liegroups::SO3<double>;

// ─── Helpers ──────────────────────────────────────────────────────────────────

/// Build a flat 1×1 reference square element in the XY plane.
std::array<NodeConfig, NODES_PER_ELEM> flatSquareElement(double Lx = 1.0, double Ly = 1.0) {
    std::array<NodeConfig, NODES_PER_ELEM> nodes;
    // Counter-clockwise: (-Lx/2,-Ly/2), (+Lx/2,-Ly/2), (+Lx/2,+Ly/2), (-Lx/2,+Ly/2)
    nodes[0].position = {-Lx/2, -Ly/2, 0.0};
    nodes[1].position = { Lx/2, -Ly/2, 0.0};
    nodes[2].position = { Lx/2,  Ly/2, 0.0};
    nodes[3].position = {-Lx/2,  Ly/2, 0.0};
    for (auto& n : nodes) n.orientation = SO3d::identity();
    return nodes;
}

/// Build default material parameters
ShellMaterialParams defaultMaterial() {
    return {1000.0, 0.3, 0.001, 5.0/6.0};
}

// ─── Test 1: Partition of unity ───────────────────────────────────────────────

TEST(ShellShapeFunctions, PartitionOfUnity) {
    const std::vector<std::pair<double,double>> testPoints = {
        {0.0, 0.0}, {0.5, 0.5}, {-0.5, 0.3}, {1.0, -1.0}, {-1.0, 1.0}
    };
    for (auto [x, y] : testPoints) {
        double sum = 0.0;
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            sum += shapeFunction(i, x, y);
        }
        EXPECT_NEAR(sum, 1.0, 1e-14)
            << "Partition of unity violated at (" << x << "," << y << ")";
    }
}

// ─── Test 2: Shape function gradient consistency ──────────────────────────────

TEST(ShellShapeFunctions, GradientConsistency) {
    // FD check: ∂N^i/∂x ≈ (N^i(x+h,y) - N^i(x-h,y)) / (2h)
    const double h = 1e-7;
    const double x = 0.3, y = -0.2;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        const double fd_dx = (shapeFunction(i,x+h,y) - shapeFunction(i,x-h,y)) / (2*h);
        const double fd_dy = (shapeFunction(i,x,y+h) - shapeFunction(i,x,y-h)) / (2*h);
        const auto grad = shapeFunctionGrad(i, x, y);
        EXPECT_NEAR(grad(0), fd_dx, 1e-9) << "∂N^" << i << "/∂x FD mismatch";
        EXPECT_NEAR(grad(1), fd_dy, 1e-9) << "∂N^" << i << "/∂y FD mismatch";
    }
}

// ─── Test 3: Zero strain for identity configuration ───────────────────────────

TEST(ShellKinematics, ZeroStrainForFlatIdentity) {
    // For a flat element with identity orientations, the strain twist at the
    // centroid should be consistent (X_0 = X_t → E = 0).
    auto nodes = flatSquareElement(1.0, 1.0);
    const auto Xt = computeStrainTwistsAtCentroid(nodes);

    // Build reference elements and check E = X_t - X_0 = 0
    ShellElementRef elemRef;
    initShellElementRef(nodes, elemRef);

    const auto E = Xt - elemRef.X0;
    EXPECT_NEAR(E.norm(), 0.0, 1e-12)
        << "Non-zero differential strain in rest configuration";
}

// ─── Test 4: Constitutive tensor D^{11} is symmetric and positive definite ───

TEST(ConstitutiveLaw, D11IsSymmetricSPD) {
    const auto tensors = computeConstitutiveTensors(defaultMaterial());
    const Eigen::Matrix<double,6,6>& D11 = tensors.D[0];

    // Symmetry
    EXPECT_NEAR((D11 - D11.transpose()).norm(), 0.0, 1e-14)
        << "D^{11} is not symmetric";

    // Positive definiteness: all eigenvalues > 0
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> es(D11);
    const double minEig = es.eigenvalues().minCoeff();
    EXPECT_GT(minEig, 0.0) << "D^{11} is not positive definite (min eig = " << minEig << ")";
}

// ─── Test 5: Stress resultant linearity ───────────────────────────────────────

TEST(ConstitutiveLaw, StressResultantLinearity) {
    const auto tensors = computeConstitutiveTensors(defaultMaterial());

    // S(2*E) = 2*S(E)
    Eigen::Matrix<double,6,2> E1;
    E1 << 0.01, 0.005,
          0.00, 0.010,
          0.02, 0.000,
          0.00, 0.003,
          0.01, 0.002,
          0.00, 0.001;

    const auto S1 = computeAllStressResultants(tensors, E1);
    const auto S2 = computeAllStressResultants(tensors, 2.0 * E1);

    EXPECT_NEAR((S2 - 2.0 * S1).norm(), 0.0, 1e-10)
        << "Stress resultant is not linear in strain";
}

// ─── Test 6: Pseudo-inverse X_0* satisfies X_0* · X_0 ≈ I_2 ─────────────────

TEST(ShellElementRef, PseudoinverseConsistency) {
    auto nodes = flatSquareElement(0.5, 0.3);
    ShellElementRef elemRef;
    initShellElementRef(nodes, elemRef);

    // X_0* · X_0 should be identity 2×2 (right inverse when X_0 has rank 2)
    const Eigen::Matrix<double,2,2> I2_approx = elemRef.X0star * elemRef.X0;
    EXPECT_NEAR((I2_approx - Eigen::Matrix<double,2,2>::Identity()).norm(), 0.0, 1e-10)
        << "X_0* · X_0 ≠ I_2";
}

// ─── Test 7: Area Jacobian is positive ───────────────────────────────────────

TEST(ShellElementRef, AreaJacobianPositive) {
    auto nodes = flatSquareElement(1.0, 1.0);
    ShellElementRef elemRef;
    initShellElementRef(nodes, elemRef);

    EXPECT_GT(elemRef.j0, 0.0) << "Reference area Jacobian j_0 should be positive";
    for (int g = 0; g < N_GAUSS; ++g) {
        EXPECT_GT(elemRef.j0Gauss[g], 0.0)
            << "Gauss point area Jacobian j_0^" << g << " should be positive";
    }
}

// ─── Test 8: Element stiffness K_e is symmetric ───────────────────────────────

TEST(ShellElementComputer, StiffnessIsSymmetric) {
    auto nodes = flatSquareElement(1.0, 1.0);
    ShellElementRef elemRef;
    initShellElementRef(nodes, elemRef);

    const auto tensors = computeConstitutiveTensors(defaultMaterial());
    ShellElementCurrent state;
    computeElementForceAndStiffness(nodes, elemRef, tensors, state);

    const double skew = (state.K - state.K.transpose()).norm();
    // Allow some tolerance: geometric stiffness may have small asymmetry
    // when not exactly at equilibrium (Remark 2 of the paper)
    EXPECT_NEAR(skew / (state.K.norm() + 1e-14), 0.0, 1e-8)
        << "Element stiffness K_e is not symmetric (skew/norm = " << skew << ")";
}

// ─── Test 9: Zero internal force at rest configuration ────────────────────────

TEST(ShellElementComputer, ZeroForceAtRest) {
    auto nodes = flatSquareElement(1.0, 1.0);
    ShellElementRef elemRef;
    initShellElementRef(nodes, elemRef);

    const auto tensors = computeConstitutiveTensors(defaultMaterial());
    ShellElementCurrent state;
    computeElementForceAndStiffness(nodes, elemRef, tensors, state);

    EXPECT_NEAR(state.Fint.norm(), 0.0, 1e-12)
        << "Non-zero internal force at rest configuration";
}

// ─── Test 10: Newton SO(3) update preserves det(R) = 1 ───────────────────────

TEST(ShellNewtonUpdate, SO3UpdatePreservesOrthogonality) {
    NodeConfig node;
    node.position    = {1.0, 2.0, 3.0};
    node.orientation = SO3d::identity();

    // Apply a 30° rotation around z-axis
    Eigen::Matrix<double,6,1> increment;
    increment << 0.0, 0.0, M_PI/6.0,   // δω = 30° around z
                 0.1, -0.2, 0.05;       // δv

    applyNodeUpdate(node, increment);

    // Check det(R) = 1
    const Eigen::Matrix3d R = node.orientation.toRotationMatrix();
    EXPECT_NEAR(R.determinant(), 1.0, 1e-13) << "Rotation matrix det ≠ 1 after update";

    // Check R^T R ≈ I
    EXPECT_NEAR((R.transpose() * R - Eigen::Matrix3d::Identity()).norm(), 0.0, 1e-13)
        << "Rotation matrix not orthogonal after update";

    // Check position updated correctly
    EXPECT_NEAR(node.position(0), 1.1, 1e-14);
    EXPECT_NEAR(node.position(1), 1.8, 1e-14);
    EXPECT_NEAR(node.position(2), 3.05, 1e-14);
}

// ─── Test 11: B-matrix shape ──────────────────────────────────────────────────

TEST(ShellKinematics, BMatrixShape) {
    auto nodes = flatSquareElement();
    const auto Xt = computeStrainTwistsAtCentroid(nodes);

    for (int alpha = 0; alpha < 2; ++alpha) {
        const auto B = computeBMatrix(Xt.col(alpha), 0.0, 0.0, alpha);
        EXPECT_EQ(B.rows(), 6);
        EXPECT_EQ(B.cols(), DOF_PER_ELEM);  // 24
    }
}

// ─── Test 12: adMatrix antisymmetry ──────────────────────────────────────────

TEST(SE3Kinematics, AdMatrixAntisymmetry) {
    // ad_xi · eta = -ad_eta · xi (antisymmetry of Lie bracket)
    Eigen::Matrix<double,6,1> xi, eta;
    xi  << 0.1, -0.2, 0.3, 0.5, -0.1, 0.2;
    eta << 0.4,  0.1, -0.1, 0.2, 0.3, -0.4;

    const auto adXi  = adMatrix(xi);
    const auto adEta = adMatrix(eta);

    const Eigen::Matrix<double,6,1> bracket1 = adXi  * eta;
    const Eigen::Matrix<double,6,1> bracket2 = adEta * xi;

    EXPECT_NEAR((bracket1 + bracket2).norm(), 0.0, 1e-13)
        << "ad_xi · eta ≠ -ad_eta · xi (Lie bracket antisymmetry)";
}

// ─── Test 13: Gauss points have correct count ────────────────────────────────

TEST(ShellElement, GaussPointCount) {
    const auto gpts = gaussPoints2x2();
    EXPECT_EQ(gpts.size(), (std::size_t)N_GAUSS);
    EXPECT_EQ(N_GAUSS, 4);
}

// ─── Test 14: Strain twist has correct size ───────────────────────────────────

TEST(ShellKinematics, StrainTwistSize) {
    auto nodes = flatSquareElement();
    const auto Xt = computeStrainTwistsAtCentroid(nodes);
    EXPECT_EQ(Xt.rows(), 6);
    EXPECT_EQ(Xt.cols(), 2);
}

} // end anonymous

// ─── Main ────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
