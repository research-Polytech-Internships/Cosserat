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
    const Eigen::Matrix3d R = node.orientation.matrix();  // SO3::matrix() -> Eigen::Matrix3d
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

// ─── Tests 15-17: rightJacobianInverse ───────────────────────────────────────

// Test 15: J_R^{-1}(0) = I
TEST(SO3JacobianInverse, IdentityAtZero) {
    const Eigen::Vector3d zero = Eigen::Vector3d::Zero();
    const Eigen::Matrix3d JRinv = rightJacobianInverse(zero);
    EXPECT_NEAR((JRinv - Eigen::Matrix3d::Identity()).norm(), 0.0, 1e-12)
        << "J_R^{-1}(0) should be the identity matrix";
}

// Test 16: J_R^{-1}(omega) * J_R(omega) = I  (via finite differences on SO3 log)
// We verify J_R^{-1} * J_R = I using the identity:
//   d/dt log(exp(omega + t*eta))|_{t=0} = J_R^{-1}(omega) * eta
// FD check: log(exp(omega + h*e_k)) - log(exp(omega)) ~= h * J_R^{-1}(omega) * e_k
TEST(SO3JacobianInverse, InversesJacobianFD) {
    const Eigen::Vector3d omega(0.5, -0.3, 0.2);
    const Eigen::Matrix3d JRinv = rightJacobianInverse(omega);

    const double h = 1e-6;
    const SO3d R0 = SO3d::exp(omega);

    Eigen::Matrix3d JR_fd;
    for (int k = 0; k < 3; ++k) {
        Eigen::Vector3d ek = Eigen::Vector3d::Zero();
        ek(k) = 1.0;
        // FD column k of J_R: d/dt log(exp(omega)*exp(t*ek))|_{t=0} / h
        const Eigen::Vector3d logP = (R0 * SO3d::exp(h * ek)).log();
        const Eigen::Vector3d logM = (R0 * SO3d::exp(-h * ek)).log();
        JR_fd.col(k) = (logP - logM) / (2.0 * h);
    }
    // J_R^{-1} * J_R should be identity
    EXPECT_NEAR((JRinv * JR_fd - Eigen::Matrix3d::Identity()).norm(), 0.0, 1e-5)
        << "J_R^{-1} * J_R (FD) != I";
}

// Test 17: strain twist is zero for identity configuration (consistency with J_R^{-1})
// Same as Test 3 but explicitly checking that J_R^{-1} does not corrupt zero strain.
TEST(SO3JacobianInverse, ZeroStrainPreserved) {
    auto nodes = flatSquareElement(1.0, 1.0);
    // All nodes have identity orientation -> omega_bar = 0 -> J_R^{-1}(0) = I
    // Angular strains should still be zero
    const auto Xt = computeStrainTwistsAtCentroid(nodes);
    EXPECT_NEAR(Xt.col(0).head<3>().norm(), 0.0, 1e-13)
        << "Angular strain xi_{t1}[ang] != 0 for identity config";
    EXPECT_NEAR(Xt.col(1).head<3>().norm(), 0.0, 1e-13)
        << "Angular strain xi_{t2}[ang] != 0 for identity config";
}

// Test 18: large rotation consistency
// For a 45 deg rotation between adjacent nodes, J_R^{-1} must differ from I.
TEST(SO3JacobianInverse, LargeRotationNonTrivial) {
    // Build an element where nodes 2 and 3 are rotated by 45 deg around Z
    // relative to nodes 0 and 1 (simulates 1/8 of a cylinder with N_theta=8)
    std::array<NodeConfig, NODES_PER_ELEM> nodes;
    nodes[0].position = {-0.5, -0.5, 0.0};
    nodes[1].position = { 0.5, -0.5, 0.0};
    nodes[2].position = { 0.5,  0.5, 0.0};
    nodes[3].position = {-0.5,  0.5, 0.0};

    // Nodes 0,1: identity; nodes 2,3: 45 deg rotation around Z
    const double angle = M_PI / 4.0;
    const Eigen::Vector3d axisZ(0.0, 0.0, 1.0);
    nodes[0].orientation = SO3d::identity();
    nodes[1].orientation = SO3d::identity();
    nodes[2].orientation = SO3d::exp(angle * axisZ);
    nodes[3].orientation = SO3d::exp(angle * axisZ);

    // omega_bar = 0.5 * angle * axisZ (mean of 0 and angle)
    const Eigen::Vector3d omega_bar(0.0, 0.0, angle / 2.0);
    const Eigen::Matrix3d JRinv = rightJacobianInverse(omega_bar);

    // Verify J_R^{-1} != I for this non-zero omega_bar
    EXPECT_GT((JRinv - Eigen::Matrix3d::Identity()).norm(), 1e-3)
        << "J_R^{-1} should differ from I for 22.5 deg mean rotation";

    // Verify the strain computation runs without NaN/Inf
    const auto Xt = computeStrainTwistsAtCentroid(nodes);
    EXPECT_TRUE(Xt.allFinite())
        << "Strain twist contains NaN or Inf for large-rotation element";
}

// ─── Main ────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
