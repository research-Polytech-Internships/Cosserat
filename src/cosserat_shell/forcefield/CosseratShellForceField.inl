/******************************************************************************
 * Cosserat Shell — CosseratShellForceField implementation                   *
 ******************************************************************************/
#pragma once

#include <cosserat_shell/forcefield/CosseratShellForceField.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/linearalgebra/BaseMatrix.h>

namespace sofa::component::cosserat::shell {

// ── Constructor ─────────────────────────────────────────────────────────────

template<class DataTypes>
CosseratShellForceField<DataTypes>::CosseratShellForceField()
    : d_youngModulus(initData(&d_youngModulus, (SReal)1000.0, "youngModulus",
        "Young's modulus E [Pa].\n"
        "Typical soft silicone: 0.1-1 MPa; rubber: 1-10 MPa."))
    , d_poissonRatio(initData(&d_poissonRatio, (SReal)0.3, "poissonRatio",
        "Poisson's ratio ν ∈ [0, 0.5).\n"
        "For nearly incompressible: ν ≈ 0.49."))
    , d_thickness(initData(&d_thickness, (SReal)0.001, "thickness",
        "Shell thickness h [m]."))
    , d_shearCorrectionFactor(initData(&d_shearCorrectionFactor, (SReal)(5.0/6.0),
        "shearCorrectionFactor",
        "Shear correction factor κ (default 5/6, Reissner-Mindlin)."))
    , d_quads(initData(&d_quads, "quads",
        "Quad connectivity: each entry = [n0, n1, n2, n3] (counter-clockwise)."))
    , d_restPosition(initData(&d_restPosition, "restPosition",
        "Reference (rest) configuration. If empty, uses initial MechanicalState position."))
{
}

// ── init ─────────────────────────────────────────────────────────────────────

template<class DataTypes>
void CosseratShellForceField<DataTypes>::init()
{
    sofa::core::behavior::ForceField<DataTypes>::init();

    // ── Build constitutive tensors ─────────────────────────────────────
    ShellMaterialParams params;
    params.E     = static_cast<Scalar>(d_youngModulus.getValue());
    params.nu    = static_cast<Scalar>(d_poissonRatio.getValue());
    params.h     = static_cast<Scalar>(d_thickness.getValue());
    params.kappa = static_cast<Scalar>(d_shearCorrectionFactor.getValue());
    m_tensors = computeConstitutiveTensors(params);

    // ── Determine rest position ────────────────────────────────────────
    const VecCoord* restPos = nullptr;
    if (d_restPosition.getValue().empty()) {
        restPos = &(this->mstate->read(sofa::core::VecCoordId::position())->getValue());
    } else {
        restPos = &d_restPosition.getValue();
    }

    // ── Build reference elements ───────────────────────────────────────
    buildReferenceElements(*restPos);

    msg_info() << "CosseratShellForceField initialized: "
               << m_elements.size() << " Q4 elements, "
               << "E=" << params.E << " Pa, nu=" << params.nu
               << ", h=" << params.h << " m.";
}

// ── reinit ────────────────────────────────────────────────────────────────────

template<class DataTypes>
void CosseratShellForceField<DataTypes>::reinit()
{
    init();
}

// ── buildReferenceElements ────────────────────────────────────────────────────

template<class DataTypes>
void CosseratShellForceField<DataTypes>::buildReferenceElements(const VecCoord& restPos)
{
    const auto& quads = d_quads.getValue();
    const std::size_t nElem = quads.size();

    m_elements.resize(nElem);
    m_elemStates.resize(nElem);

    for (std::size_t e = 0; e < nElem; ++e) {
        m_elements[e].nodeIds = {
            (int)quads[e][0], (int)quads[e][1],
            (int)quads[e][2], (int)quads[e][3]
        };

        // Build reference node configs
        std::array<NodeConfig, NODES_PER_ELEM> refNodes;
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            refNodes[i] = coordToNodeConfig(restPos[quads[e][i]]);
        }

        initShellElementRef(refNodes, m_elements[e]);
    }
}

// ── coordToNodeConfig ─────────────────────────────────────────────────────────

template<class DataTypes>
NodeConfig CosseratShellForceField<DataTypes>::coordToNodeConfig(const Coord& c) const
{
    NodeConfig node;
    // Convention: Coord = [ω_x, ω_y, ω_z, v_x, v_y, v_z]
    // For rest position: ω usually = [0,0,0] → orientation = identity
    const Eigen::Vector3d omega(c[0], c[1], c[2]);
    node.orientation = SO3d::exp(omega);
    node.position    = Eigen::Vector3d(c[3], c[4], c[5]);
    return node;
}

// ── getElementNodes ───────────────────────────────────────────────────────────

template<class DataTypes>
std::array<NodeConfig, NODES_PER_ELEM>
CosseratShellForceField<DataTypes>::getElementNodes(
    const VecCoord& x, std::size_t elemIdx) const
{
    std::array<NodeConfig, NODES_PER_ELEM> nodes;
    const auto& ids = m_elements[elemIdx].nodeIds;
    for (int i = 0; i < NODES_PER_ELEM; ++i) {
        nodes[i] = coordToNodeConfig(x[ids[i]]);
    }
    return nodes;
}

// ── addForce ──────────────────────────────────────────────────────────────────

template<class DataTypes>
void CosseratShellForceField<DataTypes>::addForce(
    const sofa::core::MechanicalParams* /*mparams*/,
    DataVecDeriv& f,
    const DataVecCoord& x,
    const DataVecDeriv& /*v*/)
{
    sofa::helper::AdvancedTimer::stepBegin("CosseratShellForceField::addForce");

    auto& fv      = sofa::helper::getWriteAccessor(f);
    const auto& xv = x.getValue();

    const std::size_t nElem = m_elements.size();

    for (std::size_t e = 0; e < nElem; ++e) {
        // ── Get current node configurations ──────────────────────────────
        const auto nodes = getElementNodes(xv, e);

        // ── Compute element force and stiffness ───────────────────────────
        computeElementForceAndStiffness(nodes, m_elements[e], m_tensors, m_elemStates[e]);

        // ── Scatter element force to global vector ─────────────────────────
        const auto& ids = m_elements[e].nodeIds;
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            const int globalIdx = ids[i];
            // Element force block for node i: Fint[6i : 6i+6]
            for (int d = 0; d < DOF_PER_NODE; ++d) {
                fv[globalIdx][d] -= static_cast<typename Deriv::value_type>(
                    m_elemStates[e].Fint(6*i + d));
            }
        }
    }

    sofa::helper::AdvancedTimer::stepEnd("CosseratShellForceField::addForce");
}

// ── addDForce ─────────────────────────────────────────────────────────────────

template<class DataTypes>
void CosseratShellForceField<DataTypes>::addDForce(
    const sofa::core::MechanicalParams* mparams,
    DataVecDeriv& df,
    const DataVecDeriv& dx)
{
    const SReal kFactor = sofa::core::mechanicalparams::kFactorIncludingRayleighDamping(mparams);

    auto& dfv      = sofa::helper::getWriteAccessor(df);
    const auto& dxv = dx.getValue();

    const std::size_t nElem = m_elements.size();

    for (std::size_t e = 0; e < nElem; ++e) {
        const auto& ids = m_elements[e].nodeIds;
        const auto& Ke  = m_elemStates[e].K;

        // Gather dx for this element: 24-vector
        Eigen::Matrix<double, DOF_PER_ELEM, 1> dx_e;
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            for (int d = 0; d < DOF_PER_NODE; ++d) {
                dx_e(6*i + d) = static_cast<double>(dxv[ids[i]][d]);
            }
        }

        // df_e = -kFactor * Ke * dx_e
        const Eigen::Matrix<double, DOF_PER_ELEM, 1> df_e = -kFactor * Ke * dx_e;

        // Scatter
        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            for (int d = 0; d < DOF_PER_NODE; ++d) {
                dfv[ids[i]][d] += static_cast<typename Deriv::value_type>(df_e(6*i + d));
            }
        }
    }
}

// ── addKToMatrix ─────────────────────────────────────────────────────────────

template<class DataTypes>
void CosseratShellForceField<DataTypes>::addKToMatrix(
    sofa::linearalgebra::BaseMatrix* matrix,
    SReal kFactor,
    unsigned int& offset)
{
    const std::size_t nElem = m_elements.size();

    for (std::size_t e = 0; e < nElem; ++e) {
        const auto& ids = m_elements[e].nodeIds;
        const auto& Ke  = m_elemStates[e].K;

        for (int i = 0; i < NODES_PER_ELEM; ++i) {
            for (int j = 0; j < NODES_PER_ELEM; ++j) {
                for (int di = 0; di < DOF_PER_NODE; ++di) {
                    for (int dj = 0; dj < DOF_PER_NODE; ++dj) {
                        const int row = offset + ids[i] * DOF_PER_NODE + di;
                        const int col = offset + ids[j] * DOF_PER_NODE + dj;
                        matrix->add(row, col,
                            -kFactor * Ke(6*i+di, 6*j+dj));
                    }
                }
            }
        }
    }
}

// ── getPotentialEnergy ────────────────────────────────────────────────────────

template<class DataTypes>
SReal CosseratShellForceField<DataTypes>::getPotentialEnergy(
    const sofa::core::MechanicalParams* /*mparams*/,
    const DataVecCoord& x) const
{
    const auto& xv = x.getValue();
    SReal W = 0.0;

    for (std::size_t e = 0; e < m_elements.size(); ++e) {
        const auto nodes = getElementNodes(xv, e);
        ShellElementCurrent state;
        computeElementForceAndStiffness(nodes, m_elements[e], m_tensors, state);

        // W += ½ ∫ ⟨E, S⟩ j_0 dA  (approximated with Gauss integration)
        // For constant strains: W_e ≈ ½ j_0 * A_e * ⟨E, S⟩
        //   where A_e = Σ_g j_0^g * w_g ≈ 4 * j0_centroid (for uniform mesh)
        const double A_e = 4.0 * m_elements[e].j0;  // approximate element area
        for (int alpha = 0; alpha < 2; ++alpha) {
            W += 0.5 * A_e * state.E.col(alpha).dot(state.S.col(alpha));
        }
    }

    return W;
}

} // namespace sofa::component::cosserat::shell
