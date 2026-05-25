/******************************************************************************
 * Cosserat Shell — SOFA ForceField Component                                 *
 *                                                                            *
 * Geometrically-exact Cosserat shell forcefield based on SE(3) Lie group    *
 * theory. Implements the internal virtual work principle (Eq. 32 of the     *
 * paper) with a 4-node isoparametric Q4 element and centroid strain          *
 * evaluation (shear-locking free).                                           *
 *                                                                            *
 * DOF layout per node: [ω_x, ω_y, ω_z, v_x, v_y, v_z]                     *
 * (angular velocity / curvature first, then translational)                   *
 *                                                                            *
 * Reference: cosserat_shell_v0.pdf                                          *
 * Author: Y. Adagolodjo (DEFROST / INRIA)                                   *
 ******************************************************************************/
#pragma once

// ── SOFA includes ──────────────────────────────────────────────────────────
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/vector.h>

// ── Shell includes ─────────────────────────────────────────────────────────
#include <cosserat_shell/fem/ShellElement.h>
#include <cosserat_shell/fem/ShellElementComputer.h>
#include <cosserat_shell/fem/ConstitutiveLaw.h>

// ── Standard ───────────────────────────────────────────────────────────────
#include <vector>

namespace sofa::component::cosserat::shell {

/**
 * @brief Geometrically-exact Cosserat shell ForceField.
 *
 * Template parameter DataTypes must provide:
 *   - Coord: 6D vector [ω₁, ω₂, ω₃, v₁, v₂, v₃] per node
 *   - Deriv: same, tangent vector
 *   - VecCoord, VecDeriv: containers
 *
 * Usage in SOFA scene:
 * ```xml
 * <CosseratShellForceField
 *     youngModulus="1000"
 *     poissonRatio="0.3"
 *     thickness="0.001"
 *     topology="@topo"
 * />
 * ```
 */
template<class DataTypes>
class CosseratShellForceField : public sofa::core::behavior::ForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CosseratShellForceField, DataTypes),
               SOFA_TEMPLATE(sofa::core::behavior::ForceField, DataTypes));

    // ── Type aliases ──────────────────────────────────────────────────────
    using Coord    = typename DataTypes::Coord;
    using Deriv    = typename DataTypes::Deriv;
    using VecCoord = typename DataTypes::VecCoord;
    using VecDeriv = typename DataTypes::VecDeriv;
    using DataVecCoord = sofa::Data<VecCoord>;
    using DataVecDeriv = sofa::Data<VecDeriv>;

    using Index = sofa::Index;

    // ── Constructor / Destructor ─────────────────────────────────────────
    CosseratShellForceField();
    ~CosseratShellForceField() override = default;

    // ── SOFA lifecycle ───────────────────────────────────────────────────

    /**
     * @brief Initialize the forcefield.
     *
     * Reads the topology (quads connectivity), initializes reference
     * element data (X_0, X_0*, j_0), and builds the constitutive tensors.
     */
    void init() override;

    /**
     * @brief Recompute stiffness (called when topology changes).
     */
    void reinit() override;

    // ── ForceField API ────────────────────────────────────────────────────

    /**
     * @brief Compute internal forces (residual of the weak form).
     *
     * Assembles F_int = -G_int evaluated at the current configuration.
     * Each node contribution: F^i = Σ_e Σ_α B_α^{i,T} · S^α · j_0 · w
     */
    void addForce(const sofa::core::MechanicalParams* mparams,
                  DataVecDeriv& f,
                  const DataVecCoord& x,
                  const DataVecDeriv& v) override;

    /**
     * @brief Compute linearized force increment (for implicit solvers).
     *
     * df += kFactor * K · dx
     * where K is the tangent stiffness assembled in addForce().
     */
    void addDForce(const sofa::core::MechanicalParams* mparams,
                   DataVecDeriv& df,
                   const DataVecDeriv& dx) override;

    /**
     * @brief Assemble contribution to the global stiffness matrix.
     */
    void addKToMatrix(sofa::linearalgebra::BaseMatrix* matrix,
                      SReal kFactor,
                      unsigned int& offset) override;

    /**
     * @brief Return potential energy W_int = ½ Σ_e ∫ ⟨E, D E⟩ j_0 dA.
     */
    SReal getPotentialEnergy(const sofa::core::MechanicalParams* mparams,
                             const DataVecCoord& x) const override;

    // ── Data fields ──────────────────────────────────────────────────────

    /// Young's modulus E [Pa]
    sofa::Data<SReal> d_youngModulus;

    /// Poisson's ratio ν ∈ [0, 0.5)
    sofa::Data<SReal> d_poissonRatio;

    /// Shell thickness h [m]
    sofa::Data<SReal> d_thickness;

    /// Shear correction factor κ (default 5/6 for Reissner-Mindlin)
    sofa::Data<SReal> d_shearCorrectionFactor;

    /// Quad connectivity: each row = [n0, n1, n2, n3] global node indices
    sofa::Data<sofa::type::vector<sofa::type::fixed_array<Index,4>>> d_quads;

    /// Rest position (reference configuration).
    /// If empty, uses the initial position from the MechanicalState.
    sofa::Data<VecCoord> d_restPosition;

    // ── Accessors (for testing / visualization) ───────────────────────────

    /// Number of elements
    [[nodiscard]] std::size_t getNumElements() const { return m_elements.size(); }

    /// Reference element data (read-only)
    [[nodiscard]] const ShellElementRef& getElementRef(std::size_t e) const {
        return m_elements[e];
    }

    /// Current element state (valid after last addForce() call)
    [[nodiscard]] const ShellElementCurrent& getElementCurrent(std::size_t e) const {
        return m_elemStates[e];
    }

private:
    // ── Internal helpers ──────────────────────────────────────────────────

    /**
     * @brief Convert a SOFA Coord (6D vector) to a NodeConfig.
     *
     * Convention: coord = [ω_x, ω_y, ω_z, v_x, v_y, v_z]
     * → position = [v_x, v_y, v_z], orientation = exp([ω_x, ω_y, ω_z])
     *
     * Note: This convention is for the *reference* configuration init only.
     * During simulation, positions and orientations are tracked separately.
     */
    NodeConfig coordToNodeConfig(const Coord& c) const;

    /**
     * @brief Build element node configurations from SOFA position vector.
     */
    std::array<NodeConfig, NODES_PER_ELEM> getElementNodes(
        const VecCoord& x, std::size_t elemIdx) const;

    /**
     * @brief Initialize/rebuild all reference element data.
     * Called once in init() or reinit().
     */
    void buildReferenceElements(const VecCoord& restPos);

    // ── Internal state ────────────────────────────────────────────────────

    /// Reference data per element (computed in init/reinit)
    std::vector<ShellElementRef> m_elements;

    /// Current element state (updated each addForce() call)
    mutable std::vector<ShellElementCurrent> m_elemStates;

    /// Constitutive tensors (rebuilt when material params change)
    ConstitutiveTensors m_tensors;

    /// Cached kFactor for addDForce
    SReal m_kFactor{1.0};
};

// ── Explicit instantiation declaration ────────────────────────────────────
// (defined in CosseratShellForceField.cpp)
#if !defined(SOFA_COSSERAT_SHELL_CPP_CosseratShellForceField)
extern template class SOFA_COSSERAT_API CosseratShellForceField<sofa::defaulttype::Vec6Types>;
#endif

// ── Factory registration ──────────────────────────────────────────────────
void registerCosseratShellForceField(sofa::core::ObjectFactory* factory);

} // namespace sofa::component::cosserat::shell
