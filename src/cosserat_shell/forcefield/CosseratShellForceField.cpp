/******************************************************************************
 * Cosserat Shell — CosseratShellForceField explicit instantiation + factory *
 ******************************************************************************/
#define SOFA_COSSERAT_SHELL_CPP_CosseratShellForceField

#include <cosserat_shell/forcefield/CosseratShellForceField.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa::component::cosserat::shell {

// ── Explicit instantiation ─────────────────────────────────────────────────
template class SOFA_COSSERAT_API CosseratShellForceField<sofa::defaulttype::Vec6Types>;

// ── SOFA factory registration ─────────────────────────────────────────────
void registerCosseratShellForceField(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData(
        "Geometrically-exact Cosserat shell force field (SE(3) Lie group).\n"
        "\n"
        "Implements the internal virtual work principle for a Cosserat shell:\n"
        "  G_int = ∫_A Σ_α ⟨S^α, ∂κ/∂ξ^α + ad_{ζ_{tα}} κ⟩ j₀ dA\n"
        "\n"
        "Element: 4-node isoparametric Q4 with centroid strain evaluation\n"
        "(singularity-free, no shear locking).\n"
        "\n"
        "DOF per node: [ω_x, ω_y, ω_z, v_x, v_y, v_z] (Vec6 MechanicalState).\n"
        "\n"
        "Required Data:\n"
        "  - youngModulus: E [Pa]\n"
        "  - poissonRatio: ν ∈ [0, 0.5)\n"
        "  - thickness:    h [m]\n"
        "  - quads:        quad connectivity [[n0,n1,n2,n3], ...]")
        .add<CosseratShellForceField<sofa::defaulttype::Vec6Types>>(true));
}

} // namespace sofa::component::cosserat::shell
