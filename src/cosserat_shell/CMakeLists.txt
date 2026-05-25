cmake_minimum_required(VERSION 3.12)
project(CosseratShell VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 20)

# ── Dependencies ──────────────────────────────────────────────────────────────
find_package(Sofa.Config REQUIRED)
sofa_find_package(Sofa.Component.StateContainer REQUIRED)
sofa_find_package(Sofa.Component.Topology.Container.Dynamic REQUIRED)

# ── Source files ──────────────────────────────────────────────────────────────
set(COSSERAT_SHELL_HEADERS
    # FEM
    fem/ShellElement.h
    fem/SE3ShellKinematics.h
    fem/ConstitutiveLaw.h
    fem/ShellElementComputer.h
    # ForceField
    forcefield/CosseratShellForceField.h
    forcefield/CosseratShellForceField.inl
)

set(COSSERAT_SHELL_SOURCES
    forcefield/CosseratShellForceField.cpp
)

# ── Library ───────────────────────────────────────────────────────────────────
add_library(CosseratShell SHARED ${COSSERAT_SHELL_HEADERS} ${COSSERAT_SHELL_SOURCES})

target_include_directories(CosseratShell
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>  # src/ → allows <cosserat_shell/...>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>  # plugin root → <liegroups/...>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(CosseratShell
    PUBLIC
        Sofa.Component.StateContainer
        Sofa.Component.Topology.Container.Dynamic
        Cosserat   # for liegroups (SE3, SO3)
)

# ── Install ───────────────────────────────────────────────────────────────────
sofa_create_package_with_targets(
    PACKAGE_NAME CosseratShell
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS CosseratShell AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "."
    INCLUDE_INSTALL_DIR CosseratShell
    RELOCATABLE "plugins"
)

# ── Tests ─────────────────────────────────────────────────────────────────────
cmake_dependent_option(COSSERAT_SHELL_BUILD_TESTS "Compile CosseratShell tests" ON
    "SOFA_BUILD_TESTS OR NOT DEFINED SOFA_BUILD_TESTS" OFF)

if(COSSERAT_SHELL_BUILD_TESTS)
    add_subdirectory(tests)
endif()
