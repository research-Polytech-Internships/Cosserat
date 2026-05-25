cmake_minimum_required(VERSION 3.12)
project(CosseratShell_tests)

find_package(Sofa.Testing REQUIRED)

set(SHELL_TEST_SOURCES
    ShellKinematicsTest.cpp
)

add_executable(${PROJECT_NAME} ${SHELL_TEST_SOURCES})
target_link_libraries(${PROJECT_NAME}
    CosseratShell
    Sofa.Testing
)

target_include_directories(${PROJECT_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../..    # src/ → <cosserat_shell/...>
        ${CMAKE_CURRENT_SOURCE_DIR}/../../.. # plugin root → <liegroups/...>
)

if(APPLE)
    target_link_libraries(${PROJECT_NAME} "-framework Cocoa")
endif()

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
