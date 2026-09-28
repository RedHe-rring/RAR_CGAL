cmake_minimum_required(VERSION 3.20)
project(RAR_CGAL VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(CGAL REQUIRED)
find_package(Eigen3 REQUIRED)

option(
    RAR_BUILD_FIELD_VIEWER
    "Build the interactive Polyscope sizing-field viewer"
    OFF
)

option(
    RAR_ENABLE_IPOPT
    "Enable IPOPT-based Chen sizing-field correction"
    ON
)

add_executable(rar_cgal
    src/main.cpp
    src/CSFSizingField.cpp
    src/ChenSizingCorrection.cpp
)

target_include_directories(rar_cgal PRIVATE include)
target_link_libraries(rar_cgal PRIVATE CGAL::CGAL Eigen3::Eigen)

if(RAR_ENABLE_IPOPT)
    find_path(
        IPOPT_INCLUDE_DIR
        NAMES IpTNLP.hpp
        PATH_SUFFIXES coin-or coin
    )
    find_library(
        IPOPT_LIBRARY
        NAMES ipopt libipopt
    )

    if(IPOPT_INCLUDE_DIR AND IPOPT_LIBRARY)
        target_include_directories(
            rar_cgal
            PRIVATE
                ${IPOPT_INCLUDE_DIR}
        )
        target_link_libraries(
            rar_cgal
            PRIVATE
                ${IPOPT_LIBRARY}
        )
        target_compile_definitions(
            rar_cgal
            PRIVATE
                RAR_HAS_IPOPT=1
        )
        message(
            STATUS
            "RAR+Chen enabled with IPOPT: ${IPOPT_LIBRARY}"
        )
    else()
        message(
            WARNING
            "IPOPT was not found. The executable will still build, "
            "but --field rar-chen will report an error at runtime. "
            "With vcpkg, install coin-or-ipopt[mumps] for the active triplet "
            "and reconfigure CMake."
        )
    endif()
endif()

if(MSVC)
    target_compile_options(rar_cgal PRIVATE /W4 /permissive-)
else()
    target_compile_options(rar_cgal PRIVATE -Wall -Wextra -Wpedantic)
endif()

include(CTest)
if(BUILD_TESTING)
    add_executable(rar_smoke_test tests/smoke.cpp)
    target_include_directories(rar_smoke_test PRIVATE include)
    target_link_libraries(rar_smoke_test PRIVATE CGAL::CGAL Eigen3::Eigen)
    add_test(NAME rar_smoke_test COMMAND rar_smoke_test)
endif()


if(RAR_BUILD_FIELD_VIEWER)
    include(FetchContent)

    FetchContent_Declare(
        polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG v2.6.1
        GIT_SHALLOW TRUE
        GIT_SUBMODULES ""
        PATCH_COMMAND
            git submodule update --init --recursive --depth 1
    )
    FetchContent_MakeAvailable(polyscope)

    add_executable(rar_field_viewer
        tools/field_viewer.cpp
    )
    target_link_libraries(
        rar_field_viewer
        PRIVATE
            polyscope
    )

    if(MSVC)
        target_compile_options(
            rar_field_viewer
            PRIVATE /W4 /permissive-
        )
    else()
        target_compile_options(
            rar_field_viewer
            PRIVATE -Wall -Wextra -Wpedantic
        )
    endif()
endif()
