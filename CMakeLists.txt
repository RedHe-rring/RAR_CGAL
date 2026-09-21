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

add_executable(rar_cgal
    src/main.cpp
)

target_include_directories(rar_cgal PRIVATE include)
target_link_libraries(rar_cgal PRIVATE CGAL::CGAL Eigen3::Eigen)

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
