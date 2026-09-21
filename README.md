# RAR_CGAL

A CGAL-based adaptive isotropic remeshing baseline intended to provide a stable backend for reproducing and studying the RAR family of methods (Dunyach et al., 2013).

## Status

V1 uses CGAL 6.1.1's `Adaptive_sizing_field` together with `Polygon_mesh_processing::isotropic_remeshing()`. CGAL documents this field as curvature-adaptive and its remesher performs edge split, collapse, flip, tangential relaxation, and projection.

This V1 should be treated as **CGAL-Adaptive / RAR-family**, not yet as a paper-faithful reimplementation of every RAR detail. A later phase will add a custom `RARSizingField` while keeping CGAL as the topology/remeshing backend.

## Requirements

- CMake >= 3.20
- C++17 compiler
- CGAL 6.1.1 or compatible recent CGAL
- Eigen3

On Windows with vcpkg:

```powershell
vcpkg install cgal:x64-windows eigen3:x64-windows
```

Configure and build:

```powershell
cmake -S . -B build ^
  -DCMAKE_TOOLCHAIN_FILE=E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake ^
  -DVCPKG_TARGET_TRIPLET=x64-windows

cmake --build build --config Release
```

## Usage

```powershell
./build/Release/rar_cgal.exe input.obj output.obj `
  --epsilon 0.001 `
  --min-edge 0.001 `
  --max-edge 0.5 `
  --iterations 5 `
  --relax-steps 3
```

Disable projection for diagnostics:

```powershell
./build/Release/rar_cgal.exe input.obj output_no_project.obj --no-project
```

## V1 data flow

```text
Input mesh
  -> read / validate
  -> triangulate if necessary
  -> CGAL Adaptive_sizing_field
  -> CGAL isotropic_remeshing
       split
       collapse
       flip
       tangential relaxation
       projection (optional)
  -> output mesh
```

## Why this architecture?

The immediate goal is to remove custom split/collapse/connectivity maintenance as a confounding source of instability. Later experiments can compare RAR, CSF, and new sizing fields while keeping exactly the same CGAL remeshing backend.

## Next phase

- custom paper-oriented `RARSizingField`
- curvature/sizing-field export for visualization
- constrained feature/boundary handling
- edge-length/target-length statistics
- problem-case regression tests
