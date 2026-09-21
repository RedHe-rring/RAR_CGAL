# RAR_CGAL

A CGAL-based adaptive isotropic remeshing testbed for studying RAR-style sizing fields while keeping a stable common remeshing backend.

## Implemented modes

The executable now keeps two adaptive sizing modes side by side:

- `cgal-adaptive`: CGAL 6.1.x `Adaptive_sizing_field`.
- `rar`: paper-oriented Dunyach et al. (2013) curvature/sizing field, passed to the same CGAL `isotropic_remeshing()` backend.

This separation is intentional. It allows direct experiments on the sizing field without changing split/collapse/flip infrastructure.

## Current methodological status

### CGAL-Adaptive

CGAL computes local principal curvatures with
`interpolated_corrected_curvatures()`, converts them to a curvature-adaptive target length, and uses the result through the `PMPSizingField` interface.

### RAR-field-CGAL

The custom RAR field implements the paper's main sizing equations:

```text
H_i     = 1/2 ||Delta x_i||
K_i     = angle_deficit / A_i
kappa_i = H_i + sqrt(max(H_i^2 - K_i, 0))

L_i = sqrt(6 epsilon / kappa_i - 3 epsilon^2)
L_i <- clamp(L_i, L_min, L_max)

L(e) = min(L_i, L_j)
```

Implementation details:

- cotangent Laplace-Beltrami discretization
- mixed Voronoi area
- split when `|e| > 4/3 L(e)`
- collapse when `|e| < 4/5 L(e)`
- midpoint placement for split vertices
- new split-vertex sizing interpolated from the two current neighbors

The local mesh operations are still performed by CGAL.

**Important:** this mode is currently named `RAR-field-CGAL`, not yet `RAR-exact`. CGAL's adaptive tangential relaxation is not identical to Eq. (6) in the 2013 paper. The next phase will implement the paper relaxation separately so that this difference can also be studied instead of silently hidden.

## Requirements

- CMake >= 3.20
- C++17 compiler
- CGAL 6.1.1 or compatible recent CGAL
- Eigen3

On Windows with vcpkg:

```bat
vcpkg install cgal:x64-windows eigen3:x64-windows
```

Configure and build:

```bat
cmake -S . -B build ^
  -DCMAKE_TOOLCHAIN_FILE=E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake ^
  -DVCPKG_TARGET_TRIPLET=x64-windows

cmake --build build --config Release
```

Run tests:

```bat
ctest --test-dir build -C Release --output-on-failure
```

## Usage

### CGAL adaptive field

```bat
build\Release\rar_cgal.exe input.obj output_cgal.obj ^
  --field cgal-adaptive ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --iterations 5
```

### RAR paper-oriented field

```bat
build\Release\rar_cgal.exe input.obj output_rar.obj ^
  --field rar ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --iterations 5
```

### Projection diagnostic

```bat
build\Release\rar_cgal.exe input.obj output_no_project.obj ^
  --field rar ^
  --no-project
```

The RAR mode prints the initial curvature and target-length min/mean/max statistics to make abnormal fields easier to detect.

## Exporting the initial fields

Both adaptive modes can export the **pre-remeshing** curvature and target-length fields:

```bat
build\Release\rar_cgal.exe input.obj output_cgal.obj ^
  --field cgal-adaptive ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --export-field diagnostics/cgal

build\Release\rar_cgal.exe input.obj output_rar.obj ^
  --field rar ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --export-field diagnostics/rar
```

Each command produces:

```text
diagnostics/cgal_field.csv
diagnostics/cgal_field.ply

diagnostics/rar_field.csv
diagnostics/rar_field.ply
```

The CSV columns are:

```text
vertex_id,x,y,z,curvature,target_length
```

The PLY keeps the original triangle connectivity and stores two scalar vertex properties:

```text
curvature
target_length
```

For `cgal-adaptive`, the exported curvature is exactly the quantity used by CGAL's sizing formula: the maximum absolute value of the principal curvatures returned by `interpolated_corrected_curvatures()`.

For `rar`, the exported curvature is the cotangent / mixed-Voronoi estimate used by our RAR sizing implementation.

This makes the two fields directly comparable on the same input vertices before any split/collapse operation changes the mesh.

## Visualizing exported fields

The repository includes:

```text
tools/colorize_ply.py
```

It maps one scalar vertex property in an ASCII PLY file to a
blue -> cyan -> green -> yellow -> red heatmap while preserving the mesh faces.

By default, it uses the full scalar range:

```text
minimum value -> blue
maximum value -> red
```

There is no percentile clipping; this corresponds to the 0th-100th percentile range.

Color the RAR target-length field:

```bat
python tools\colorize_ply.py ^
  diagnostics\rar_field.ply ^
  diagnostics\rar_target_length_rgb.ply ^
  --property target_length
```

Color the RAR curvature field:

```bat
python tools\colorize_ply.py ^
  diagnostics\rar_field.ply ^
  diagnostics\rar_curvature_rgb.ply ^
  --property curvature
```

Do the same for CGAL:

```bat
python tools\colorize_ply.py ^
  diagnostics\cgal_field.ply ^
  diagnostics\cgal_target_length_rgb.ply ^
  --property target_length
```

Reverse the color direction when desired:

```bat
python tools\colorize_ply.py ^
  diagnostics\rar_field.ply ^
  diagnostics\rar_target_length_rgb_invert.ply ^
  --property target_length ^
  --invert
```

For direct side-by-side comparison, it is often better to force the same display range for both fields:

```bat
python tools\colorize_ply.py diagnostics\cgal_field.ply diagnostics\cgal_rgb.ply ^
  --property target_length --min 0.001 --max 0.05

python tools\colorize_ply.py diagnostics\rar_field.ply diagnostics\rar_rgb.ply ^
  --property target_length --min 0.001 --max 0.05
```

Using the same `--min/--max` avoids a misleading visualization where two different
numeric ranges are each independently stretched to the full heatmap.

## Architecture

```text
                            CGAL isotropic_remeshing
                                      |
                      +---------------+---------------+
                      |                               |
            CGAL Adaptive field              RAR paper field
       interpolated corrected curvature      cotangent H/K
                      |                               |
                      +---------------+---------------+
                                      |
                         same split/collapse/flip
                         same CGAL relaxation
                         same CGAL projection
```

This is the useful comparison for the current stage: change the field, hold the local remeshing implementation fixed.

## Next phase

1. compile and regression-test `--field rar` and field export on the Windows/CGAL 6.1.1 target;
2. visualize `curvature` and `target_length` for CGAL-Adaptive vs RAR on the same difficult meshes;
3. add edge `length / target_length` distribution statistics;
4. implement RAR Eq. (6) tangential relaxation as a separate mode;
5. add feature/boundary constraints after the smooth-surface baseline is stable.
