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

## Automatic output naming

The output mesh argument is optional.

If you provide an output path, it is used exactly as given:

```bat
build\Release\rar_cgal.exe input.obj my_result.obj ^
  --field rar ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --iterations 5
```

If you omit the output path:

```bat
build\Release\rar_cgal.exe input.obj ^
  --field rar ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --iterations 5 ^
  --relax-steps 3
```

the program writes the result next to the input mesh using a parameter-aware filename such as:

```text
input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar.obj
```

The automatic mesh filename records the common experiment parameters first and the field type last:

```text
epsilon
min edge length
max edge length
iteration count
relaxation-step count
projection on/off
field
```

Putting `field` last is intentional: when filenames are sorted lexicographically, results with identical remeshing parameters but different sizing fields stay adjacent.

The default field files are placed inside a folder whose name matches the output mesh stem, so the field files themselves can stay concise.

Decimal points are encoded as `p` so filenames remain shell-friendly
(for example, `0.001 -> 0p001`). The original input extension is preserved.

## Exporting the initial fields

Field export is **enabled by default** for both adaptive modes.

If the input is `input.obj` and the remeshing parameters are:

```text
field=rar
epsilon=0.001
min-edge=0.001
max-edge=0.5
iterations=5
relax-steps=3
projection=on
```

the mesh is written normally, while field artifacts are grouped in a same-named folder:

```text
input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar.obj

input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar/
├── field.ply
└── field.csv
```

This keeps the experiment parameters in one place—the model/folder name—without repeating them on every diagnostic file.

You do not need to pass `--export-field`.

To disable field export:

```bat
build\Release\rar_cgal.exe input.obj ^
  --field rar ^
  --no-export-field
```

To override the field-output stem manually:

```bat
build\Release\rar_cgal.exe input.obj ^
  --field rar ^
  --export-field diagnostics/custom_rar
```

which produces:

```text
diagnostics/custom_rar.csv
diagnostics/custom_rar.ply
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

With the default field layout, an experiment now looks like:

```text
input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar.obj

input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar/
├── field.ply
└── field.csv
```

Color the target-length field directly inside that folder:

```bat
python tools\colorize_ply.py ^
  input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__field-rar\field.ply ^
  --property target_length ^
  --invert
```

Because the input is now simply `field.ply`, the auto-generated visualization name is also concise:

```text
field__color-target_length__invert-on__vmin-0p001__vmax-0p5.ply
```

and it remains inside the experiment folder.

For curvature:

```bat
python tools\colorize_ply.py ^
  <experiment-folder>\field.ply ^
  --property curvature
```

For direct CGAL-vs-RAR comparison, force the same display range:

```bat
python tools\colorize_ply.py <cgal-folder>\field.ply ^
  --property target_length ^
  --min 0.001 ^
  --max 0.05 ^
  --invert

python tools\colorize_ply.py <rar-folder>\field.ply ^
  --property target_length ^
  --min 0.001 ^
  --max 0.05 ^
  --invert
```

The visualization filename records the scalar property, invert state, and actual display range, while the parent folder records the remeshing parameters.

If you want a custom visualization filename, provide it as the second positional argument:

```bat
python tools\colorize_ply.py <experiment-folder>\field.ply my_visualization.ply ^
  --property target_length ^
  --invert
```

## Interactive field viewer

For day-to-day field inspection, the repository also provides a small C++ viewer:

```text
rar_field_viewer.exe
```

It uses Polyscope and opens the `field.ply` generated by `rar_cgal`. The main remeshing executable does not depend on the viewer; the viewer is an optional build target so normal experiments remain lightweight.

Configure once with the viewer enabled:

```bat
cmake -S . -B build ^
  -DCMAKE_TOOLCHAIN_FILE=E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake ^
  -DVCPKG_TARGET_TRIPLET=x64-windows ^
  -DRAR_BUILD_FIELD_VIEWER=ON
```

Polyscope is fetched automatically by CMake the first time the viewer is configured.

Build:

```bat
cmake --build build --config Release --target rar_field_viewer
```

Open a field:

```bat
build\Release\rar_field_viewer.exe ^
  <experiment-folder>\field.ply
```

On Windows, a `field.ply` can also be dragged onto `rar_field_viewer.exe`.

The viewer registers three switchable scalar modes:

```text
target_length
curvature
refinement_demand
```

`refinement_demand` is a normalized inverted target-length visualization:

```text
1 = smallest target length = strongest refinement demand
0 = largest target length  = weakest refinement demand
```

It is enabled by default because this direction is usually the most intuitive for inspecting where the remesher wants more triangles.

In the Polyscope quantity panel you can interactively:

- switch among the three scalar modes;
- drag or type the color-map minimum and maximum;
- change the colormap;
- show/hide the colorbar and histogram;
- inspect values by clicking mesh elements.

In the mesh options you can also adjust edge width, flat/smooth shading, and back-face display. These controls are useful for distinguishing real holes from flipped/back-facing triangles.

The Python `tools/colorize_ply.py` utility remains available for scripted figure export, but it is no longer required for normal interactive inspection.

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
