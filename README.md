# RAR_CGAL

A CGAL-based adaptive isotropic remeshing testbed for studying RAR-style sizing fields while keeping a stable common remeshing backend.

## Implemented modes

The executable keeps several adaptive sizing modes side by side:

- `cgal-adaptive`: CGAL 6.1.x `Adaptive_sizing_field`.
- `cgal-adaptive-radius`: curvature-radius-guarded CGAL-style field used directly by the remesher, without Chen correction.
- `cgal-adaptive-chen`: the same CGAL adaptive field followed by the Chen gradient-constrained correction.
- `cgal-adaptive-radius-chen`: a curvature-radius-guarded CGAL-style field followed by the same Chen correction.
- `rar`: paper-oriented Dunyach et al. (2013) curvature/sizing field, passed to the same CGAL `isotropic_remeshing()` backend.
- `rar-chen`: the RAR field followed by the same Chen correction.
- `csf`: code-oriented curvature-smoothed-field sizing (Lv et al.), passed to the same CGAL backend.

This separation is intentional. It allows direct experiments on the sizing field without changing split/collapse/flip infrastructure.

## Current methodological status

### CGAL-Adaptive

CGAL computes local principal curvatures with
`interpolated_corrected_curvatures()`, converts them to a curvature-adaptive target length, and uses the result through the `PMPSizingField` interface.

### CGAL-Adaptive radius guard, with or without Chen

The `cgal-adaptive-radius` and `cgal-adaptive-radius-chen` modes share the same
initial sizing field:

```text
kappa = max(abs(k_min), abs(k_max))
r = 1 / kappa
h0 = sqrt(6 epsilon r - 3 epsilon^2),  epsilon <= 0.183 r
h0 = r,                                epsilon >  0.183 r
h0 <- clamp(h0, min_edge, max_edge)
```

The `0.183` factor is the requested rounded transition ratio; the exact point
where the formula equals `r` is approximately `0.183503`. At the implemented
threshold the formula is approximately `0.998765 r`, and immediately above it
the field falls back to `r`. For smaller epsilon the formula can produce
`h0 < r`; `r` is not a lower bound. For zero or non-finite curvature, `h0` is
`max_edge`. The `cgal-adaptive-radius` mode uses `h0` directly for splitting,
collapsing, and relaxation. It exports `curvature` and `target_length`.
The `cgal-adaptive-radius-chen` mode passes `h0` as the per-vertex upper bound
to the Chen gradient-constrained correction, so the corrected value cannot
exceed this raw value.

When `--epsilon` is omitted, epsilon is selected automatically for every
epsilon-based field. The program computes
`kappa = max(abs(k_min), abs(k_max))` at the input vertices, discards zero and
non-finite values, and takes the vertex-area-weighted 90th percentile. It then
uses

```text
epsilon = 0.183 / weighted_p90(kappa)
```

where each vertex weight is one third of the total area of its incident input
triangles. Consequently, only the highest-curvature approximately 10% of the
weighted surface enters the `h0 = r` branch. This makes the default scale with
the input geometry while avoiding the median rule's approximately 50% branch
transition. If no valid curvature sample exists, the fallback is
`epsilon = 0.001`. An explicit `--epsilon <value>` always overrides the
automatic value. The resolved epsilon is printed and used in automatic output
filenames.

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

### CSF-field-CGAL

The CSF sizing-field construction is now ported directly from the public author implementation in
`vvvwo/Adaptively-Isotropic-Remeshing`, commit
`53cbd9afd429e67bd1736e63594797cea702ceb0`
(`Mesh_Geometric.cpp` + `AdpIsotropic.cpp`).

The previous Dirichlet / sparse harmonic solve has been removed completely. The current field follows the author-code sequence:

```text
area-weighted vertex normals
-> mean normal-angle curvature N_Value
-> cotangent weights with tan(|angle|) clamped to [0.1, 10]
-> divide weights by edge length and normalize
-> 3 in-place smoothing sweeps, lambda = 0.5
-> 3 synchronous edge-length-weighted neighbor averages
-> 40-bin histogram
-> thresholds li1/li2/li3/li4 exactly from the author-code index rules
-> multipliers {1.8, 1.4, 1.0, 0.8, 0.6}
-> target_L = mean_edge_length * mesh_scale * multiplier
```

The field construction is author-code-faithful; the remeshing backend is still deliberately held fixed as CGAL `isotropic_remeshing()`. Therefore this mode is named `CSF-field-CGAL`, not an exact reproduction of the author's complete remesher.

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

### CGAL adaptive field + Chen correction

```bat
build\Release\rar_cgal.exe input.obj output_cgal_chen.obj ^
  --field cgal-adaptive-chen ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --beta 1.2 ^
  --iterations 5
```

The initial CGAL adaptive target lengths are copied to a vertex property map,
projected by the same Chen optimization used by `rar-chen`, and the corrected
field is then consumed by the same CGAL local remeshing backend. The exported
field contains `curvature`, `raw_target_length`, and `target_length`.

### RAR paper-oriented field

```bat
build\Release\rar_cgal.exe input.obj output_rar.obj ^
  --field rar ^
  --epsilon 0.001 ^
  --min-edge 0.001 ^
  --max-edge 0.5 ^
  --iterations 5
```

### CSF field

```bat
build\\Release\\rar_cgal.exe input.obj output_csf.obj ^
  --field csf ^
  --mesh-scale 1.0 ^
  --iterations 5
```

For the current CSF experiments, useful first-pass scales are `1.2`, `0.5`, and `0.1`.

### Projection diagnostic

```bat
build\Release\rar_cgal.exe input.obj output_no_project.obj ^
  --field rar ^
  --no-project
```

The RAR mode prints the initial curvature and target-length min/mean/max statistics to make abnormal fields easier to detect.

### Sharp-feature preservation

Sharp-feature detection and preservation are disabled by default. Enable them
with `--preserve-features`; the default dihedral-angle threshold is 50 degrees:

```bat
build\Release\rar_cgal.exe input.obj output_features.obj ^
  --field rar ^
  --preserve-features ^
  --feature-angle 50
```

Edges whose adjacent face normals differ by more than the threshold are passed
to CGAL as constrained edges. Feature endpoints and junctions are fixed, while
regular degree-2 feature vertices can relax along their feature polylines.
`--no-preserve-features` explicitly restores the default disabled state. Mesh
boundary edges remain subject to CGAL's built-in boundary handling regardless
of this switch.

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

the program writes the result next to the input mesh using a parameter-aware filename. Field-specific parameters are included only when they actually affect that method:

```text
RAR / CGAL-Adaptive:
input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar.obj

CSF:
input__scale-1p2__it-5__relax-3__proj-on__features-off__field-csf.obj
```

The naming skeleton is:

```text
input
+ field-specific parameters
+ common CGAL-backend parameters
+ field tag
```

Putting `field` last is intentional. Irrelevant parameters are not written into filenames; for example, CSF does not carry an `epsilon` tag.

The default field files are placed inside a folder whose name matches the output mesh stem, so the field files themselves can stay concise.

Decimal points are encoded as `p` so filenames remain shell-friendly
(for example, `0.001 -> 0p001`). The original input extension is preserved.

## Exporting the initial fields

Field export is **enabled by default** for all adaptive modes.

If the input is `input.obj` and the remeshing parameters are:

```text
field=rar
epsilon=0.001
min-edge=0.001
max-edge=0.5
iterations=5
relax-steps=3
projection=on
preserve-features=off
```

the mesh is written normally, while field artifacts are grouped in a same-named folder:

```text
input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar.obj

input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar/
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

For RAR and CGAL-Adaptive, the CSV columns are:

```text
vertex_id,x,y,z,curvature,target_length
```

For CSF, one additional diagnostic is exported:

```text
vertex_id,x,y,z,raw_curvature,curvature,target_length
```

Here `curvature` is the smoothed CSF value, while `raw_curvature` is the pre-smoothing normal-angle estimate. The PLY stores the same scalar properties and preserves triangle connectivity.

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
input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar.obj

input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar/
├── field.ply
└── field.csv
```

Color the target-length field directly inside that folder:

```bat
python tools\colorize_ply.py ^
  input__eps-0p001__lmin-0p001__lmax-0p5__it-5__relax-3__proj-on__features-off__field-rar\field.ply ^
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
                 +----------------------+----------------------+
                 |                      |                      |
       CGAL Adaptive field          RAR field              CSF field
       corrected principal       cotangent H/K       smoothed curvature
            curvature                                     + histogram
                 |                      |                      |
                 +----------------------+----------------------+
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
5. compare feature preservation on/off across sharp and smooth inputs.

## RAR + Chen sizing-field correction

The `rar-chen` field keeps the existing RAR curvature-to-size mapping and adds
the gradient-constrained sizing correction from Chen et al. before the same CGAL
isotropic remeshing backend is called. The `cgal-adaptive-chen` mode applies the
same correction to CGAL's adaptive field, which makes it possible to test Chen
as a field-independent post-processing step.

For each input vertex, the raw target length is denoted by `h0`, and `A_i` is
one third of the total area of the incident input triangles. The correction
solves

```text
minimize    sum_i A_i / h_i^2
subject to  ||grad h||_T^2 <= log(beta)^2   for every input triangle T
            min_edge <= h_i <= h0_i
```

The objective approximates the element count of an isotropic surface mesh.
Because every corrected value is bounded above by its raw value, the optimizer
can satisfy the gradation constraint only by retaining or reducing local target
lengths; it never coarsens beyond `h0`.

When `h0` violates the gradation bound, IPOPT uses an upper-bound continuation
with warm starts. It begins from a feasible scaled variation of `h0`, then
increases the allowed upper bound geometrically by a factor of 10 until the
full raw field is reached. If `h0` already satisfies the gradient bound, it is
used directly because it is the componentwise largest feasible field for this
objective. IPOPT variable-bound relaxation is disabled so clamping the returned
field to `[min_edge, h0]` does not introduce gradient violations on very small
triangles.

If a later continuation stage fails, the last successful, independently
validated stage is used as a fallback. Its field is still feasible for the
original bounds, but does not complete the optimization at `f=1` and may cause
substantially denser remeshing. Before field export/remeshing, a warning reports
the failed IPOPT status, failed fraction, retained fraction, gradient and
element-count objective. The final summary also identifies the fallback.
If the first stage fails, there is no successful checkpoint and the run still
reports an error. Non-finite values, material bound violations, or a failed
gradient check also cause an error before the corrected field is written.
Validation allows bound roundoff of `1e-12 * max(1, abs(h0_i))`, clamps it to the
stage bounds, and then checks `||grad h||^2 <= log(beta)^2 + 1e-6`, matching
IPOPT's configured acceptable constraint tolerance.

Thus the experiment changes only the initial sizing field; split/collapse/relax
still use the same `PMP::isotropic_remeshing()` backend as `--field rar`.

IPOPT is optional at build time but required to run `--field rar-chen`.
With vcpkg, install it for the same triplet used by CGAL:

```powershell
vcpkg install coin-or-ipopt[mumps]:x64-windows
```

Then reconfigure CMake so that the IPOPT include directory and library are
detected.

Example:

```powershell
build\Release\rar_cgal.exe model.obj ^
  --field rar-chen ^
  --epsilon 0.01 ^
  --min-edge 0.001 ^
  --max-edge 0.05 ^
  --beta 1.2 ^
  --iterations 5 ^
  --relax-steps 3
```

The exported RAR+Chen field contains three scalar properties:
`curvature`, `raw_target_length`, and `target_length`. The last one is the
Chen-corrected field actually consumed by the remesher.

