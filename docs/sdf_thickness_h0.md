# SDF thickness initial sizing field

This worktree provides `--field cgal-adaptive-sdf-thickness` and
`--field cgal-adaptive-sdf-thickness-chen`. Both construct the same initial
target edge-length field and use the existing `PropertyMapSizingField` with
the selected `cgal` or `custom` remeshing backend. The second mode applies
Chen's gradient-constrained correction to `h0` before remeshing; the original
mode does not call Chen.

## Workspace and baseline

- Worktree: `E:/adaiso/compare/rar_csf_our/RAR_CGAL_sdf_thickness`
- Branch: `exp/sdf-thickness-h0`
- Parent checkpoint: `7aafa08`, annotated tag `custom-remesh-diagnostics-baseline`
- Parent validation: Release build, CTest 5/5 (including backend/observer parity).
- Original custom worktree and its local data remain at `RAR_CGAL_custom_remesh`.
- Radius ablation commits on `exp/chen-element-count-objective` have not been
  merged. The existing `cgal-adaptive-radius` mode is available in this branch.

This stage implements and tests the modes. No Wine/Camel/Shark runs, parameter
sweeps, quality comparisons, or experiment results are included.

## Definition

For each initial vertex, with `kappa = max(abs(k1), abs(k2))`:

```
epsilon_local = min(epsilon, 0.183 / kappa, thickness / 4)
h_formula    = sqrt(6 * epsilon_local / kappa - 3 * epsilon_local^2)
h0           = clamp(h_formula, min_edge, max_edge)
```

Curvature uses CGAL `interpolated_corrected_curvatures`, matching the estimator
used by the existing Radius mode. For exactly zero curvature, the formula's
limit is infinite and `h0 = max_edge`. Negative or non-finite curvature is an
error. The factored numerical evaluation avoids squaring large length values.
Unrepresentable local tolerances are reported as errors.

Raw SDF and the scalar sizing formula scale with physical units. The unchanged
CGAL principal-curvature estimator is not strictly scale covariant: its
eigensolve includes a finite normal-direction penalty depending on mean edge
length. On the synthetic slab scaled by 7, the maximum relative curvature and
integrated `h0` residuals are about `1.87e-5` in this installed CGAL version.
The test records these residuals rather than assuming a scale-invariance
guarantee from CGAL. It strictly checks raw thickness scaling, scalar formula
scaling, and integrated targets against the actual estimated curvature.

The literal 0.183 is retained from the proposed formula. When that term limits
epsilon, `h_formula = sqrt(6*0.183 - 3*0.183^2)/kappa`, approximately
`0.998766/kappa`. The old Radius mode returns exactly `1/kappa` above its
threshold. Thus disabling the thickness effect does not make these two modes
identical. The old Radius implementation is retained for comparison.

The thickness ratio is fixed at `1/4` and recorded in the SDF metadata.
`thickness/4` caps geometric tolerance, not edge length. In particular, a flat
thin region can still have `h0 = max_edge`. This formula and the remesher do not
provide a certified approximation-error, wall-preservation, or no-intersection
guarantee. Applying `min_edge` can also raise the target above the formula.

In the `-chen` mode, the existing Chen solver receives this clamped `h0`,
the selected `--beta` (default `1.2`), and `min_edge`. It reduces target
lengths as needed to meet its face-gradient limit and leaves the original
mode untouched. The solver requires an IPOPT-enabled build. The initial
SDF and curvature values are computed once; Chen does not recompute them.

## SDF and input contract

`CGAL::sdf_values(mesh, face_map, cone_angle, rays, false)` supplies raw SDF.
CGAL performs ray sampling, intersection, and ray-length outlier removal.
Postprocessing is disabled because it includes normalization to `[0,1]`.
No segmentation or additional bilateral smoothing is performed.

SDF is an estimate of local object diameter, in the same length units as the
mesh. It is not the exact minimum wall thickness. Transfer to a vertex is:

```
t(v) = sum(area(f) * sdf(f), f incident on v) / sum(area(f), f incident on v)
```

This averaging is a declared first-version choice; it can soften abrupt
thickness changes near rims and corners. Values are computed once on the input.
During remeshing the existing adapter assigns new split vertices the mean of
their neighbors' target lengths. SDF and curvature are not recomputed each
iteration, and initial diagnostic properties are not evolving fields.

The first version requires a nonempty, closed triangle mesh with finite
coordinates, manifold vertices, nondegenerate faces, no self-intersections,
and orientation bounding a positive-volume solid. The CLI can triangulate
polygonal input using its existing path. It does not repair or reorient it.
SDF values of `-1`, zero, or non-finite values cause an error reporting the
number of invalid faces. There is no silent fallback or imputation.

The previous Radius report identified self-intersections in the Camel/Shark
inputs; those inputs will need separate attention before using this strict mode.

## Parameters and output

- `--sdf-rays N`: positive integer, default 25.
- `--sdf-cone-angle DEGREES`: full cone opening angle, strictly in `(0,180)`,
  default 120. The implementation converts degrees to radians for CGAL.
- The two SDF options are rejected with other field modes.
- `--beta B` applies to the `-chen` mode and must be finite and greater than 1.
- Existing epsilon selection and min/max-edge options apply. Supplying epsilon
  explicitly makes future comparisons easier to reproduce.
- `--iterations 0` constructs and exports the initial field without remeshing.
- `--remesh-backend custom --remesh-diagnostics` enables inherited snapshots.
- Automatic output names use compact parameter labels such as
  `__sdf-h0__e-...__mn-...__mx-...__r25__a120__i5__q3__p1__f0__b-custom`.
  The Chen mode uses `__sdf-h0-chen` and adds `__beta-...`.
  Numeric values retain the existing 10-significant-digit formatter; enabling
  feature preservation adds its angle. This shortens the result directory too,
  since field files and remesh snapshots live below it.

With field prefix `result/field`, initial diagnostics are:

| File | Contents |
| --- | --- |
| `field.csv`, `field.ply` | Initial vertex curvature, target length, thickness, vertex area, local epsilon, dimensionless epsilon*kappa, limiter and clipping flags. In Chen mode, `target_length` is corrected and `raw_target_length` is the original `h0`. |
| `field_sdf_faces.csv` | Original face descriptor, centroid, area, raw thickness, validity flag |
| `field_sdf.json` | SDF parameters, units, transfer/invalid-value policy, Chen status and correction summary, area-weighted thickness-limited fraction |

Limiter flags include ties; they do not form mutually exclusive categories.
Clipping flags are independent of epsilon-limit flags. Face validity is 1 on
successful exports because any invalid raw SDF stops field construction.
`--no-export-field` suppresses these field files; stage snapshots can still be
requested independently. The remesh snapshot manifest includes SDF parameters.

## Build and functional validation

Run in this worktree (the build directory belongs to this checkout):

```powershell
cmake -S . -B build-chen -DCMAKE_TOOLCHAIN_FILE=E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DRAR_BUILD_FIELD_VIEWER=OFF -DRAR_ENABLE_IPOPT=ON
cmake --build build-chen --config Release --target rar_cgal rar_sdf_thickness_test rar_smoke_test rar_chen_continuation_test rar_custom_remesh_test -j 2
ctest --test-dir build-chen -C Release --output-on-failure
```

The new tests use small synthetic fixtures to check the formula and clipping,
flat limit, physical scaling, analytic slab thickness with near-normal rays,
face-to-vertex transfer, invalid geometry/parameters, initial exports, CLI
naming, automatic epsilon, and parity between backends and observer settings.
The CLI tests also check Chen mode selection, initial and corrected exports,
metadata, bounds, and backend parity.
These are functional checks, not comparative remeshing experiments.

Reference: [CGAL 6.1.1 SDF manual](https://doc.cgal.org/6.1.1/Surface_mesh_segmentation/index.html).
