# Sizing Field Diagnostics for RAR_CGAL

**Scope:** Offline, read-only diagnosis of an exported *initial* sizing field and the *final* remeshed triangle surface. The tool **does not modify the mesh or target field**. It is designed for the `exp/chen-element-count-objective` branch and other branches that produce `field.ply` + `field.csv` with `vertex_id,x,y,z,target_length`.

## Installation (Windows PowerShell)

```powershell
cd E:\adaiso\compare\rar_csf_our\RAR_CGAL
python -m pip install numpy scipy vtk matplotlib
```

Place `diagnose_sizing.py` in `tools\diagnose_sizing.py` (or invoke its own absolute path).

## Run the remesher and export a field

```powershell
.\build\Release\rar_cgal.exe .\data\plate2.ply .\data\plate2_diagnostic_out.ply `
  --field rar-chen --epsilon 0.01 --min-edge 0.001 --max-edge 0.5 `
  --beta 1.2 --iterations 5 --relax-steps 3 `
  --export-field .\data\plate2_diagnostic_field
```

Expected artifacts: `data\plate2_diagnostic_out.ply`, `data\plate2_diagnostic_field.ply`, `data\plate2_diagnostic_field.csv`. If you do not specify a field prefix, your branch will create the `field.ply`/`field.csv` inside a same-named experiment folder instead.

## Analyze (PowerShell)

```powershell
python .\tools\diagnose_sizing.py `
  --field-ply .\data\plate2_diagnostic_field.ply `
  --field-csv .\data\plate2_diagnostic_field.csv `
  --output-mesh .\data\plate2_diagnostic_out.ply `
  --epsilon 0.01 `
  --out-dir .\data\plate2_diagnostic_results
```

Run all surface samples for a small mesh:

```powershell
python .\tools\diagnose_sizing.py `
  --field-ply .\data\plate2_diagnostic_field.ply `
  --field-csv .\data\plate2_diagnostic_field.csv `
  --output-mesh .\data\plate2_diagnostic_out.ply `
  --epsilon 0.01 --max-edges 0 --max-vertices 0 --max-faces 0 `
  --out-dir .\data\plate2_diagnostic_results_full
```

**Note:** The sampled geometric distance is only a diagnostic and is **not** a certified Hausdorff bound. `--epsilon` must use the same coordinate units. An RAR input epsilon does not become an absolute output error guarantee just because you pass it here.

## Files produced

| File | Content | Usage |
|---|---|---|
| `summary.json` | Mesh counts, edge-ratio percentiles, mismatch fractions, two-way sample errors | Compare runs |
| `edges.csv` | Midpoint, actual edge length, interpolated target, ratio/log ratio, local screening proxy | Numerical analysis |
| `size_mismatch_edge_points.ply` | Colored point cloud, one point per sampled output edge | Most faithful visualization of sampled edge values; open alongside original output mesh |
| `size_mismatch_mesh.ply` | Output mesh; incident sampled-edge log residual averaged to vertices | Smooth display, gray = unsampled, averages may cancel |
| `chen_correction_input.ply` | Background mesh colored by `h/h0` | See Chen's change to original field; only if raw_target_length exists |
| `error_output_to_input.ply` | Output colored by sampled output-to-input deviations | Find geometric departures |
| `error_input_to_output.ply` | Input colored by sampled input-to-output deviations | Find omitted features / displaced regions |
| `candidate_screening_edge_points.ply` | Orange = **screening-only** coarsening candidate, gray = no flag | Pick regions for controlled remeshing tests |
| `size_histogram.png` | log(actual length / target) histogram | Compare field realization |
| `error_histogram.png` | Sampled bidirectional distance histogram | See one-way vs two-way imbalance |
| `mismatch_vs_error.png` | Scatter of edge-size mismatch and local error proxy | Distinguish different phenomena |
| `error_input_to_output.csv`, `error_output_to_input.csv` | Sample locations and nearest-surface distances | Audit or produce more plots |

### Color meaning

For `size_mismatch_edge_points.ply` or `size_mismatch_mesh.ply`: **blue = actual edge shorter than interpolated target; neutral = similar; red = actual edge longer**. The diverging color range is `[-3 log(match-factor), +3 log(match-factor)]` and values outside are clipped. Default match factor is 1.2. It **does not** equal CGAL's internal split/collapse thresholds.

For `chen_correction_input.ply`: dark = stronger correction, yellow/light = little or no correction (using colormap `viridis`, values from 0 to 1).

For error maps: colors represent sampled distance divided by the supplied `epsilon`, clipped to [0,2]. If epsilon is omitted, the scale uses the 95th percentile of both sampled directions and should **not** be compared numerically between runs.

For `candidate_screening_edge_points.ply`: orange indicates **at least one** heuristic is met: (a) edge considerably shorter than target; (b) edge shorter than global median length. Both require at least 3 sampled input reference points within 3 target lengths and a nearby local error proxy less than `candidate-budget * epsilon` (default 0.3). These are **hypotheses**, NOT proofs that edges or triangles can safely be deleted.

## Essential mathematical definitions

- For each sampled output edge midpoint `m`, find the closest point `q` on the input field triangle mesh, find its triangle `T=(a,b,c)`, and compute `h(q)=lambda_a h_a+lambda_b h_b+lambda_c h_c`.
- Save `r = |edge| / h(q)`, `d=log(r)`. `r<1` indicates output edge is shorter than interpolated background field; `r>1` indicates longer.
- At initial field vertices, `c = h/h0` if `raw_target_length` exists.
- Sample vertices and triangle centroids on both surfaces, using exact closest-triangle queries for those sample locations. Report the maximum observed distances in both directions.
- The local error screening uses sampled reference-to-output error values within `3 h(q)` of the edge midpoint and the midpoint-to-reference distance. It is not a local Hausdorff error measure. A tiny/thin feature may go unsampled; increasing sampling is essential.

## Scientific limitations

1. **Closest-triangle correspondence is purely geometric.** It can select the wrong sheet around folds, thin structures, close-by components or sharp edges. The sampled output edge midpoint may be off-surface; `distance_to_background_surface` helps flag unreliable matches. Better production implementations should track projection/source triangle identifiers or enforce feature-aware correspondence.
2. **CGAL internal target sizes can change during local operations.** For the RAR/PropertyMap adapters in the examined branch, `is_too_long` compares edge length with `(4/3) * min(h_u,h_v)` and `is_too_short` with `(4/5) * min(h_u,h_v)`, based on **current** vertex sizing values. This Python tool instead interpolates the **initial background** field at the edge midpoint. Therefore it diagnoses the geometric *field-realization gap*, not the exact internal reasons why CGAL accepted an edge. To diagnose exact acceptance, add a C++ export of the **final sizing per output vertex/edge**, its constrained status, and split/collapse decisions.
3. **Limited sampling is non-exhaustive.** By default only up to 30,000 edges and 15,000 vertices + 15,000 face centroids *per direction* are sampled. In very large meshes, default edge sample is face-stratified rather than exactly uniform over all edges. All ratios and candidate percentages are sampled statistics; do not claim these are global counts.
4. **Sample error is not a Hausdorff certificate.** Geometry can deviate between sample points. Strengthen with more samples, triangle-level bounds, and the true input surface (CAD if available) for defensible research results.
5. **Small residual ≠ optimal field; small error ≠ redundant triangles.** Confirmation requires locally enlarging permitted target size, re-running the same CGAL pipeline, and checking two-way geometry error, quality and feature/topology invariants.
6. **Large meshes may be slow or memory-intensive.** VTK locators index both triangle meshes; even sample-limited runs still load both full meshes, and the colored mesh output writes every triangle. Begin with small/medium models. For multi-million-face meshes, a CGAL C++ diagnostics exporter and batch triangle distance evaluation is preferable.

## Method for research experiments

1. Pick a small input model with observed anomalous RAR+Chen field distribution.
2. Run three methods with the same input/projection/features/backend: `rar`, `rar-chen`, `cgal-adaptive-chen`.
3. For each run, inspect `chen_correction_input.ply`, `size_mismatch_edge_points.ply`, both directional error maps and `summary.json`.
4. Choose a patch with small geometry error but many small elements **and explicitly exclude features/boundaries for the first experiment**.
5. Increase the target-length cap by a factor such as 1.2 in that patch with a smoothly decaying transition to adjacent regions. Rerun the same remesher with identical iteration/projection/feature settings.
6. Record triangle counts, two-way error (prefer accurate computation), minimum angle P05, aspect ratio P95, constrained-feature retention and timing. If all remain within criteria and N falls, you have *empirical slack* in that patch; not yet a new optimization algorithm.

## Tests

```powershell
python .\tools\test_diagnostics.py
```

Two smoke tests generate tiny planar meshes: coincident surfaces must give zero distance, and planes separated by 0.1 must give a two-way distance of 0.1. They also verify ratio and raw/Chen field alignment.
