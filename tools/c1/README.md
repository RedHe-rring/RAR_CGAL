# C1: curvature-ordered error-aware sizing (research prototype)

C1 adds an **automatic curvature-bin target sizing** experiment to RAR_CGAL.
It leaves the existing RAR, CSF and Chen remeshing modes unchanged.
It is not the original CSF algorithm, and is **not Hausdorff-error certified**.

## Method

1. Estimate vertex normal-angle and normal-angle/edge-length proxies; smooth for curvature ordering.
2. Estimate per-vertex H from raw curvature, epsilon and a radius-guarded RAR-like heuristic.
3. Construct up to 64 base curvature bins, then merge consecutive bins into groups.
4. For each group, H_g = min H_i; find spatial group neighbors from original mesh edges.
5. Compute h_g = min_j H_j beta^distance(g,j) using multi-source Dijkstra.
6. Greedily split groups by estimated complexity C = sum A_g/h_g^2; stop when gain < 2% or K reaches 12.
7. Refine boundaries by one base-bin shift; export field/diagnostics and optionally remesh with CGAL.

All grouping and stopping parameters can be overridden by CLI. C1 groups may
contain disconnected surface regions. Beta controls *group-adjacent target
ratios*, not the continuous field gradient or actual triangle quality.

## Install and test

From RAR_CGAL root (Windows PowerShell):

```powershell
python -m pip install -r .\tools\c1\requirements.txt
python -m pytest .\tools\c1\tests -q
cmake --build .\build --config Release --target rar_cgal
```

## Generate the C1 sizing field

```powershell
python .\tools\c1\c1_pipeline.py --input .\data\wine_glass24.ply --output .\data\wine_glass24_c1 --epsilon 0.01 --min-edge 0.001 --max-edge 0.5 --beta 1.4 --bins 64 --max-groups 12
```

## Generate sizing field and remesh in one call

```powershell
python .\tools\c1\c1_pipeline.py --input .\data\wine_glass24.ply --output .\data\wine_glass24_c1 --epsilon 0.01 --beta 1.4 --remesh-exe .\build\Release\rar_cgal.exe --iterations 5 --relax-steps 3
```

## Use an existing field.csv directly

```powershell
.\build\Release\rar_cgal.exe .\data\wine_glass24.ply .\data\wine_c1.ply --field c1-external --sizing-csv .\data\wine_glass24_c1\field.csv --iterations 5 --relax-steps 3 --no-export-field
```

## Fixed-five-bin ablation

```powershell
python .\tools\c1\c1_pipeline.py --input .\data\wine_glass24.ply --output .\data\wine_fixed5 --epsilon 0.01 --mode fixed5
```

Use the existing --field csf for the CSF-paper-inspired baseline,
keeping the common CGAL remeshing backend. The fixed5 mode is not the
author's original five-bin histogram scheme; it only fixes K=5.

## Outputs

- field.csv: original vertex coordinates and index, target length, group, and curvature signals.
- groups.csv: group curvature extrema, area, H_g, h_g, and vertex count.
- summary.json: K, estimated number of faces, optimization history, conflicts and mesh metrics.
- target_length_heat.ply / groups_heat.ply: vertex color preview (unless --no-heatmap).
- c1_remeshed.ply / remesher.log: only when --remesh-exe is supplied.

## Known limitations

- H is a heuristic local size proxy, not a proven Hausdorff error bound.
- A min-edge conflict means the minimum size has overridden an even smaller local request.
- The optimizer assumes a target size constant inside every curvature group.
- Target values at newly split CGAL vertices come from neighbor interpolation.
- The C++ CSV adapter checks row count, vertex IDs and original coordinates.
- Python unit tests do not establish CGAL compatibility; compile and test with CGAL 6.1.1.
