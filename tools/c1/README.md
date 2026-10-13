# C1: error-aware adaptive curvature binning (research prototype)

C1 extends the **curvature-bin sizing idea** of Lv et al. (TVCG 2024). It is
*not* a reproduction of the paper's exact remeshing implementation and does
*not* certify a Hausdorff error bound.

## Method

1. Estimate vertex-normal angle and angle/edge-length curvature proxies.
2. Smooth each proxy with nonnegative cotangent adjacency weights.
3. Use raw curvature, epsilon, and a radius-guarded RAR-like formula to
   estimate a per-vertex upper target size H. Values below min-edge are
   clipped and **counted as conflicts**, not silently certified.
4. Partition vertices along the smoothed curvature axis (up to 64 base bins).
5. For each candidate partition, assign H_g = min(H_i in G_g), compute
   area A_g, and construct the spatial adjacency graph between bins.
6. Solve h_g = min_j H_j beta**distance(g,j) by multi-source shortest path:
   componentwise largest feasible sizes for graph-adjacent ratio <= beta.
7. Greedily add curvature boundaries maximizing reduction of the
   equilateral-mesh complexity proxy sum A_g/h_g**2. Stop if gain < 2%,
   or K=12; locally refine boundaries by single-bin shifts.
8. Export field.csv, groups.csv, summary.json, and vertex-color PLY previews.
   Optionally run common CGAL isotropic remeshing using --field c1-external.

A group can contain spatially *disconnected* regions. Graph gradation is a
per-group proxy, not a continuous sizing-field gradient constraint.
New vertices created by CGAL receive neighborhood-interpolated targets.

## Requirements

Python >=3.10, numpy, scipy, trimesh, matplotlib; pytest for tests.
For remeshing: this repo's CGAL 6.1.1 executable, rebuilt after updating.

From the RAR_CGAL repository root:

`powershell
python -m pip install -r .\tools\c1\requirements.txt
python .\tools\c1\c1_pipeline.py --input .\data\wine_glass24.ply --output .\data\wine_glass24_c1 --epsilon 0.01 --beta 1.4 --min-edge 0.001 --max-edge 0.5
`

To run the CGAL backend as well, append:

`powershell
  --remesh-exe .\build\Release\rar_cgal.exe --iterations 5 --relax-steps 3
`

Or use the generated field directly:

`powershell
.\build\Release\rar_cgal.exe .\data\wine_glass24.ply .\data\wine_c1.ply --field c1-external --sizing-csv .\data\wine_glass24_c1\field.csv --iterations 5 --relax-steps 3 --no-export-field
`

Fixed-five-bin ablation (automatic sizing but fixed number of groups):

`powershell
python .\tools\c1\c1_pipeline.py --input .\data\wine_glass24.ply --output .\data\wine_fixed5 --epsilon 0.01 --mode fixed5
`

Tests:

`powershell
python -m pytest .\tools\c1\tests -q
`

## Outputs

- field.csv: ordered original vertex coordinates, target_length, group_id,
  curvature signals, local upper H.
- groups.csv: group curvature range, area, H_g, h_g, vertex count.
- summary.json: K, estimated faces, sizing range, conflicts, greedy history,
  and optional output metrics.
- target_length_heat.ply and groups_heat.ply: vertex-color visualizations.
- With --remesh-exe: c1_remeshed.ply and remesher.log.

**Limitations**: H is a local proxy (not a fidelity guarantee). beta bounds
target ratios between *groups sharing an original mesh edge*, not triangle
quality. The remeshing backend interpolates new vertex targets and uses the
shared CGAL split/collapse/flip/relaxation logic.

The C++ adapter validates CSV row count, index order and coordinates to
prevent applying targets to the wrong input mesh.
