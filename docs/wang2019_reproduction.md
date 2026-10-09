# Wang et al. (TVCG 2019): angle-driven isotropic remeshing

Reference: Yiqun Wang et al., *Isotropic Surface Remeshing without Large and Small Angles*, IEEE TVCG 25(7):2430–2442 (2019), DOI [10.1109/TVCG.2018.2837115](https://doi.org/10.1109/TVCG.2018.2837115).
[Official project page](https://nlpr.ia.ac.cn/ivc/project/acute3d/).

## Reproduction status (IMPORTANT)

This branch provides an **independent experimental postprocessor** written against CGAL 6.1.x; **it is not the author's source code or a faithful full-pipeline reproduction yet**. We start with the paper's main angle-driven concept, then independently check each operation against its figures and supplementary document. Do not report this executable as *Wang2019 exact* in comparison tables.

Implemented:
- Large-angle triangles: insert an **unprojected** edge midpoint, split adjacent quads, then evaluate four outer-edge flip candidates and apply the legal flip minimizing angle-deviation energy (may fall back to split only).
- Selection of small-angle triangles and shortest-edge collapse, conditioned on Euler's link condition, simulated triangle normals/areas, and local minimum-angle improvement.
- Two valence-reducing edge-flip passes per round (target valence 6 for unconstrained interiors), following the Algorithm 1 stage ordering.
- Area-weighted neighboring **triangle-centroid** tangential smoothing (Eq. 1-inspired, area weights are our choice), performed after both operation stages, with nearest-point projection to the **initial input mesh** using a CGAL AABB tree.
- Boundary / sharp-corner freezing (a conservative stand-in for the paper's feature-specific operations).
- Per-round angle/violation and vertex-count diagnostics, rejected-operation reason counters, per-operation V/F assertions, and synthetic C++ smoke tests.
- Default strict-N transactional rounds: if legal collapses cannot balance accepted insertions, restore that round's input mesh. This safeguard can be disabled with `--allow-drift`.

Still missing for a paper-level reproduction:
1. Match **all** cases of Fig. 4 precisely: currently only a simple split plus a locally best legal outer-edge flip is attempted; boundary/feature configurations and the exact affected-angle objective remain incomplete.
2. Match the paper's exact k adaptation, stopping conditions, smoothing weights, and efficient local 2–3-ring scheduling. Our default `--budget 0.02` is the fraction of bad triangles, while `--k` specifies an explicit stage budget.
3. Reproduce initial sizing-field construction for an explicitly requested target N. Strict-N currently preserves the input's live vertex count by rolling back an unbalanced round; it is not a faithful recovery strategy when collapse options run out.
4. Implement the paper's own feature handling rather than freezing nearby vertices.
5. Add quantitative approximation-error / intersection checks; nearest-point projection and local normal guards do **not** guarantee global geometric fidelity or absence of self-intersection.
6. Run the actual paper's benchmark meshes and compare distributions/timings against published figures.

## Build (Windows, CGAL 6.1.1)

Use the existing project's CMake/vcpkg setup. Build **only** the standalone target; unlike the existing RAR+Chen executable this target does not need IPOPT at compile/link time:

```powershell
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -DRAR_ENABLE_IPOPT=OFF
cmake --build build --config Release --target wang2019_remesh wang2019_smoke_test
ctest --test-dir build -C Release -R wang2019_smoke_test --output-on-failure
```

## Run as a CGAL/RAR postprocessor

```powershell
.\build\Release\rar_cgal.exe .\data\plate2.ply .\data\plate2_cgal.ply --field cgal-adaptive --epsilon 0.001 --iterations 5
.\build\Release\wang2019_remesh.exe .\data\plate2_cgal.ply .\data\plate2_wang_proto.ply --min-angle 30 --max-angle 90 --rounds 10 --budget 0.02 --smooth 3
```

Or directly provide a triangulated PLY/OBJ mesh. Run `wang2019_remesh --help` for options. This program **does not** invoke standard CGAL `isotropic_remeshing()` internally; it uses CGAL's Euler split/collapse/flip building blocks instead.

## Evaluation and interpretation

Compare **the same initial mesh** before/after postprocessing:
- N / number of triangles, runtime, smallest and largest interior angles, counts of angles <30° and >90°, 5% minimum-angle percentile and 95% aspect ratio.
- Surface-to-surface error (ideally symmetric HD and/or sampled CD), connected components and boundary edges, geometric self-intersections.
- Separate results for feature protection on and off.

The operation budget bounds how many splits/collapses are *attempted successfully* per pass and is **not** a target vertex count. The initial input mesh is the projection reference. A remeshed input yields projection to that remeshed surface, not to its own earlier ground-truth geometry.

Implementation location: `include/rar/Wang2019Remesher.h`, standalone CLI `src/wang2019_main.cpp`. No edits to existing RAR/Chen algorithm code.

## Diagnostics and safeguards (updated prototype)

```powershell
.\build\Release\wang2019_remesh.exe .\data\wine_glass24.ply .\data\wine_glass24_wang19.ply --min-angle 30 --max-angle 90 --rounds 10 --budget 0.02 --smooth 3
```

- Operation budgets are derived from the number of currently **bad triangles**, not total input vertices. Use `--k 50` to specify a fixed number explicitly.
- `--allow-drift` disables strict live-vertex-count rollback and permits independent insertions/collapses for debugging. Default mode commits only balanced rounds.
- `V` and `F` now count **live** Surface_mesh elements (`number_of_vertices()`, `number_of_faces()`), *not* the allocated descriptor slots returned by `num_vertices()` / `num_faces()`.
- Every accepted interior split must add +1 live vertex and +2 live triangles; every interior collapse must remove -1 live vertex and -2 live triangles. A mismatch throws before saving output.
- A detailed rejection breakdown (`boundary`, `feature`, `geometry`, `no-improvement`, `topology`) is printed per round.
- The default strict-N mechanism copies the full mesh once per round. It is a correctness-first implementation with significant memory/time overhead on large inputs.
