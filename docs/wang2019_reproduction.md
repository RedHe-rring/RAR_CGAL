# Wang et al. (TVCG 2019): angle-driven isotropic remeshing

Reference: Yiqun Wang et al., *Isotropic Surface Remeshing without Large and Small Angles*, IEEE TVCG 25(7):2430–2442 (2019), DOI [10.1109/TVCG.2018.2837115](https://doi.org/10.1109/TVCG.2018.2837115).
[Official project page](https://nlpr.ia.ac.cn/ivc/project/acute3d/).

## Reproduction status (IMPORTANT)

This branch provides an **independent experimental postprocessor** written against CGAL 6.1.x; **it is not the author's source code or a faithful full-pipeline reproduction yet**. We start with the paper's main angle-driven concept, then independently check each operation against its figures and supplementary document. Do not report this executable as *Wang2019 exact* in comparison tables.

Implemented:
- Selection of large-angle triangles and insertion at the midpoint of their longest edge, with two adjacent quads re-triangulated using CGAL Euler operations.
- Selection of small-angle triangles and shortest-edge collapse, conditioned on Euler's link condition, simulated triangle normals/areas, and local minimum-angle improvement.
- Optional valence-reducing edge flips (target valence 6 for unconstrained interiors).
- Tangential Laplacian smoothing, with nearest-point projection to the **initial input mesh** using a CGAL AABB tree.
- Boundary / sharp-corner freezing (a conservative stand-in for the paper's feature-specific operations).
- Per-round angle/violation and vertex-count diagnostics, operation budgets, stopping checks, and synthetic C++ smoke test.

Still missing for a paper-level reproduction:
1. Match **all** cases of the vertex insertion figure, including alternative local re-triangulations and special boundary/feature configurations.
2. Match the exact selection priorities, parameter k, stopping conditions, and valence smoothing formulas with the paper and supplement.
3. Reproduce initial mesh construction / fixed vertex-count target N, and track accepted insertion/removal counts to control N as in the paper.
4. Implement the paper's own feature handling rather than freezing nearby vertices.
5. Add quantitative approximation-error / intersection checks; nearest-point projection and local normal guards do **not** guarantee global geometric fidelity or absence of self-intersection.
6. Run the actual paper's benchmark meshes and compare distributions/timings against published figures.

## Build (Windows, CGAL 6.1.1)

Use the existing project's CMake/vcpkg setup. Build **only** the standalone target; unlike the existing RAR+Chen executable this target does not need IPOPT at compile/link time:

\`\`\`powershell
cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=E:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake -DVCPKG_TARGET_TRIPLET=x64-windows -DRAR_ENABLE_IPOPT=OFF
cmake --build build --config Release --target wang2019_remesh wang2019_smoke_test
ctest --test-dir build -C Release -R wang2019_smoke_test --output-on-failure
\`\`\`

## Run as a CGAL/RAR postprocessor

\`\`\`powershell
.\build\Release\rar_cgal.exe .\data\plate2.ply .\data\plate2_cgal.ply --field cgal-adaptive --epsilon 0.001 --iterations 5
.\build\Release\wang2019_remesh.exe .\data\plate2_cgal.ply .\data\plate2_wang_proto.ply --min-angle 30 --max-angle 90 --rounds 10 --budget 0.02 --smooth 3
\`\`\`

Or directly provide a triangulated PLY/OBJ mesh. Run \`wang2019_remesh --help\` for options. This program **does not** invoke standard CGAL \`isotropic_remeshing()\` internally; it uses CGAL's Euler split/collapse/flip building blocks instead.

## Evaluation and interpretation

Compare **the same initial mesh** before/after postprocessing:
- N / number of triangles, runtime, smallest and largest interior angles, counts of angles <30° and >90°, 5% minimum-angle percentile and 95% aspect ratio.
- Surface-to-surface error (ideally symmetric HD and/or sampled CD), connected components and boundary edges, geometric self-intersections.
- Separate results for feature protection on and off.

The operation budget bounds how many splits/collapses are *attempted successfully* per pass and is **not** a target vertex count. The initial input mesh is the projection reference. A remeshed input yields projection to that remeshed surface, not to its own earlier ground-truth geometry.

Implementation location: \`include/rar/Wang2019Remesher.h\`, standalone CLI \`src/wang2019_main.cpp\`. No edits to existing RAR/Chen algorithm code.
