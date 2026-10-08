# Local remeshing operators

This directory contains a local fork of the installed CGAL implementation,
initially preserving its algorithm. The exact source version is
`6.1.1-I-900` (`CGAL_VERSION_NR=1060110900`), not an unqualified 6.1 release.
`SOURCE.json` records the original paths and SHA-256 hashes of the source bytes.
The original copyright and SPDX notices are retained in each header.

## Entry points and editing locations

| File | What to edit |
| --- | --- |
| `Remesh.h` | Initialization, named parameters, iteration loop and operator order |
| `IncrementalRemesher.h` | `split_long_edges`, `collapse_short_edges`, `flip_edges_for_valence_and_shape`, constraints and `project_to_surface` |
| `TangentialRelaxation.h` | Tangential displacement, sizing weights, constraint relaxation and normal checks |

The entry point is `rar::custom_remesh::isotropic_remeshing()`. All five project
call sites use `rar::dispatch_remeshing()` in `../RemeshBackend.h`. The sizing
field is independent of the backend.

The initial fork changes namespaces, header guards, local includes and call
qualification only. Namespace-local using directives retain upstream dependency
lookup without importing names into the enclosing `rar` namespace. Calls to
local remeshing/relaxation overloads are explicitly qualified to prevent ADL
from selecting installed CGAL operators. `Uniform_sizing_field` stays the CGAL
type, preserving the uniform-weight relaxation branch. The flip geometry helpers
`should_flip` and `is_badly_shaped` remain CGAL dependencies; edit/fork those too
if an experiment needs to change their predicates.

CGAL mesh types, Euler operations, geometric predicates, property maps, AABB
trees and other utilities are still dependencies. This is not a standalone
replacement for CGAL. Review compatibility and rerun parity checks when changing
the installed CGAL version.

One call constructs one remesher and reference-surface structure for the entire
iteration loop. Do not replace this with repeated calls for single iterations:
that can change the reference surface used for projection.

## Usage

From the `RAR_CGAL_custom_remesh` worktree:

```powershell
.\build\Release\rar_cgal.exe input.off --field rar --remesh-backend cgal
.\build\Release\rar_cgal.exe input.off --field rar --remesh-backend custom
```

`cgal` is the default. Automatic mesh filenames include `__backend-cgal` or
`__backend-custom`; diagnostic directories inherit that distinction. Explicit
output paths and `--export-field` prefixes remain exactly what the caller asks
for, so choose different paths when running a comparison.

```powershell
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

`rar_custom_remesh_test` runs 24 baseline parity cases: closed ellipsoid, open
nonuniform grid, sharp cube; uniform and property-map sizing; projection on/off;
one iteration without relaxation and five iterations with three relaxations.
It compares counts, halfedge connectivity, positions (distance tolerance 1e-12
on unit-scale fixtures), and feature edge flags, and checks fixed corners and
mesh validity. It also checks empty input and backend-specific naming. These
checks remain active in Release builds.

`rar_remesh_backend_cli_test` compares exported OFF files through the actual CLI
for all non-Chen fields, and all Chen fields when IPOPT is available. It also
checks default selection, automatic naming and invalid/missing backend values.
Generated test meshes stay in the worktree's ignored build directory.

After deliberately changing an operator, baseline parity may fail as expected.
Replace the relevant equality expectation with a documented experiment and
behavior checks; retain validity and feature-preservation checks. Keep a clean
commit of the behavior-preserving fork for future reference.
