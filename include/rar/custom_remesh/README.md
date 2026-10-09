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

## Remeshing snapshots

The behavior-preserving fork is committed as `939f90e`, with the annotated local
tag `custom-remesh-cgal-parity-baseline`. Optional observation hooks were added
after that checkpoint; they do not change the remeshing decisions.

```powershell
.\build\Release\rar_cgal.exe input.off result.off --field rar --remesh-backend custom --remesh-diagnostics
.\build\Release\rar_cgal.exe input.off result.off --field rar --remesh-backend custom --remesh-diagnostics --diagnostic-level operation
```

`--remesh-diagnostics` enables output; `--diagnostic-level stage|operation`
selects the granularity, defaulting to `stage`. The level alone does not enable
output. Diagnostics require `--remesh-backend custom`; the installed CGAL backend
rejects this option instead of silently switching implementations.

All snapshots are complete meshes in ASCII OFF with 17-digit coordinates.
They sit beside the field artifacts, without a new output root:

```text
result.off
result/
  field.ply
  field.csv
  remesh_diagnostics/
    run_0001/
      manifest.json
      events.jsonl
      status.json
      input.off
      iter_001/
        01_split.off
        02_collapse.off
        03_flip.off
        04_relax_001.off
        04_relax_002.off
        04_relax_003.off
        04_relax.off
        05_project.off
        result.off
      iter_002/...
```

With `--export-field some/folder/field`, the diagnostics root is
`some/folder/remesh_diagnostics/`. With `--no-export-field`, snapshots still use
that explicit prefix's parent, or the normal output-mesh stem directory if no
prefix was provided. Each run exclusively creates the next `run_NNNN` directory;
earlier runs are preserved. The actual path is printed to the console.

`stage` records the input at entry to remeshing (after any input triangulation
and field preparation), each completed split/collapse/flip/relax/project stage,
each complete relaxation sweep, and every completed remeshing iteration.
Disabled stages have a `skipped` event and no redundant snapshot. Zero remeshing
iterations produce only the input snapshot. Zero relaxation steps produce no
sweep or relaxation-stage snapshot.

`operation` additionally writes `iter_NNN/split/op_NNNNNN.off`, and equivalent
`collapse/` and `flip/` files after each accepted operation. Split snapshots are
written after retriangulation, constraint/sizing updates and queue maintenance.
A collapse snapshot includes its associated degeneracy repairs; intermediate
repair operations are not exposed as separate transactions. Reverted flips
are counted as rejected and do not create snapshots. Relaxation stays at full
sweep granularity, not individual vertex moves.

The optional `RemeshObserver` is passed explicitly through the local driver,
remesher and relaxation functions. It has no global state. The writer copies the
live mesh and compacts only the copy, preserving the live descriptors and queues.
Writing errors propagate to the CLI. `status.json` remains `incomplete` unless
all remeshing snapshots and event records have been written successfully;
`complete` refers to the remeshing diagnostics, not a later final mesh export.

`events.jsonl` records chronological event numbers, round/stage/operation numbers,
mesh counts, paths, acceptance counters and operation descriptors. `attempted`
counts candidates popped from the split/collapse queues, or edges visited by
the flip loop; `rejected = attempted - accepted`. Initial queue filtering is not
part of that count. IDs describe the live mesh at that event, not the compact
OFF vertex indices and not persistent identities across later operations.

`manifest.json` records the input file SHA-1, effective sizing/remeshing parameters,
paths, compiler/configuration, CGAL version, and the source commit/worktree-dirty
state captured at build time. Build metadata is refreshed on every CMake build;
untracked files also count as a dirty worktree. In-memory library calls with no
input file record a null checksum. These traces aid debugging; formal comparison
experiments still follow `notes/versioning_and_experiment_workflow.md`.

The operation level writes a complete mesh for every accepted edge operation,
so use the stage level for routine inspection and operation level for small
reproductions. Diagnostics are off by default.

## Validation

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

`rar_remesh_diagnostics_test` additionally compares diagnostics on/off for closed,
open and feature-constrained meshes at both levels, checks disabled stages, reads
every emitted snapshot and verifies topology/counts, verifies all three accepted
edge-operation kinds, repeat-run isolation and write-error propagation. CLI tests
also verify field-directory placement, input checksums, unchanged final output
and diagnostic-option validation.

After deliberately changing an operator, baseline parity may fail as expected.
Replace the relevant equality expectation with a documented experiment and
behavior checks; retain validity and feature-preservation checks. Keep a clean
commit of the behavior-preserving fork for future reference.
