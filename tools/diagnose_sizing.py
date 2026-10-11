#!/usr/bin/env python3
"""Offline diagnostics for RAR_CGAL field.ply + field.csv and a remeshed surface.

Dependencies: numpy, vtk, matplotlib
This estimates (not certifies) geometric deviation by sampled point-to-triangle distances.
The interpolated target size is a diagnostic proxy, not a query to the CGAL sizing API.
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def read_surface(path):
    path = Path(path)
    reader_by_ext = {
        ".ply": vtk.vtkPLYReader,
        ".obj": vtk.vtkOBJReader,
        ".stl": vtk.vtkSTLReader,
    }
    ext = path.suffix.lower()
    if ext not in reader_by_ext:
        raise ValueError(f"Unsupported mesh extension {ext!r}; use .ply, .obj or .stl")
    rd = reader_by_ext[ext]()
    rd.SetFileName(str(path))
    rd.Update()
    triangulate = vtk.vtkTriangleFilter()
    triangulate.SetInputConnection(rd.GetOutputPort())
    triangulate.Update()
    mesh = vtk.vtkPolyData()
    mesh.ShallowCopy(triangulate.GetOutput())
    if mesh.GetNumberOfPoints() == 0 or mesh.GetNumberOfPolys() == 0:
        raise ValueError(f"Empty or invalid triangle mesh: {path}")
    vertices = vtk_to_numpy(mesh.GetPoints().GetData()).astype(np.float64, copy=True)
    polygons = mesh.GetPolys()
    offsets = vtk_to_numpy(polygons.GetOffsetsArray())
    conn = vtk_to_numpy(polygons.GetConnectivityArray())
    if not np.all(np.diff(offsets) == 3):
        raise ValueError("Triangulation failed: found non-triangular cells")
    faces = conn.reshape((-1, 3)).astype(np.int64, copy=True)
    return mesh, vertices, faces


def read_field_csv(path, field_vertices):
    n = len(field_vertices)
    arrays = {}
    seen = np.zeros(n, dtype=bool)
    with open(path, "r", encoding="utf-8-sig", newline="") as f:
        rd = csv.DictReader(f)
        need = {"vertex_id", "x", "y", "z", "target_length"}
        if not need.issubset(rd.fieldnames or []):
            raise ValueError(f"Field CSV missing columns: {need - set(rd.fieldnames or [])}")
        properties = [k for k in rd.fieldnames if k not in ("vertex_id", "x", "y", "z")]
        arrays = {name: np.full(n, np.nan) for name in properties}
        coords = np.empty((n, 3), dtype=np.float64)
        for row in rd:
            i = int(row["vertex_id"])
            if i < 0 or i >= n or seen[i]:
                raise ValueError(f"Invalid or duplicate vertex_id: {i}")
            seen[i] = True
            coords[i] = (float(row["x"]), float(row["y"]), float(row["z"]))
            for name in properties:
                arrays[name][i] = float(row[name])
    if not seen.all():
        raise ValueError("Field CSV and PLY vertex counts / IDs disagree")
    diag = float(np.linalg.norm(np.ptp(field_vertices, axis=0)))
    err = float(np.max(np.linalg.norm(coords - field_vertices, axis=1)))
    if err > max(1e-9, diag * 1e-7):
        raise ValueError(f"Field CSV and PLY vertices misaligned (max difference: {err:.6g})")
    if not np.all(np.isfinite(arrays["target_length"])) or np.any(arrays["target_length"] <= 0):
        raise ValueError("target_length must be finite and strictly positive")
    return arrays


def build_locator(mesh):
    loc = vtk.vtkStaticCellLocator()
    loc.SetDataSet(mesh)
    loc.BuildLocator()
    return loc


def nearest(loc, point):
    closest = [0.0, 0.0, 0.0]
    cell = vtk.mutable(0)
    sub = vtk.mutable(0)
    dist2 = vtk.mutable(0.0)
    loc.FindClosestPoint(point.tolist(), closest, cell, sub, dist2)
    return int(cell), np.asarray(closest, dtype=np.float64), math.sqrt(max(0.0, float(dist2)))


def barycentric(point, triangle):
    a, b, c = triangle
    x, y = b - a, c - a
    p = point - a
    aa, ab, bb = np.dot(x, x), np.dot(x, y), np.dot(y, y)
    pa, pb = np.dot(p, x), np.dot(p, y)
    det = aa * bb - ab * ab
    if det <= max(1e-30, aa * bb * 1e-13):
        w = np.zeros(3)
        w[int(np.argmin(np.sum((triangle - point) ** 2, axis=1)))] = 1.0
        return w, True
    v = (bb * pa - ab * pb) / det
    w = (aa * pb - ab * pa) / det
    weights = np.maximum(np.array([1 - v - w, v, w]), 0.0)
    weights /= weights.sum()
    return weights, False


def subset_indices(n, limit, seed):
    if limit == 0 or n <= limit:
        return np.arange(n, dtype=np.int64)
    return np.sort(np.random.default_rng(seed).choice(n, size=limit, replace=False))


def sample_edges(faces, max_edges, seed):
    # If a cap is set, select faces uniformly before constructing edges to avoid
    # materializing all 3F edges of a multi-million-face mesh in Python.
    if max_edges and len(faces) * 3 > max_edges:
        selected = subset_indices(len(faces), max_edges, seed)
        tri = faces[selected]
    else:
        tri = faces
    edges = np.vstack([tri[:, (0, 1)], tri[:, (1, 2)], tri[:, (2, 0)]])
    edges.sort(axis=1)
    edges = np.unique(edges, axis=0)
    if max_edges and len(edges) > max_edges:
        edges = edges[subset_indices(len(edges), max_edges, seed + 1)]
    return edges


def stats(arr):
    a = np.asarray(arr, dtype=float)
    a = a[np.isfinite(a)]
    if not len(a):
        return {"count": 0}
    return {
        "count": int(len(a)), "mean": float(a.mean()),
        "p05": float(np.percentile(a, 5)), "median": float(np.median(a)),
        "p95": float(np.percentile(a, 95)), "max": float(a.max()),
    }


def write_csv(path, cols):
    keys = list(cols)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for row in zip(*(cols[k] for k in keys)):
            w.writerow([f"{float(x):.10g}" if isinstance(x, (np.floating, float)) else x for x in row])


def colors(values, lo, hi, cmap="coolwarm"):
    import matplotlib
    colormap = matplotlib.colormaps[cmap]
    values = np.asarray(values)
    m = np.isfinite(values)
    out = np.full((len(values), 3), 160, dtype=np.uint8)  # gray = not sampled
    if np.any(m):
        t = np.clip((values[m] - lo) / (hi - lo), 0, 1)
        out[m] = (np.asarray(colormap(t))[:, :3] * 255).astype(np.uint8)
    return out


def write_binary_colored_ply(path, xyz, rgb, faces=None):
    # Portable format: binary little-endian XYZ double + RGB uchar.
    path = Path(path)
    if len(xyz) != len(rgb):
        raise ValueError("Color count differs from vertex count")
    vtype = np.dtype([("x", "<f8"), ("y", "<f8"), ("z", "<f8"),
                      ("red", "u1"), ("green", "u1"), ("blue", "u1")])
    verts = np.empty(len(xyz), dtype=vtype)
    for k, c in zip(("x", "y", "z"), range(3)):
        verts[k] = xyz[:, c]
    for k, c in zip(("red", "green", "blue"), range(3)):
        verts[k] = rgb[:, c]
    nfaces = 0 if faces is None else len(faces)
    head = ("ply\nformat binary_little_endian 1.0\n"
            f"element vertex {len(xyz)}\nproperty double x\nproperty double y\n"
            "property double z\nproperty uchar red\nproperty uchar green\n"
            f"property uchar blue\nelement face {nfaces}\n"
            "property list uchar int vertex_indices\nend_header\n")
    with open(path, "wb") as fh:
        fh.write(head.encode("ascii"))
        verts.tofile(fh)
        if nfaces:
            if np.max(faces) >= np.iinfo(np.int32).max:
                raise ValueError("PLY 32-bit face vertex index overflow")
            ftype = np.dtype([("n", "u1"), ("a", "<i4"), ("b", "<i4"), ("c", "<i4")])
            records = np.empty(nfaces, dtype=ftype)
            records["n"] = 3
            for k, c in zip(("a", "b", "c"), range(3)):
                records[k] = faces[:, c]
            records.tofile(fh)


def aggregate_edge_vertex_values(nverts, edges, values):
    ids = edges.reshape(-1)
    sums = np.bincount(ids, weights=np.repeat(values, 2), minlength=nverts)
    count = np.bincount(ids, minlength=nverts)
    out = np.full(nverts, np.nan)
    np.divide(sums, count, out=out, where=count > 0)
    return out


def error_samples(src_vertices, src_faces, other_locator, max_vertices, max_faces, seed):
    """Sample source vertices and face centroids; closest distances to other surface.
    Return sample positions, distances, sparse per-vertex max error for visualization.
    """
    vi = subset_indices(len(src_vertices), max_vertices, seed)
    fi = subset_indices(len(src_faces), max_faces, seed + 3)
    p_vertices = src_vertices[vi]
    p_faces = src_vertices[src_faces[fi]].mean(axis=1)
    samples = np.vstack((p_vertices, p_faces))
    distances = np.empty(len(samples))
    for k, q in enumerate(samples):
        _, _, distances[k] = nearest(other_locator, q)
    vertex_map = np.full(len(src_vertices), np.nan)
    vertex_map[vi] = distances[:len(vi)]
    if len(fi):
        # Include face-centroid error as a conservative diagnostic on its vertices.
        indices = src_faces[fi].ravel()
        values = np.repeat(distances[len(vi):], 3)
        tmp = np.full(len(src_vertices), -np.inf)
        np.maximum.at(tmp, indices, values)
        valid = np.isfinite(tmp)
        vertex_map[valid] = np.fmax(vertex_map[valid], tmp[valid])
    return samples, distances, vertex_map


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--field-ply", type=Path, required=True)
    ap.add_argument("--field-csv", type=Path, required=True)
    ap.add_argument("--output-mesh", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--epsilon", type=float, help="Reference geometric error budget in mesh units (optional)")
    ap.add_argument("--max-edges", type=int, default=30000, help="0 = all unique edges")
    ap.add_argument("--max-vertices", type=int, default=15000, help="0 = all vertices per direction")
    ap.add_argument("--max-faces", type=int, default=15000, help="0 = all face centroids per direction")
    ap.add_argument("--match-factor", type=float, default=1.2, help="Visualization screening ratio; NOT CGAL split/collapse criteria")
    ap.add_argument("--candidate-budget", type=float, default=0.3,
                    help="Candidate local sampled error utilization cutoff, if epsilon is supplied")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    if args.epsilon is not None and args.epsilon <= 0:
        ap.error("--epsilon must be positive")
    if args.match_factor <= 1:
        ap.error("--match-factor must be > 1")
    if min(args.max_edges, args.max_vertices, args.max_faces) < 0:
        ap.error("sample caps must be >= 0; 0 means all")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    input_poly, input_v, input_f = read_surface(args.field_ply)
    output_poly, output_v, output_f = read_surface(args.output_mesh)
    fields = read_field_csv(args.field_csv, input_v)
    input_loc = build_locator(input_poly)
    output_loc = build_locator(output_poly)
    print(f"Background mesh: {len(input_v):,} V / {len(input_f):,} F")
    print(f"Remeshed output: {len(output_v):,} V / {len(output_f):,} F")
    print("Computing sampled edge length/field residuals...")
    edges = sample_edges(output_f, args.max_edges, args.seed)
    edge_vecs = output_v[edges[:, 1]] - output_v[edges[:, 0]]
    lengths = np.linalg.norm(edge_vecs, axis=1)
    mids = (output_v[edges[:, 1]] + output_v[edges[:, 0]]) / 2
    h = np.empty(len(edges))
    raw_h = np.full(len(edges), np.nan)
    field_distance = np.empty(len(edges))
    degenerate_count = 0
    for k, p in enumerate(mids):
        cid, cp, field_distance[k] = nearest(input_loc, p)
        ids = input_f[cid]
        weight, bad = barycentric(cp, input_v[ids])
        degenerate_count += int(bad)
        h[k] = float(np.dot(weight, fields["target_length"][ids]))
        if "raw_target_length" in fields:
            raw_h[k] = float(np.dot(weight, fields["raw_target_length"][ids]))
    ratios = lengths / h
    logratios = np.log(np.maximum(ratios, np.finfo(float).tiny))
    # Color both the underlying mesh (when edges sampled, unsampled = gray)
    # and a point cloud at every sampled edge midpoint (often easier to inspect).
    color_bound = math.log(args.match_factor) * 3
    d_vertex = aggregate_edge_vertex_values(len(output_v), edges, logratios)
    write_binary_colored_ply(args.out_dir / "size_mismatch_mesh.ply", output_v,
                             colors(d_vertex, -color_bound, color_bound), output_f)
    write_binary_colored_ply(args.out_dir / "size_mismatch_edge_points.ply", mids,
                             colors(logratios, -color_bound, color_bound))
    if "raw_target_length" in fields:
        raw = fields["raw_target_length"]
        if np.any(raw <= 0):
            raise ValueError("raw_target_length must be strictly positive")
        correction_ratio = fields["target_length"] / raw
        write_binary_colored_ply(args.out_dir / "chen_correction_input.ply", input_v,
                                 colors(correction_ratio, 0, 1, "viridis"), input_f)
    else:
        correction_ratio = None

    print("Computing sampled two-direction surface distances...")
    out_pts, out_d, out_map = error_samples(output_v, output_f, input_loc,
                                           args.max_vertices, args.max_faces, args.seed + 10)
    in_pts, in_d, in_map = error_samples(input_v, input_f, output_loc,
                                        args.max_vertices, args.max_faces, args.seed + 20)
    write_csv(args.out_dir / "error_output_to_input.csv", {
        "x": out_pts[:, 0], "y": out_pts[:, 1], "z": out_pts[:, 2], "distance": out_d,
    })
    write_csv(args.out_dir / "error_input_to_output.csv", {
        "x": in_pts[:, 0], "y": in_pts[:, 1], "z": in_pts[:, 2], "distance": in_d,
    })
    # Screening only: nearby sampled reverse distances plus edge-midpoint forward
    # distance. This is not a local Hausdorff bound. Require nearby input samples.
    tree = cKDTree(in_pts)
    knn = min(8, len(in_pts))
    neighbor_d, neighbor_id = tree.query(mids, k=knn)
    if knn == 1:
        neighbor_d, neighbor_id = neighbor_d[:, None], neighbor_id[:, None]
    vicinity = neighbor_d <= 3.0 * h[:, None]
    n_nearby = np.sum(vicinity, axis=1)
    near_errors = np.where(vicinity, in_d[neighbor_id], -np.inf)
    near_max = np.maximum(0.0, near_errors.max(axis=1))
    local_error_proxy = np.maximum(field_distance, near_max)
    candidate = np.zeros(len(edges), dtype=bool)
    candidate_short = np.zeros(len(edges), dtype=bool)
    candidate_small = np.zeros(len(edges), dtype=bool)
    if args.epsilon is not None:
        budget_ok = (n_nearby >= 3) & (local_error_proxy < args.candidate_budget * args.epsilon)
        candidate_short = budget_ok & (ratios < 1 / args.match_factor)
        # A second, independent hypothesis: even if the target is perfectly
        # realized, globally small elements in low-error regions may deserve
        # a controlled local-coarsening experiment. Median is just a screen.
        candidate_small = budget_ok & (lengths < np.median(lengths))
        candidate = candidate_short | candidate_small
    write_csv(args.out_dir / "edges.csv", {
        "x": mids[:, 0], "y": mids[:, 1], "z": mids[:, 2],
        "edge_length": lengths, "target_length_interpolated": h,
        "length_over_target": ratios, "log_length_over_target": logratios,
        "raw_target_length_interpolated": raw_h,
        "distance_to_background_surface": field_distance,
        "reverse_samples_within_3h": n_nearby,
        "local_error_proxy": local_error_proxy,
        "candidate_shorter_than_target_0_or_1": candidate_short.astype(int),
        "candidate_small_low_error_0_or_1": candidate_small.astype(int),
        "candidate_screening_only_0_or_1": candidate.astype(int),
    })
    write_binary_colored_ply(args.out_dir / "candidate_screening_edge_points.ply", mids,
                             np.where(candidate[:, None], np.array([239, 128, 44], dtype=np.uint8),
                                      np.array([160, 160, 160], dtype=np.uint8)))
    ref = args.epsilon if args.epsilon is not None else float(np.percentile(np.r_[in_d, out_d], 95))
    ref = max(ref, 1e-15)
    write_binary_colored_ply(args.out_dir / "error_output_to_input.ply", output_v,
                             colors(out_map / ref, 0, 2, "inferno"), output_f)
    write_binary_colored_ply(args.out_dir / "error_input_to_output.ply", input_v,
                             colors(in_map / ref, 0, 2, "inferno"), input_f)

    indicator = np.log(args.match_factor)
    summary = {
        "inputs": {"field_ply": str(args.field_ply), "field_csv": str(args.field_csv),
                   "output_mesh": str(args.output_mesh)},
        "sampling": {"edge_count": len(edges), "output_vertex_and_centroid_count": len(out_d),
                     "input_vertex_and_centroid_count": len(in_d), "seed": args.seed,
                     "max_edges": args.max_edges, "max_faces": args.max_faces,
                     "max_vertices": args.max_vertices},
        "mesh": {"input_vertices": len(input_v), "input_faces": len(input_f),
                 "output_vertices": len(output_v), "output_faces": len(output_f)},
        "size": {"ratio": stats(ratios), "log_ratio": stats(logratios),
                 "fraction_shorter_than_1_over_factor": float(np.mean(logratios < -indicator)),
                 "fraction_longer_than_factor": float(np.mean(logratios > indicator)),
                 "fraction_within_factor": float(np.mean(np.abs(logratios) <= indicator)),
                 "nearest_background_distance": stats(field_distance),
                 "degenerate_interpolation_faces": degenerate_count,
                 "candidate_shorter_than_target_count": int(np.sum(candidate_short)),
                 "candidate_small_low_error_count": int(np.sum(candidate_small)),
                 "candidate_screening_count_NOT_redundant_faces": int(np.sum(candidate)),
                 "candidate_screening_fraction_of_sampled_edges": float(np.mean(candidate)),
                 "chen_h_over_h0": stats(correction_ratio) if correction_ratio is not None else None},
        "geometric_error_sampled_NOT_Hausdorff": {
            "out_to_in": stats(out_d), "in_to_out": stats(in_d),
            "two_way_sampled_max": float(max(out_d.max(), in_d.max())),
            "epsilon": args.epsilon,
            "two_way_max_over_epsilon": float(max(out_d.max(), in_d.max()) / args.epsilon)
            if args.epsilon is not None else None,
        },
        "warnings": [
            "Target length uses barycentric interpolation on the closest input triangle; this is a proxy for the CGAL sizing API.",
            "Distance sampling is NOT a certified Hausdorff bound and can miss narrow features.",
            "Gray on a colored mesh means unsampled. Mismatch point cloud shows each sampled edge.",
            "A small length/target ratio or error does NOT prove mesh elements are redundant.",
            "Candidate screening uses only sparse nearby input->output samples within 3h and edge-midpoint output->input distance; it can miss severe errors.",
            "Internally RAR/PropertyMap sizing compares edge length with 4/5 or 4/3 times the MINIMUM endpoint sizing value; this tool uses interpolated background h at midpoint instead.",
        ],
    }
    bbox_d = max(float(np.linalg.norm(np.ptp(input_v, axis=0))), 1e-15)
    if np.median(field_distance) > 0.05 * bbox_d:
        summary["warnings"].append("Large median input/output distance: check units and alignment")
    with open(args.out_dir / "summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.hist(logratios, bins=60)
    ax.axvline(-indicator, linestyle="--")
    ax.axvline(+indicator, linestyle="--")
    ax.axvline(0, linestyle=":")
    ax.set_xlabel("log(actual edge length / interpolated target length)")
    ax.set_ylabel("Sampled edge count")
    ax.set_title("Sizing mismatch distribution")
    fig.tight_layout()
    fig.savefig(args.out_dir / "size_histogram.png", dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    ax.hist([in_d, out_d], bins=60, label=["input -> output", "output -> input"])
    if args.epsilon is not None:
        ax.axvline(args.epsilon, linestyle="--")
    ax.set_xlabel("Point-to-triangle surface distance (mesh units)")
    ax.set_ylabel("Sample count")
    ax.set_title("Sampled bidirectional geometric deviation")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out_dir / "error_histogram.png", dpi=150)
    plt.close(fig)
    if args.epsilon is not None:
        fig, ax = plt.subplots(figsize=(6.6, 4.2))
        ax.scatter(logratios, local_error_proxy / args.epsilon, s=3, alpha=0.25)
        ax.axhline(args.candidate_budget, linestyle="--")
        ax.axvline(-indicator, linestyle="--")
        ax.axvline(0, linestyle=":")
        ax.set_xlabel("log(actual / target)")
        ax.set_ylabel("Sampled local error proxy / epsilon")
        ax.set_title("Mismatch vs. local error (screening only)")
        fig.tight_layout()
        fig.savefig(args.out_dir / "mismatch_vs_error.png", dpi=150)
        plt.close(fig)
    print("Sampled edge ratio:", stats(ratios))
    print("Sampled two-way max geometric deviation:", summary["geometric_error_sampled_NOT_Hausdorff"]["two_way_sampled_max"])
    print(f"Written diagnostics to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
