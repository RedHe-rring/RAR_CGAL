#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Remeshing metrics: Hu 2017 quality metrics + Chamfer Distance (CD)

Metrics:
    #V                  output vertex count
    Q_min               minimum triangle quality
    Q_avg               average triangle quality
    theta_min_deg       global minimum interior angle
    theta_min_avg_deg   average of per-face minimum angles
    theta_max_deg       global maximum interior angle
    tri_lt_30_pct       percentage of triangles with min angle < 30 deg
    V567_pct            percentage of vertices with valence 5, 6, or 7
    CD                   normalized symmetric squared Chamfer Distance
    CD_x1e4             CD * 1e4, convenient for tables
    runtime_s           remeshing runtime supplied by user (optional)

CD definition:
    CD(A,B) =
        1/2 * mean_{a in A} min_{b in B} ||a-b||^2
      + 1/2 * mean_{b in B} min_{a in A} ||b-a||^2

A and B are independently sampled, area-uniform point sets on the two
triangle surfaces. The reported `CD` is normalized by the squared
bounding-box diagonal of the INPUT mesh:

    CD = CD_abs2 / D_bbox^2

This makes CD comparable across models with different scales.

Dependencies:
    pip install numpy trimesh scipy

Example:
    python hu2017_metrics.py input.ply remesh.ply

    python hu2017_metrics.py input.ply remesh.ply \
        --cd-samples 100000 \
        --runtime-seconds 31.96 \
        --csv results.csv --append-csv
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import trimesh
from scipy.spatial import cKDTree


EPS = 1e-15


def progress(message: str) -> None:
    now = time.strftime("%H:%M:%S")
    print(f"[{now}] {message}", flush=True)


def load_triangle_mesh(path: str | Path) -> trimesh.Trimesh:
    """Load a triangular mesh without cleanup that could change #V."""
    obj = trimesh.load(str(path), process=False)

    if isinstance(obj, trimesh.Scene):
        geoms = [
            g for g in obj.geometry.values()
            if isinstance(g, trimesh.Trimesh) and len(g.faces) > 0
        ]
        if not geoms:
            raise ValueError(f"No triangle mesh found in scene: {path}")
        mesh = trimesh.util.concatenate(geoms)
    elif isinstance(obj, trimesh.Trimesh):
        mesh = obj
    else:
        raise TypeError(f"Unsupported mesh type from {path}: {type(obj)}")

    if mesh.faces.ndim != 2 or mesh.faces.shape[1] != 3:
        raise ValueError(f"{path} is not a triangular mesh.")

    return mesh


def _angle_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise angle in degrees. Degenerate vectors yield 0 degrees."""
    na = np.linalg.norm(a, axis=1)
    nb = np.linalg.norm(b, axis=1)
    denom = na * nb

    c = np.ones(len(a), dtype=np.float64)
    valid = denom > EPS
    c[valid] = np.einsum("ij,ij->i", a[valid], b[valid]) / denom[valid]
    c = np.clip(c, -1.0, 1.0)

    out = np.degrees(np.arccos(c))
    out[~valid] = 0.0
    return out


def triangle_metrics(mesh: trimesh.Trimesh) -> dict:
    """
    Hu-style triangle quality / angle statistics.

    Triangle quality:
        Q = 2*sqrt(3)*S / (p*h)

    where S is area, p is half-perimeter, and h is longest edge.
    Q = 1 for an equilateral triangle.
    """
    V = np.asarray(mesh.vertices, dtype=np.float64)
    F = np.asarray(mesh.faces, dtype=np.int64)
    tri = V[F]

    v0, v1, v2 = tri[:, 0], tri[:, 1], tri[:, 2]

    l01 = np.linalg.norm(v1 - v0, axis=1)
    l12 = np.linalg.norm(v2 - v1, axis=1)
    l20 = np.linalg.norm(v0 - v2, axis=1)

    area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    half_perimeter = 0.5 * (l01 + l12 + l20)
    longest = np.maximum(np.maximum(l01, l12), l20)

    denom = half_perimeter * longest
    Q = np.zeros(len(F), dtype=np.float64)
    valid = denom > EPS
    Q[valid] = 2.0 * math.sqrt(3.0) * area[valid] / denom[valid]
    Q = np.clip(Q, 0.0, 1.0)

    a0 = _angle_between(v1 - v0, v2 - v0)
    a1 = _angle_between(v0 - v1, v2 - v1)
    a2 = _angle_between(v0 - v2, v1 - v2)
    angles = np.column_stack((a0, a1, a2))

    face_min = np.min(angles, axis=1)
    face_max = np.max(angles, axis=1)
    degenerate = area <= EPS

    return {
        "Q_min": float(np.min(Q)) if len(Q) else float("nan"),
        "Q_avg": float(np.mean(Q)) if len(Q) else float("nan"),
        "theta_min_deg": float(np.min(face_min)) if len(face_min) else float("nan"),
        "theta_min_avg_deg": float(np.mean(face_min)) if len(face_min) else float("nan"),
        "theta_max_deg": float(np.max(face_max)) if len(face_max) else float("nan"),
        "tri_lt_30_pct": (
            float(np.mean(face_min < 30.0) * 100.0)
            if len(face_min) else float("nan")
        ),
        "degenerate_faces": int(np.sum(degenerate)),
    }


def unique_edges_from_faces(faces: np.ndarray) -> np.ndarray:
    e = np.vstack((
        faces[:, [0, 1]],
        faces[:, [1, 2]],
        faces[:, [2, 0]],
    ))
    e.sort(axis=1)
    return np.unique(e, axis=0)


def v567_metric(mesh: trimesh.Trimesh) -> dict:
    """
    V567 = percentage of vertices whose topological valence is 5, 6, or 7.
    Valence = number of unique neighboring vertices.
    """
    F = np.asarray(mesh.faces, dtype=np.int64)
    nV = len(mesh.vertices)
    edges = unique_edges_from_faces(F)

    valence = np.bincount(edges.ravel(), minlength=nV)
    good = np.isin(valence, [5, 6, 7])

    return {
        "V567_pct": float(np.mean(good) * 100.0) if nV else float("nan"),
        "valence_min": int(valence.min()) if nV else 0,
        "valence_max": int(valence.max()) if nV else 0,
    }


def sample_surface_uniform(
    mesh: trimesh.Trimesh,
    count: int,
    seed: int,
) -> np.ndarray:
    """
    Area-uniform random sampling on a triangle mesh.

    1) choose triangles proportional to triangle area
    2) sample barycentric coordinates uniformly inside each chosen triangle
    """
    if count <= 0:
        raise ValueError("--cd-samples must be > 0.")

    V = np.asarray(mesh.vertices, dtype=np.float64)
    F = np.asarray(mesh.faces, dtype=np.int64)
    tri = V[F]

    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    area_sum = float(np.sum(areas))

    if area_sum <= EPS:
        raise ValueError("Mesh has zero total surface area.")

    prob = areas / area_sum
    rng = np.random.default_rng(seed)
    face_ids = rng.choice(len(F), size=count, replace=True, p=prob)
    chosen = tri[face_ids]

    # Uniform barycentric surface sampling.
    r1 = np.sqrt(rng.random(count))
    r2 = rng.random(count)

    w0 = 1.0 - r1
    w1 = r1 * (1.0 - r2)
    w2 = r1 * r2

    return (
        w0[:, None] * chosen[:, 0]
        + w1[:, None] * chosen[:, 1]
        + w2[:, None] * chosen[:, 2]
    )


def _query_tree_with_progress(
    tree: cKDTree,
    points: np.ndarray,
    label: str,
    chunk_size: int,
) -> np.ndarray:
    n = len(points)
    out = np.empty(n, dtype=np.float64)
    n_chunks = (n + chunk_size - 1) // chunk_size
    start_all = time.perf_counter()

    for chunk_id, start in enumerate(range(0, n, chunk_size), start=1):
        end = min(start + chunk_size, n)
        dist, _ = tree.query(points[start:end], k=1, workers=-1)
        out[start:end] = dist

        elapsed = time.perf_counter() - start_all
        fraction = end / n
        eta = elapsed * (1.0 - fraction) / max(fraction, EPS)

        progress(
            f"{label}: chunk {chunk_id}/{n_chunks}, "
            f"{end:,}/{n:,} ({fraction*100:.1f}%), ETA~{eta:.1f}s"
        )

    return out


def chamfer_metrics(
    input_mesh: trimesh.Trimesh,
    remesh: trimesh.Trimesh,
    cd_samples: int = 100_000,
    seed: int = 0,
    chunk_size: int = 50_000,
) -> dict:
    """
    Symmetric squared Chamfer Distance between area-uniform surface samples.

        CD_abs2 =
            0.5 * mean_a min_b ||a-b||^2
          + 0.5 * mean_b min_a ||b-a||^2

        CD = CD_abs2 / D_bbox^2

    where D_bbox is the INPUT mesh bounding-box diagonal.
    """
    t0 = time.perf_counter()

    progress(f"Sampling input surface: {cd_samples:,} points ...")
    ts = time.perf_counter()
    points_input = sample_surface_uniform(input_mesh, cd_samples, seed)
    input_sampling_s = time.perf_counter() - ts
    progress(f"Input sampling finished in {input_sampling_s:.2f}s")

    progress(f"Sampling remesh surface: {cd_samples:,} points ...")
    ts = time.perf_counter()
    points_remesh = sample_surface_uniform(remesh, cd_samples, seed + 1)
    remesh_sampling_s = time.perf_counter() - ts
    progress(f"Remesh sampling finished in {remesh_sampling_s:.2f}s")

    progress("Building KD-tree for remesh samples ...")
    ts = time.perf_counter()
    tree_remesh = cKDTree(points_remesh)
    tree_remesh_s = time.perf_counter() - ts
    progress(f"Remesh KD-tree built in {tree_remesh_s:.2f}s")

    progress("Computing input -> remesh nearest-neighbor distances ...")
    ts = time.perf_counter()
    d_input_to_remesh = _query_tree_with_progress(
        tree_remesh, points_input, "input -> remesh", chunk_size
    )
    input_to_remesh_s = time.perf_counter() - ts

    progress("Building KD-tree for input samples ...")
    ts = time.perf_counter()
    tree_input = cKDTree(points_input)
    tree_input_s = time.perf_counter() - ts
    progress(f"Input KD-tree built in {tree_input_s:.2f}s")

    progress("Computing remesh -> input nearest-neighbor distances ...")
    ts = time.perf_counter()
    d_remesh_to_input = _query_tree_with_progress(
        tree_input, points_remesh, "remesh -> input", chunk_size
    )
    remesh_to_input_s = time.perf_counter() - ts

    cd_i2r = float(np.mean(d_input_to_remesh ** 2))
    cd_r2i = float(np.mean(d_remesh_to_input ** 2))
    cd_abs2 = 0.5 * (cd_i2r + cd_r2i)

    V = np.asarray(input_mesh.vertices, dtype=np.float64)
    bbox_diag = float(np.linalg.norm(V.max(axis=0) - V.min(axis=0)))
    if bbox_diag <= EPS:
        raise ValueError("Input mesh bounding-box diagonal is zero.")

    cd = cd_abs2 / (bbox_diag * bbox_diag)

    return {
        "bbox_diag": bbox_diag,
        "CD": float(cd),
        "CD_x1e4": float(cd * 1.0e4),
        "CD_abs2": float(cd_abs2),
        "CD_input_to_remesh_abs2": cd_i2r,
        "CD_remesh_to_input_abs2": cd_r2i,
        "cd_samples_per_mesh": int(cd_samples),
        "input_sampling_time_s": float(input_sampling_s),
        "remesh_sampling_time_s": float(remesh_sampling_s),
        "input_tree_time_s": float(tree_input_s),
        "remesh_tree_time_s": float(tree_remesh_s),
        "input_to_remesh_time_s": float(input_to_remesh_s),
        "remesh_to_input_time_s": float(remesh_to_input_s),
        "cd_eval_time_s": float(time.perf_counter() - t0),
    }


def compute_all(
    input_path: str | Path,
    remesh_path: str | Path,
    cd_samples: int,
    seed: int,
    chunk_size: int,
    runtime_seconds: float | None,
    skip_cd: bool,
) -> dict:
    t0 = time.perf_counter()

    progress(f"Loading input mesh: {input_path}")
    ts = time.perf_counter()
    input_mesh = load_triangle_mesh(input_path)
    input_load_s = time.perf_counter() - ts
    progress(
        f"Input loaded: V={len(input_mesh.vertices):,}, "
        f"F={len(input_mesh.faces):,} in {input_load_s:.2f}s"
    )

    progress(f"Loading remesh: {remesh_path}")
    ts = time.perf_counter()
    remesh = load_triangle_mesh(remesh_path)
    remesh_load_s = time.perf_counter() - ts
    progress(
        f"Remesh loaded: V={len(remesh.vertices):,}, "
        f"F={len(remesh.faces):,} in {remesh_load_s:.2f}s"
    )

    out = {
        "input": str(input_path),
        "remesh": str(remesh_path),
        "num_vertices": int(len(remesh.vertices)),
        "num_faces": int(len(remesh.faces)),
        "input_load_time_s": float(input_load_s),
        "remesh_load_time_s": float(remesh_load_s),
    }

    progress("Computing triangle quality and angle statistics ...")
    ts = time.perf_counter()
    out.update(triangle_metrics(remesh))
    out["triangle_metrics_time_s"] = float(time.perf_counter() - ts)
    progress(f"Triangle metrics finished in {out['triangle_metrics_time_s']:.2f}s")

    progress("Computing V567 connectivity regularity ...")
    ts = time.perf_counter()
    out.update(v567_metric(remesh))
    out["v567_time_s"] = float(time.perf_counter() - ts)
    progress(f"V567 finished in {out['v567_time_s']:.2f}s")

    if skip_cd:
        progress("Skipping CD (--skip-cd).")
    else:
        progress(
            f"Starting CD: samples={cd_samples:,}/mesh, "
            f"chunk_size={chunk_size:,}"
        )
        out.update(chamfer_metrics(
            input_mesh=input_mesh,
            remesh=remesh,
            cd_samples=cd_samples,
            seed=seed,
            chunk_size=chunk_size,
        ))

    out["runtime_s"] = None if runtime_seconds is None else float(runtime_seconds)
    out["metric_eval_time_s"] = float(time.perf_counter() - t0)
    progress(f"All requested metrics finished in {out['metric_eval_time_s']:.2f}s")

    return out


def _fmt(x, digits=6):
    if x is None:
        return "-"
    if isinstance(x, float):
        if math.isnan(x):
            return "nan"
        return f"{x:.{digits}f}"
    return str(x)


def print_metrics_table(result: dict) -> None:
    rows = [
        ("#V", result.get("num_vertices")),
        ("Q_min", result.get("Q_min")),
        ("theta_min (deg)", result.get("theta_min_deg")),
        ("theta_max (deg)", result.get("theta_max_deg")),
        ("CD", result.get("CD")),
        ("CD x 1e4", result.get("CD_x1e4")),
        ("theta < 30deg (%)", result.get("tri_lt_30_pct")),
        ("V567 (%)", result.get("V567_pct")),
        ("Time (s)", result.get("runtime_s")),
    ]

    width = max(len(k) for k, _ in rows)

    print("\nRemeshing metrics")
    print("-" * (width + 26))
    for key, value in rows:
        print(f"{key:<{width}} : {_fmt(value)}")

    print("\nExtra diagnostics")
    print("-" * (width + 26))
    extras = [
        ("#F", result.get("num_faces")),
        ("Q_avg", result.get("Q_avg")),
        ("avg face min angle", result.get("theta_min_avg_deg")),
        ("degenerate faces", result.get("degenerate_faces")),
        ("CD raw squared", result.get("CD_abs2")),
        ("metric eval time (s)", result.get("metric_eval_time_s")),
    ]
    for key, value in extras:
        print(f"{key:<{width}} : {_fmt(value)}")


def write_json(result: dict, path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(result, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


def write_csv(result: dict, path: str | Path, append: bool = False) -> None:
    main_keys = [
        "input",
        "remesh",
        "num_vertices",
        "Q_min",
        "Q_avg",
        "theta_min_deg",
        "theta_min_avg_deg",
        "theta_max_deg",
        "tri_lt_30_pct",
        "V567_pct",
        "CD",
        "CD_x1e4",
        "CD_abs2",
        "cd_samples_per_mesh",
        "runtime_s",
        "num_faces",
        "degenerate_faces",
        "metric_eval_time_s",
    ]

    row = {k: result.get(k) for k in main_keys}
    path = Path(path)
    file_exists = path.exists() and path.stat().st_size > 0
    mode = "a" if append else "w"

    with open(path, mode, newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=main_keys)
        if not append or not file_exists:
            writer.writeheader()
        writer.writerow(row)


def parse_args():
    p = argparse.ArgumentParser(
        description="Compute Hu-style remeshing quality metrics plus Chamfer Distance."
    )
    p.add_argument("input_mesh", help="Original/input triangular mesh.")
    p.add_argument("remesh", help="Remeshed triangular mesh.")

    p.add_argument(
        "--cd-samples",
        type=int,
        default=100_000,
        help="Area-uniform surface samples PER MESH for CD (default: 100000).",
    )
    p.add_argument(
        "--chunk-size",
        type=int,
        default=50_000,
        help="KD-tree query chunk size (default: 50000).",
    )
    p.add_argument("--seed", type=int, default=0)

    p.add_argument(
        "--runtime-seconds",
        type=float,
        default=None,
        help="Remeshing runtime to store in the Time column.",
    )
    p.add_argument(
        "--skip-cd",
        action="store_true",
        help="Skip CD. By default CD is always computed.",
    )

    p.add_argument("--json", dest="json_path", default=None)
    p.add_argument("--csv", dest="csv_path", default=None)
    p.add_argument(
        "--append-csv",
        action="store_true",
        help="Append one experiment row to --csv instead of overwriting it.",
    )
    return p.parse_args()


def main():
    args = parse_args()

    result = compute_all(
        input_path=args.input_mesh,
        remesh_path=args.remesh,
        cd_samples=args.cd_samples,
        seed=args.seed,
        chunk_size=args.chunk_size,
        runtime_seconds=args.runtime_seconds,
        skip_cd=args.skip_cd,
    )

    print_metrics_table(result)

    if args.json_path:
        write_json(result, args.json_path)
        print(f"\nJSON saved to: {args.json_path}")

    if args.csv_path:
        write_csv(result, args.csv_path, append=args.append_csv)
        action = "appended to" if args.append_csv else "saved to"
        print(f"CSV {action}: {args.csv_path}")


if __name__ == "__main__":
    main()
