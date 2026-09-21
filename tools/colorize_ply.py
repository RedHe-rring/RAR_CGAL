#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Colorize one scalar vertex property in an ASCII PLY file.

The script preserves the original mesh topology and vertex properties, and
writes/overwrites standard RGB vertex properties using a blue -> cyan ->
green -> yellow -> red heatmap.

By default, the scalar range is the full [min, max] range of the selected
property (0th to 100th percentile, i.e. no percentile clipping).

If the output path is omitted, a parameter-aware filename is generated from
the input field filename plus the visualization settings.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Sequence, Tuple


STANDARD_RGB_PROPERTIES = {"red", "green", "blue"}


def clamp01(t: float) -> float:
    return max(0.0, min(1.0, t))


def heat_color(t: float) -> Tuple[int, int, int]:
    """Blue -> cyan -> green -> yellow -> red."""
    t = clamp01(t)

    stops = [
        (0.00, (0, 0, 255)),
        (0.25, (0, 255, 255)),
        (0.50, (0, 255, 0)),
        (0.75, (255, 255, 0)),
        (1.00, (255, 0, 0)),
    ]

    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t <= t1:
            u = 0.0 if t1 == t0 else (t - t0) / (t1 - t0)
            return tuple(
                int(round(c0[k] + u * (c1[k] - c0[k])))
                for k in range(3)
            )

    return stops[-1][1]


def filename_number(value: float) -> str:
    """Encode a number in the same filename-friendly style as the C++ tool."""
    s = format(value, ".10g").replace(".", "p")

    out = []
    for ch in s:
        if ch == "-":
            out.append("m")
        elif ch == "+":
            out.append("p")
        else:
            out.append(ch)
    return "".join(out)


def make_auto_output_path(
    input_path: Path,
    property_name: str,
    invert: bool,
    lo: float,
    hi: float,
) -> Path:
    safe_property = "".join(
        ch if ch.isalnum() or ch in {"_", "-"} else "-"
        for ch in property_name
    )

    name = (
        f"{input_path.stem}"
        f"__color-{safe_property}"
        f"__invert-{'on' if invert else 'off'}"
        f"__vmin-{filename_number(lo)}"
        f"__vmax-{filename_number(hi)}"
        ".ply"
    )

    return input_path.with_name(name)


def parse_ascii_ply(
    path: Path,
) -> Tuple[List[str], List[str], int, List[str], List[List[str]]]:
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if not lines or lines[0].strip() != "ply":
        raise ValueError(f"Not a PLY file: {path}")

    try:
        end_header = next(
            i for i, line in enumerate(lines)
            if line.strip() == "end_header"
        )
    except StopIteration as exc:
        raise ValueError(f"Missing end_header: {path}") from exc

    header = lines[: end_header + 1]
    body = lines[end_header + 1 :]

    if not any(
        line.strip() == "format ascii 1.0"
        for line in header
    ):
        raise ValueError("Only ASCII PLY format is supported.")

    vertex_count = None
    vertex_properties: List[str] = []
    in_vertex_element = False

    for line in header:
        parts = line.strip().split()
        if not parts:
            continue

        if parts[0] == "element":
            in_vertex_element = (
                len(parts) >= 3 and parts[1] == "vertex"
            )
            if in_vertex_element:
                vertex_count = int(parts[2])
            continue

        if (
            in_vertex_element
            and parts[0] == "property"
            and len(parts) >= 3
        ):
            if parts[1] == "list":
                raise ValueError(
                    "List properties are not supported "
                    "inside the vertex element."
                )
            vertex_properties.append(parts[-1])

    if vertex_count is None:
        raise ValueError("PLY has no vertex element.")

    if len(body) < vertex_count:
        raise ValueError(
            "PLY body is shorter than the declared vertex count."
        )

    vertex_rows = [
        body[i].split()
        for i in range(vertex_count)
    ]

    tail = body[vertex_count:]

    expected_columns = len(vertex_properties)
    for i, row in enumerate(vertex_rows):
        if len(row) != expected_columns:
            raise ValueError(
                f"Vertex {i} has {len(row)} columns, "
                f"expected {expected_columns}."
            )

    return (
        header,
        tail,
        vertex_count,
        vertex_properties,
        vertex_rows,
    )


def scalar_range(
    values: Sequence[float],
    explicit_min: float | None,
    explicit_max: float | None,
) -> Tuple[float, float]:
    if not values:
        raise ValueError("No scalar values found.")

    lo = min(values) if explicit_min is None else explicit_min
    hi = max(values) if explicit_max is None else explicit_max

    if hi < lo:
        raise ValueError(
            f"Invalid scalar range: min={lo}, max={hi}"
        )

    return lo, hi


def rewrite_vertex_header_with_rgb(
    header: Sequence[str],
    vertex_properties: Sequence[str],
) -> List[str]:
    existing_rgb = STANDARD_RGB_PROPERTIES.intersection(
        vertex_properties
    )

    new_header: List[str] = []
    in_vertex_element = False
    rgb_inserted = False

    for line in header:
        parts = line.strip().split()

        if parts and parts[0] == "element":
            if in_vertex_element and not rgb_inserted:
                new_header.extend(
                    [
                        "property uchar red",
                        "property uchar green",
                        "property uchar blue",
                    ]
                )
                rgb_inserted = True

            in_vertex_element = (
                len(parts) >= 3 and parts[1] == "vertex"
            )
            new_header.append(line)
            continue

        if in_vertex_element and parts and parts[0] == "property":
            prop_name = parts[-1]
            if prop_name in existing_rgb:
                continue

        if line.strip() == "end_header" and not rgb_inserted:
            new_header.extend(
                [
                    "property uchar red",
                    "property uchar green",
                    "property uchar blue",
                ]
            )
            rgb_inserted = True

        new_header.append(line)

    return new_header


def colorize_ply(
    input_path: Path,
    output_path: Path | None,
    property_name: str,
    invert: bool,
    explicit_min: float | None,
    explicit_max: float | None,
) -> Path:
    (
        header,
        tail,
        _vertex_count,
        vertex_properties,
        vertex_rows,
    ) = parse_ascii_ply(input_path)

    if property_name not in vertex_properties:
        raise ValueError(
            f"Vertex property '{property_name}' not found. "
            f"Available properties: {', '.join(vertex_properties)}"
        )

    scalar_index = vertex_properties.index(property_name)

    try:
        values = [
            float(row[scalar_index])
            for row in vertex_rows
        ]
    except ValueError as exc:
        raise ValueError(
            f"Property '{property_name}' contains "
            "non-numeric values."
        ) from exc

    lo, hi = scalar_range(
        values,
        explicit_min,
        explicit_max,
    )

    if output_path is None:
        output_path = make_auto_output_path(
            input_path=input_path,
            property_name=property_name,
            invert=invert,
            lo=lo,
            hi=hi,
        )

    new_header = rewrite_vertex_header_with_rgb(
        header,
        vertex_properties,
    )

    new_vertex_rows: List[str] = []

    for value, row in zip(values, vertex_rows):
        kept = [
            token
            for i, token in enumerate(row)
            if vertex_properties[i]
            not in STANDARD_RGB_PROPERTIES
        ]

        if hi == lo:
            t = 0.5
        else:
            t = (value - lo) / (hi - lo)

        t = clamp01(t)
        if invert:
            t = 1.0 - t

        r, g, b = heat_color(t)
        kept.extend([str(r), str(g), str(b)])
        new_vertex_rows.append(" ".join(kept))

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_lines = (
        new_header
        + new_vertex_rows
        + tail
    )

    output_path.write_text(
        "\n".join(output_lines) + "\n",
        encoding="utf-8",
    )

    print(f"Input:    {input_path}")
    print(f"Output:   {output_path}")
    print(f"Property: {property_name}")
    print(f"Range:    [{lo}, {hi}]")
    print(f"Invert:   {invert}")

    return output_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Map one ASCII PLY vertex scalar property "
            "to RGB heatmap colors."
        )
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Input ASCII PLY file.",
    )
    parser.add_argument(
        "output",
        type=Path,
        nargs="?",
        default=None,
        help=(
            "Optional output PLY. If omitted, a parameter-aware "
            "filename is generated beside the input file."
        ),
    )
    parser.add_argument(
        "--property",
        required=True,
        dest="property_name",
        help=(
            "Vertex scalar property to colorize, "
            "e.g. curvature or target_length."
        ),
    )
    parser.add_argument(
        "--invert",
        action="store_true",
        help="Reverse the heatmap direction.",
    )
    parser.add_argument(
        "--min",
        dest="explicit_min",
        type=float,
        default=None,
        help=(
            "Optional explicit lower bound. "
            "Default: actual minimum."
        ),
    )
    parser.add_argument(
        "--max",
        dest="explicit_max",
        type=float,
        default=None,
        help=(
            "Optional explicit upper bound. "
            "Default: actual maximum."
        ),
    )

    return parser


def main() -> None:
    args = build_parser().parse_args()

    colorize_ply(
        input_path=args.input,
        output_path=args.output,
        property_name=args.property_name,
        invert=args.invert,
        explicit_min=args.explicit_min,
        explicit_max=args.explicit_max,
    )


if __name__ == "__main__":
    main()
