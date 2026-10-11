"""Deterministic integration smoke tests for diagnose_sizing.py."""
import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import vtk

ROOT = Path(__file__).resolve().parent


def grid(n, z=0.0):
    coords = np.linspace(0, 1, n + 1)
    verts = np.array([[x, y, z] for y in coords for x in coords], float)
    faces = []
    for j in range(n):
        for i in range(n):
            a = j * (n + 1) + i
            b = a + 1
            c = a + (n + 1)
            d = c + 1
            faces.extend([(a, b, c), (b, d, c)])
    return verts, np.asarray(faces, dtype=int)


def ply(path, verts, faces):
    with open(path, 'w', encoding='utf-8') as f:
        f.write('ply\nformat ascii 1.0\n')
        f.write(f'element vertex {len(verts)}\nproperty double x\nproperty double y\nproperty double z\n')
        f.write(f'element face {len(faces)}\nproperty list uchar int vertex_indices\nend_header\n')
        for v in verts: f.write(' '.join(map(str, v))+'\n')
        for q in faces: f.write('3 '+' '.join(map(str, q))+'\n')


def run_test(offset):
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        vi, fi = grid(2)
        vo, fo = grid(4, offset)
        ply(base/'field.ply', vi, fi)
        ply(base/'output.ply', vo, fo)
        with open(base/'field.csv', 'w', newline='') as f:
            w=csv.writer(f)
            w.writerow(['vertex_id','x','y','z','curvature','raw_target_length','target_length'])
            for i,p in enumerate(vi):
                w.writerow([i,*p,0.0,2.0,1.0])
        subprocess.run([sys.executable, str(ROOT/'diagnose_sizing.py'),
                        '--field-ply',str(base/'field.ply'), '--field-csv',str(base/'field.csv'),
                        '--output-mesh',str(base/'output.ply'), '--out-dir',str(base/'diag'),
                        '--epsilon','0.2','--max-edges','0','--max-vertices','0',
                        '--max-faces','0'], check=True, capture_output=True, text=True)
        diag=json.loads((base/'diag'/'summary.json').read_text())
        assert diag['size']['ratio']['max'] < 0.36, diag['size']['ratio']
        assert abs(diag['size']['chen_h_over_h0']['median']-.5) < 1e-12
        val=diag['geometric_error_sampled_NOT_Hausdorff']['two_way_sampled_max']
        assert abs(val - offset) < 1e-8, (val,offset)
        for f in ('size_mismatch_mesh.ply','size_mismatch_edge_points.ply',
                  'chen_correction_input.ply','error_output_to_input.ply',
                  'error_input_to_output.ply','edges.csv','summary.json',
                  'size_histogram.png','error_histogram.png'):
            assert (base/'diag'/f).exists(), f
        color_reader = vtk.vtkPLYReader()
        color_reader.SetFileName(str(base/'diag'/'size_mismatch_mesh.ply'))
        color_reader.Update()
        assert color_reader.GetOutput().GetNumberOfPoints() == len(vo)
        assert color_reader.GetOutput().GetNumberOfPolys() == len(fo)
        assert color_reader.GetOutput().GetPointData().GetScalars() is not None
        return diag


if __name__=='__main__':
    for z in (0., 0.1):
        result=run_test(z)
        print(f'PASS planar grid z-offset={z}: ratio_p95={result["size"]["ratio"]["p95"]:.5f}, '
              f'geometric deviation={result["geometric_error_sampled_NOT_Hausdorff"]["two_way_sampled_max"]:.6g}')
