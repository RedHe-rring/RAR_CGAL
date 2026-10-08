"""Analytic fixtures for metrics, surface distance, region weights and volume gating."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
import numpy as np
import trimesh
from analyze import closest, membership
from mesh_metrics import quality, locator

class AnalysisTests(unittest.TestCase):
    def test_equilateral_and_surface_distance(self):
        mesh=trimesh.Trimesh(vertices=[[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]],
                             faces=[[0,1,2],[0,3,1],[0,2,3],[1,3,2]],process=False)
        q=quality(mesh)
        self.assertAlmostEqual(q['triQ05'],1)
        self.assertAlmostEqual(q['AR95'],1)
        self.assertAlmostEqual(q['minAngP05'],60)
        plane=trimesh.Trimesh(vertices=[[0,0,0],[1,0,0],[0,1,0]],faces=[[0,1,2]],process=False)
        loc,poly=locator(plane)
        distance,ids,positions=closest(loc,np.array([[.2,.3,0],[.2,.3,2.]]))
        np.testing.assert_allclose(distance,[0,2])  # closest point is in face interior
        weights=membership(plane,np.eye(3),ids,positions)
        np.testing.assert_allclose(weights,[[.5,.2,.3],[.5,.2,.3]])

    def test_closed_is_not_sufficient_for_volume(self):
        exe=Path(__file__).resolve().parents[2]/'build/Release/rar_validate_volume.exe'
        if not exe.exists(): self.skipTest('Build rar_validate_volume for geometric volume checks')
        with tempfile.TemporaryDirectory(prefix='radius-volume-') as directory:
            path=Path(directory)/'fixture.ply'; box=trimesh.creation.box();box.export(path)
            checked=json.loads(subprocess.check_output([str(exe),str(path)],text=True))
            self.assertTrue(checked['bounds_volume']);self.assertFalse(checked['self_intersects'])
            shifted=box.copy();shifted.apply_translation([.5,.5,.5])
            trimesh.util.concatenate([box,shifted]).export(path)
            checked=json.loads(subprocess.check_output([str(exe),str(path)],text=True))
            self.assertTrue(checked['closed']);self.assertTrue(checked['self_intersects'])
            self.assertFalse(checked['bounds_volume'])

if __name__=='__main__':unittest.main()
