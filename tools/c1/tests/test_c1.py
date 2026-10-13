import sys
from argparse import Namespace
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import trimesh
from c1_pipeline import C1Optimizer,error_upper,shortest_upper,run,mesh_quality


def test_graph_gradation_chain():
    h=shortest_upper(np.array([.1,.5,2]),np.array([[0,1],[1,2]]),1.4)
    np.testing.assert_allclose(h,[.1,.14,.196],atol=1e-12)


def test_graph_gradation_nonlocal():
    h=shortest_upper(np.array([.1,.5,2]),np.array([[0,1],[1,2],[0,2]]),1.4)
    np.testing.assert_allclose(h,[.1,.14,.14],atol=1e-12)


def test_split_cannot_raise_objective():
    opt=C1Optimizer(np.array([0,1,2,3]),np.ones(4),
                    np.array([.1,.2,.4,.8]),
                    np.array([[0,1],[1,2],[2,3]]),1.4)
    costs=[opt.evaluate(c).cost for c in [(),(1,),(1,2),(1,2,3)]]
    assert all(a>=b-1e-8 for a,b in zip(costs,costs[1:]))
    result,history=opt.solve(max_groups=4,gain_fraction=.02,refine_rounds=4)
    assert 1<=len(result.group_h)<=4 and history


def test_local_upper_clamped_and_reports_conflicts():
    k=np.array([0,.01,1,100.])
    h,_,radius,conflict=error_upper(k,.01,.001,.5)
    assert np.all(h>=.001) and np.all(h<=.5)
    assert radius>=1 and h[0]==.5


def test_equilateral_quality_normalized():
    tri=trimesh.Trimesh(vertices=np.array([
        [0,0,0],[1,0,0],[.5,np.sqrt(3)/2,0]
    ]),faces=np.array([[0,1,2]]),process=False)
    assert abs(mesh_quality(tri)['triangle_quality_mean']-1)<1e-12


def test_end_to_end(tmp_path):
    mesh=trimesh.creation.icosphere(subdivisions=2,radius=.5)
    path=tmp_path/'sphere.ply'
    mesh.export(path)
    args=Namespace(input=str(path),output=str(tmp_path/'out'),epsilon='.003',
        min_edge=.001,max_edge=.5,beta=1.4,bins=32,max_groups=8,
        min_gain=.02,refine_rounds=3,smooth_steps=5,smooth_rate=.5,
        group_signal='curvature',mode='c1',remesh_exe='',iterations=5,
        relax_steps=3,no_project=False,preserve_features=False,
        feature_angle=50,samples=100,seed=123,no_heatmap=True)
    summary=run(args)
    csv=tmp_path/'out'/'field.csv'
    assert csv.exists() and len(csv.read_text().splitlines())==len(mesh.vertices)+1
    assert summary['K']<=8
    assert (tmp_path/'out'/'summary.json').exists()
