#!/usr/bin/env python3
"""C1: error-aware curvature-ordered sizing quantization for triangle meshes.

Research prototype; neither the local sizing proxy nor the CGAL remeshing
backend provides a certified Hausdorff error bound.
"""
from __future__ import annotations

import argparse
import csv
import heapq
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.spatial import cKDTree


@dataclass
class PartitionResult:
    cuts: tuple[int, ...]
    group_per_bin: np.ndarray
    group_h: np.ndarray
    group_upper: np.ndarray
    group_area: np.ndarray
    adjacency: np.ndarray
    cost: float


def vertex_areas(v: np.ndarray, f: np.ndarray):
    cross = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
    face_area = .5 * np.linalg.norm(cross, axis=1)
    area = np.zeros(len(v), dtype=np.float64)
    for j in range(3):
        np.add.at(area, f[:, j], face_area / 3)
    return area, face_area


def normal_curvature(v: np.ndarray, f: np.ndarray, edges: np.ndarray,
                     smooth_steps: int, smooth_rate: float):
    """Angle proxy for binning and normal-angle / edge-length for sizing.

    CSF-inspired, not an exact reproduction of Lv et al. or principal curvature.
    """
    p0, p1, p2 = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    cross = np.cross(p1-p0, p2-p0)
    normals = np.zeros_like(v)
    for j in range(3):
        np.add.at(normals, f[:, j], cross)
    norm = np.linalg.norm(normals, axis=1)
    valid = norm > 1e-30
    normals[valid] /= norm[valid, None]
    u, w = edges[:, 0], edges[:, 1]
    length = np.linalg.norm(v[u]-v[w], axis=1)
    cos = np.einsum('ij,ij->i', normals[u], normals[w])
    angles = np.arccos(np.clip(cos, -1., 1.))
    angles[~(valid[u] & valid[w])] = 0.
    s = np.zeros(len(v))
    kappa = np.zeros(len(v))
    np.maximum.at(s, u, angles)
    np.maximum.at(s, w, angles)
    edge_k = np.divide(angles, length, out=np.zeros_like(angles),
                       where=length > 1e-30)
    np.maximum.at(kappa, u, edge_k)
    np.maximum.at(kappa, w, edge_k)

    den = np.maximum(np.linalg.norm(cross, axis=1), 1e-30)
    cot0 = np.einsum('ij,ij->i', p1-p0, p2-p0) / den
    cot1 = np.einsum('ij,ij->i', p2-p1, p0-p1) / den
    cot2 = np.einsum('ij,ij->i', p0-p2, p1-p2) / den
    aa = np.concatenate([f[:, 1], f[:, 2], f[:, 0]])
    bb = np.concatenate([f[:, 2], f[:, 0], f[:, 1]])
    ww = np.maximum(.5*np.concatenate([cot0, cot1, cot2]), 0.)
    n = len(v)
    adj = coo_matrix((np.r_[ww, ww], (np.r_[aa, bb], np.r_[bb, aa])),
                     shape=(n,n)).tocsr()
    weak = coo_matrix((np.full(2*len(edges), 1e-8),
                       (np.r_[u,w], np.r_[w,u])), shape=(n,n)).tocsr()
    adj += weak
    degree = np.asarray(adj.sum(axis=1)).ravel()

    def smooth(x):
        y = x.copy()
        vmax, vmin = int(np.argmax(x)), int(np.argmin(x))
        for _ in range(smooth_steps):
            avg = np.divide(adj @ y, degree, out=y.copy(), where=degree > 0)
            y = (1-smooth_rate)*y + smooth_rate*avg
            if x[vmax] > x[vmin] + 1e-15:
                y[vmax], y[vmin] = x[vmax], x[vmin]
        return y

    return s, smooth(s), kappa, smooth(kappa)


def weighted_quantile(x, w, q):
    order = np.argsort(x)
    x, w = x[order], w[order]
    good = w > 0
    x, w = x[good], w[good]
    if not len(x):
        return 0.
    return float(x[np.searchsorted(np.cumsum(w), q*np.sum(w), side='left')])


def error_upper(kappa, epsilon, hmin, hmax):
    """Heuristic radius-guarded proxy; NOT a certified error bound."""
    kappa = np.asarray(kappa)
    upper = np.full(len(kappa), hmax)
    good = np.isfinite(kappa) & (kappa > 1e-12)
    formula = good & (epsilon*kappa <= .183)
    upper[formula] = np.sqrt(np.maximum(
        0., 6*epsilon/kappa[formula] - 3*epsilon**2))
    radius = good & ~formula
    upper[radius] = 1/kappa[radius]
    conflicts = int(np.count_nonzero(upper < hmin*(1-1e-9)))
    return np.clip(upper, hmin, hmax), formula, int(radius.sum()), conflicts


def group_bins(signal, requested):
    cuts = np.unique(np.quantile(signal, np.linspace(0,1,requested+1)[1:-1]))
    ids = np.searchsorted(cuts, signal, side='right').astype(np.int32)
    used = np.unique(ids)
    if len(used) < int(ids.max())+1:
        ids = np.searchsorted(used, ids).astype(np.int32)
    return ids, cuts


def shortest_upper(upper, links, beta):
    """Exact componentwise greatest h <= upper with adjacent ratio <= beta."""
    n = len(upper)
    graph = [[] for _ in range(n)]
    for a,b in links:
        a,b = int(a),int(b)
        if a == b:
            continue
        graph[a].append(b)
        graph[b].append(a)
    y = np.log(upper.copy())
    queue = [(float(x),i) for i,x in enumerate(y)]
    heapq.heapify(queue)
    step = math.log(beta)
    while queue:
        d,i = heapq.heappop(queue)
        if d > y[i] + 1e-13:
            continue
        for j in graph[i]:
            nd = d+step
            if nd < y[j]-1e-13:
                y[j] = nd
                heapq.heappush(queue,(nd,j))
    return np.exp(y)


class C1Optimizer:
    def __init__(self, base_id, areas, H, edges, beta):
        self.base_id = base_id
        self.beta = beta
        self.B = int(base_id.max())+1
        self.bin_area = np.bincount(base_id, weights=areas, minlength=self.B)
        self.bin_upper = np.full(self.B,np.inf)
        np.minimum.at(self.bin_upper,base_id,H)
        pairs = base_id[edges]
        pairs = pairs[pairs[:,0] != pairs[:,1]]
        self.base_links = (np.unique(np.sort(pairs,axis=1),axis=0) if len(pairs)
                           else np.empty((0,2),dtype=np.int32))
        self.cache = {}

    def evaluate(self,cuts):
        cuts = tuple(sorted(int(i) for i in cuts))
        if cuts in self.cache:
            return self.cache[cuts]
        starts = np.array((0,)+cuts,dtype=np.int32)
        groups = np.searchsorted(cuts,np.arange(self.B),side='right').astype(np.int32)
        area = np.add.reduceat(self.bin_area,starts)
        upper = np.minimum.reduceat(self.bin_upper,starts)
        pairs = groups[self.base_links]
        pairs = pairs[pairs[:,0]!=pairs[:,1]]
        links = (np.unique(np.sort(pairs,axis=1),axis=0) if len(pairs)
                 else np.empty((0,2),dtype=np.int32))
        h = shortest_upper(upper,links,self.beta)
        cost = float(np.sum(area/(h*h)))
        r = PartitionResult(cuts,groups,h,upper,area,links,cost)
        self.cache[cuts] = r
        return r

    def solve(self,max_groups=12,gain_fraction=.02,refine_rounds=8,fixed_groups=None):
        if fixed_groups:
            boundary = np.linspace(0,self.B,min(fixed_groups,self.B)+1)
            cuts = tuple(sorted(set(int(round(x)) for x in boundary[1:-1])-{0,self.B}))
            result = self.evaluate(cuts)
            return result,[dict(K=len(cuts)+1,cost=result.cost,operation='fixed')]
        result = self.evaluate(())
        history = [dict(K=1,cost=result.cost,operation='initial')]
        while len(result.cuts)+1 < min(max_groups,self.B):
            current = set(result.cuts)
            best = None
            for s in range(1,self.B):
                if s in current:
                    continue
                trial = self.evaluate(tuple(sorted(current|{s})))
                if best is None or trial.cost < best.cost:
                    best = trial
            if best is None:
                break
            gain = (result.cost-best.cost)/max(result.cost,1e-30)
            if gain < gain_fraction-1e-12:
                break
            result = best
            history.append(dict(K=len(result.cuts)+1,cost=result.cost,
                                relative_gain=float(gain),operation='split'))
        for _ in range(refine_rounds):
            best = result
            for i,cut in enumerate(result.cuts):
                for delta in (-1,1):
                    trial_cuts = list(result.cuts)
                    trial_cuts[i] = cut+delta
                    if not 0 < trial_cuts[i] < self.B or len(set(trial_cuts))!=len(trial_cuts):
                        continue
                    trial = self.evaluate(tuple(sorted(trial_cuts)))
                    if trial.cost < best.cost-1e-9*max(1.,best.cost):
                        best = trial
            if best is result:
                break
            result = best
            history.append(dict(K=len(result.cuts)+1,cost=result.cost,
                                operation='boundary-shift'))
        return result,history


def mesh_quality(mesh):
    if not len(mesh.faces):
        return {}
    tri = mesh.vertices[mesh.faces]
    edge = np.column_stack([np.linalg.norm(tri[:,1]-tri[:,2],axis=1),
                            np.linalg.norm(tri[:,2]-tri[:,0],axis=1),
                            np.linalg.norm(tri[:,0]-tri[:,1],axis=1)])
    a,b,c = edge.T
    ang = np.rad2deg(np.arccos(np.clip(np.column_stack([
        (b*b+c*c-a*a)/np.maximum(2*b*c,1e-30),
        (c*c+a*a-b*b)/np.maximum(2*c*a,1e-30),
        (a*a+b*b-c*c)/np.maximum(2*a*b,1e-30)]),-1,1)))
    area = .5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    per2 = edge.sum(axis=1)/2
    # Q=1 for an equilateral triangle; Lv et al.'s area/perimeter/longest-edge metric.
    quality = 2*np.sqrt(3)*area/np.maximum(per2*edge.max(axis=1),1e-30)
    return dict(vertices=int(len(mesh.vertices)),faces=int(len(mesh.faces)),
                min_angle_p05_deg=float(np.percentile(ang.min(axis=1),5)),
                min_angle_min_deg=float(ang.min()),
                triangle_quality_mean=float(quality.mean()),
                triangle_quality_p05=float(np.percentile(quality,5)))


def sampled_chamfer(a,b,n,seed):
    np.random.seed(seed)
    p,_ = trimesh.sample.sample_surface(a,n)
    q,_ = trimesh.sample.sample_surface(b,n)
    da = cKDTree(p).query(q,workers=-1)[0]
    db = cKDTree(q).query(p,workers=-1)[0]
    return dict(sample_count_each=n,chamfer_squared=float((da**2).mean()+(db**2).mean()),
                sampled_distance_p95=float(np.percentile(np.r_[da,db],95)),
                note="Sampled point-cloud estimates; NOT a Hausdorff guarantee.")


def write_heatmap(mesh,values,path,discrete=False):
    from matplotlib import colormaps
    if discrete:
        rgba = colormaps['tab20']((values.astype(int)%20)/20.)
    else:
        lo,hi = np.percentile(values,[2,98])
        rgba = colormaps['viridis'](np.clip((values-lo)/max(hi-lo,1e-15),0,1))
    colored = trimesh.Trimesh(vertices=mesh.vertices,faces=mesh.faces,process=False,
                              vertex_colors=(rgba*255).astype(np.uint8))
    colored.export(path)


def run(args):
    outdir = Path(args.output)
    outdir.mkdir(parents=True,exist_ok=True)
    mesh = trimesh.load(args.input,force='mesh',process=False)
    if not isinstance(mesh,trimesh.Trimesh) or not len(mesh.faces):
        raise ValueError("Input must be a nonempty triangular mesh")
    v = np.asarray(mesh.vertices,dtype=np.float64)
    f = np.asarray(mesh.faces,dtype=np.int64)
    if f.ndim!=2 or f.shape[1]!=3 or not np.isfinite(v).all():
        raise ValueError("Invalid triangular mesh")
    if np.any(np.bincount(f.ravel(),minlength=len(v))==0):
        raise ValueError("Isolated vertices: clean mesh first")
    edges = np.asarray(mesh.edges_unique,dtype=np.int64)
    area,face_area = vertex_areas(v,f)
    if np.any(face_area<1e-22):
        raise ValueError("Degenerate triangles: clean mesh first")
    if not 0 < args.min_edge < args.max_edge:
        raise ValueError("Require 0 < min-edge < max-edge")
    if not args.beta>1 or args.bins<1 or args.max_groups<1:
        raise ValueError("Require beta>1, bins>=1 and max-groups>=1")
    s,ss,k,ks = normal_curvature(v,f,edges,args.smooth_steps,args.smooth_rate)
    if args.epsilon=='auto':
        good = k>1e-12
        k90 = weighted_quantile(k[good],area[good],.9)
        eps = .183/k90 if k90>1e-12 else .001*max(float(np.linalg.norm(np.ptp(v,axis=0))),1e-3)
    else:
        eps = float(args.epsilon)
    if not np.isfinite(eps) or eps<=0:
        raise ValueError("epsilon must be finite positive")
    H,valid_formula,radius_count,conflicts = error_upper(k,eps,args.min_edge,args.max_edge)
    signal = ks if args.group_signal=='curvature' else ss
    bin_id,_ = group_bins(signal,args.bins)
    opt = C1Optimizer(bin_id,area,H,edges,args.beta)
    result,history = opt.solve(args.max_groups,args.min_gain,args.refine_rounds,
                               fixed_groups=5 if args.mode=='fixed5' else None)
    group = result.group_per_bin[bin_id]
    h = result.group_h[group]
    if np.max(h-H)>1e-7:
        raise AssertionError("Group target exceeds H proxy")
    if len(result.adjacency):
        x,y = result.adjacency.T
        ratio = np.max(np.maximum(result.group_h[x]/result.group_h[y],
                                  result.group_h[y]/result.group_h[x]))
        if ratio>args.beta+1e-9:
            raise AssertionError("Group gradation violated")
    with (outdir/'field.csv').open('w',newline='',encoding='utf-8') as fp:
        writer=csv.writer(fp)
        writer.writerow(['vertex_id','x','y','z','target_length','group_id','s_raw',
                         's_smooth','kappa_proxy','H_upper','kappa_smooth','group_signal'])
        for i,p in enumerate(v):
            writer.writerow([i,*[format(t,'.17g') for t in p],format(h[i],'.17g'),
                             int(group[i]),format(s[i],'.17g'),format(ss[i],'.17g'),
                             format(k[i],'.17g'),format(H[i],'.17g'),
                             format(ks[i],'.17g'),format(signal[i],'.17g')])
    with (outdir/'groups.csv').open('w',newline='',encoding='utf-8') as fp:
        writer=csv.writer(fp)
        writer.writerow(['group_id','curvature_min','curvature_max','area',
                         'H_upper','h_target','vertex_count'])
        for i in range(len(result.group_h)):
            mask = group==i
            writer.writerow([i,float(signal[mask].min()),float(signal[mask].max()),
                             float(result.group_area[i]),float(result.group_upper[i]),
                             float(result.group_h[i]),int(mask.sum())])
    if not args.no_heatmap:
        write_heatmap(mesh,h,outdir/'target_length_heat.ply')
        write_heatmap(mesh,group,outdir/'groups_heat.ply',True)
    summary=dict(model=str(args.input),mode=args.mode,epsilon=eps,
                 epsilon_source=args.epsilon,beta=args.beta,K=len(result.group_h),
                 cuts_base_bins=list(result.cuts),base_bins=opt.B,
                 total_area=float(area.sum()),
                 estimated_faces=float((4/np.sqrt(3))*result.cost),
                 target_min=float(h.min()),target_max=float(h.max()),
                 target_ratio=float(h.max()/h.min()),
                 curvature_proxy_max=float(k.max()),
                 radius_fallback_vertices=radius_count,
                 min_edge_conflict_vertices=conflicts,
                 smooth_iterations=args.smooth_steps,history=history,
                 field_csv=str(outdir/'field.csv'),input_quality=mesh_quality(mesh),
                 limitations=['H is a geometric proxy, NOT certified Hausdorff bound',
                              'Target is constant per curvature bin',
                              'beta is not a triangle quality guarantee'])
    if args.remesh_exe:
        exe=Path(args.remesh_exe)
        if not exe.is_file():
            raise FileNotFoundError(f"CGAL executable missing: {exe}")
        outmesh=outdir/'c1_remeshed.ply'
        cmd=[str(exe.resolve()),str(Path(args.input).resolve()),str(outmesh.resolve()),
             '--field','c1-external','--sizing-csv',str((outdir/'field.csv').resolve()),
             '--iterations',str(args.iterations),'--relax-steps',str(args.relax_steps),
             '--no-export-field']
        if args.no_project: cmd.append('--no-project')
        if args.preserve_features:
            cmd+=['--preserve-features','--feature-angle',str(args.feature_angle)]
        res=subprocess.run(cmd,text=True,capture_output=True)
        (outdir/'remesher.log').write_text(res.stdout+'\n'+res.stderr,encoding='utf-8')
        summary['remesh_returncode']=res.returncode
        summary['remesh_command']=cmd
        if res.returncode==0 and outmesh.is_file():
            output=trimesh.load(outmesh,force='mesh',process=False)
            summary['output_quality']=mesh_quality(output)
            if args.samples>0:
                summary['approx_geometry']=sampled_chamfer(mesh,output,args.samples,args.seed)
        else:
            summary['remesh_error']='See remesher.log'
    (outdir/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({key:summary[key] for key in
                      ('mode','K','epsilon','estimated_faces','target_min','target_max',
                       'min_edge_conflict_vertices')},indent=2))
    return summary


def cli():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--epsilon',default='auto')
    p.add_argument('--min-edge',type=float,default=.001)
    p.add_argument('--max-edge',type=float,default=.5)
    p.add_argument('--beta',type=float,default=1.4)
    p.add_argument('--bins',type=int,default=64)
    p.add_argument('--group-signal',choices=['curvature','normal-angle'],default='curvature')
    p.add_argument('--max-groups',type=int,default=12)
    p.add_argument('--min-gain',type=float,default=.02)
    p.add_argument('--refine-rounds',type=int,default=8)
    p.add_argument('--smooth-steps',type=int,default=8)
    p.add_argument('--smooth-rate',type=float,default=.5)
    p.add_argument('--mode',choices=['c1','fixed5'],default='c1')
    p.add_argument('--remesh-exe',default='')
    p.add_argument('--iterations',type=int,default=5)
    p.add_argument('--relax-steps',type=int,default=3)
    p.add_argument('--no-project',action='store_true')
    p.add_argument('--preserve-features',action='store_true')
    p.add_argument('--feature-angle',type=float,default=50)
    p.add_argument('--samples',type=int,default=5000)
    p.add_argument('--seed',type=int,default=12345)
    p.add_argument('--no-heatmap',action='store_true')
    a=p.parse_args()
    if a.smooth_steps<0 or not 0<a.smooth_rate<=1:
        p.error('smooth-steps >= 0 and 0 < smooth-rate <= 1 required')
    return a


if __name__=='__main__':
    run(cli())
