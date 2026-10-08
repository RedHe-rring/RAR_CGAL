"""Fixed-seed, triangle-surface error metrics and area-weighted field diagnostics."""
import argparse
import json
import re
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import trimesh
import vtk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mesh_metrics import quality, locator

REGIONS = ['Unchanged','Cap','Recovery']
COLORS = ['#aeb8c4','#2768c4','#df7132']

def closest(loc, points):
    cell=vtk.vtkGenericCell(); cp=[0.,0.,0.]; cid=vtk.mutable(0); sub=vtk.mutable(0); d2=vtk.mutable(0.)
    dist=np.empty(len(points)); ids=np.empty(len(points),dtype=int); positions=np.empty_like(points)
    for i,x in enumerate(points):
        loc.FindClosestPoint(x.tolist(),cp,cell,cid,sub,d2)
        dist[i]=np.sqrt(float(d2)); ids[i]=int(cid); positions[i]=cp
    return dist,ids,positions

def membership(mesh, vertex_membership, face_ids, points):
    bary=trimesh.triangles.points_to_barycentric(mesh.triangles[face_ids], points)
    bary=np.clip(bary,0,1); bary/=bary.sum(axis=1)[:,None]
    return np.einsum('ni,nij->nj',bary,vertex_membership[mesh.faces[face_ids]])

def wquantile(values, weights, qs):
    order=np.argsort(values); v=np.asarray(values)[order]; w=np.asarray(weights)[order]
    return np.interp(np.asarray(qs)*w.sum(),np.cumsum(w),v)

def field_summary(df):
    a=df.vertex_area.to_numpy(); total=a.sum(); out={'area':float(total),'V':len(df)}
    for cutoff in [.183,1.8165,2.0]:
        out[f'area_t_gt_{cutoff}']=float(a[df.epsilon_kappa.to_numpy()>cutoff].sum()/total)
    for key in REGIONS:
        out[f'area_{key}']=float(a[(df.mechanism==key).to_numpy()].sum()/total)
    for col in ['guard_branch','clipped_unchanged','raw_at_min','raw_at_max','radius_at_min','radius_at_max']:
        out['area_'+col]=float(np.dot(a,df[col])/total)
    for col in ['h_raw','h_radius','h_cap','h_recovery','epsilon_kappa','delta','ratio']:
        v=df[col].to_numpy()
        out[col]={'min':float(v.min()),'max':float(v.max()),'area_mean':float(np.dot(a,v)/total),
                  'area_p01_p05_p50_p95_p99':wquantile(v,a,[.01,.05,.5,.95,.99]).tolist()}
        if col.startswith('h_'): out[col]['element_proxy']=float(np.sum(a/v**2))
    out['area_mean_abs_delta']=float(np.dot(a,np.abs(df.delta))/total)
    out['roundoff_unchanged_vertices']=int(((df.mechanism=='Unchanged') & (df.delta!=0)).sum())
    return out

def plot_field(name, mesh, df, dest):
    a=df.vertex_area.to_numpy(); t=df.epsilon_kappa.to_numpy(); r=df.ratio.to_numpy()
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    axes[0].hist(np.maximum(t,1e-8), bins=np.geomspace(1e-8,max(t.max(),3),90),weights=a/a.sum())
    axes[0].set_xscale('log'); axes[0].set_xlabel('epsilon * kappa')
    for cutoff in [.183,1.8165,2]: axes[0].axvline(cutoff,lw=.8,color='crimson')
    axes[1].hist(r,bins=np.geomspace(max(r.min()*.999,1e-5),max(r.max()*1.001,1.01),90),weights=a/a.sum())
    axes[1].set_xscale('log'); axes[1].set_xlabel('h_radius / h_raw'); axes[1].axvline(1,color='crimson',lw=.8)
    for ax in axes: ax.set_ylabel('Input area fraction per bin'); ax.grid(alpha=.2)
    fig.suptitle(name+' | area-weighted initial field'); fig.tight_layout();fig.savefig(dest/'histograms.png',dpi=180);plt.close(fig)
    labels=np.array([REGIONS.index(x) for x in df.mechanism]); rgba=np.array([matplotlib.colors.to_rgba(x) for x in COLORS])
    vertex_colors=(rgba[labels]*255).astype(np.uint8)
    colored=mesh.copy();colored.visual.vertex_colors=vertex_colors;colored.export(dest/'mechanism_regions.ply')
    face_colors=rgba[labels[mesh.faces]].mean(axis=1)
    fig=plt.figure(figsize=(13,5))
    for i,(az,el) in enumerate([(35,20),(125,15)]):
        ax=fig.add_subplot(1,2,i+1,projection='3d')
        surf=ax.plot_trisurf(*mesh.vertices.T,triangles=mesh.faces,shade=False,linewidth=0,antialiased=False)
        surf.set_facecolors(face_colors);ax.view_init(el,az);ax.set_box_aspect(np.ptp(mesh.vertices,axis=0));ax.set_axis_off()
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(color=c,label=r) for c,r in zip(COLORS,REGIONS)],loc='lower center',ncol=3)
    fig.suptitle(name+' | initial Cap / Recovery regions');fig.savefig(dest/'regions.png',dpi=180);plt.close(fig)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('results',type=Path);args=ap.parse_args();root=args.results.resolve()
    manifest=json.loads((root/'manifest.json').read_text(encoding='utf-8'))
    rows=[];local_rows=[];fields={}; failures=[]
    for source in manifest['data']['inputs']:
        name=Path(source['path']).stem; dest=root/name; dest.mkdir(exist_ok=True)
        mesh=trimesh.load(source['path'],force='mesh',process=False)
        runs=[r for r in manifest['runs'] if r['model']==name]
        good=[r for r in runs if r['status']=='success']
        if not good: continue
        df=pd.read_csv(Path(good[0]['path'])/'initial.csv')
        if len(df)!=len(mesh.vertices) or not np.allclose(df[['x','y','z']],mesh.vertices,rtol=0,atol=1e-14):
            raise ValueError('Initial diagnostics do not match original input vertices: '+name)
        for r in good:
            other=pd.read_csv(Path(r['path'])/'initial.csv')
            if not df.equals(other): raise ValueError('Initial fields differ across modes: '+name)
        fields[name]=field_summary(df)
        fields[name]['input']={'watertight':bool(mesh.is_watertight),'winding_consistent':bool(mesh.is_winding_consistent),
                               'is_volume':bool(mesh.is_volume),'V':len(mesh.vertices),'F':len(mesh.faces),'volume':float(mesh.volume)}
        plot_field(name,mesh,df,dest)
        vertex_membership=np.stack([(df.mechanism==r).to_numpy().astype(float) for r in REGIONS],axis=1)
        region_areas=np.einsum('n,nj->j',df.vertex_area.to_numpy(),vertex_membership)
        initial_counts=vertex_membership[mesh.faces].mean(axis=1).sum(axis=0)
        ref_loc,ref_poly=locator(mesh)
        ref_points,ref_ids=trimesh.sample.sample_surface(mesh,80000,seed=731)
        ref_members=membership(mesh,vertex_membership,ref_ids,ref_points)
        for run in runs:
            base={'model':name,'mode':run['mode'],'status':run['status'],'wall_seconds':run['wall_seconds']}
            if run['status']!='success': rows.append(base);continue
            cell=Path(run['path'])
            try:
                out=trimesh.load(cell/'output.ply',force='mesh',process=False)
                q=quality(out); row={**base,**q}
                row['watertight']=bool(out.is_watertight);row['winding_consistent']=bool(out.is_winding_consistent)
                row['components']=len(trimesh.graph.connected_components(out.face_adjacency,nodes=np.arange(len(out.faces))))
                row['volume_valid']=bool(mesh.is_volume and out.is_volume)
                row['relative_volume_change']=(float(out.volume/mesh.volume-1) if row['volume_valid'] else None)
                stdout=(cell/'stdout.log').read_text(encoding='utf-8')
                match=re.search(r'Remeshing runtime:\s*([\d.eE+-]+)',stdout)
                row['runtime_seconds']=float(match.group(1)) if match else None
                out_loc,out_poly=locator(out)
                out_points,_=trimesh.sample.sample_surface(out,80000,seed=343)
                fwd,_,_=closest(out_loc,ref_points); rev,_,_=closest(ref_loc,out_points);both=np.r_[fwd,rev]
                row.update(chamfer_squared_averaged=float((np.mean(fwd**2)+np.mean(rev**2))/2),
                           chamfer_unsquared_averaged=float(both.mean()),sampled_p99=float(np.percentile(both,99)),
                           sampled_max=float(both.max()),ref_out_p99=float(np.percentile(fwd,99)),out_ref_p99=float(np.percentile(rev,99)))
                md,ids,cp=closest(ref_loc,out.triangles_center)
                members=membership(mesh,vertex_membership,ids,cp);counts=members.sum(axis=0)
                normal_dot=np.einsum('ij,ij->i',out.face_normals,mesh.face_normals[ids])
                row['mapping_opposite_normal_fraction']=float(np.mean(normal_dot<0))
                row['mapping_centroid_distance_p99']=float(np.percentile(md,99))
                for j,region in enumerate(REGIONS):
                    w=ref_members[:,j]; area=region_areas[j]
                    local_rows.append({'model':name,'mode':run['mode'],'region':region,'input_area':area,
                        'input_face_equivalents':initial_counts[j],'output_face_equivalents':counts[j],
                        'density_per_input_area':counts[j]/area if area>0 else None,
                        'density_ratio_to_input':counts[j]/initial_counts[j] if initial_counts[j]>0 else None,
                        'ref_out_weighted_mean':np.dot(w,fwd)/w.sum() if w.sum()>0 else None,
                        'ref_out_weighted_p99':wquantile(fwd,w,[.99])[0] if w.sum()>0 else None})
                np.savez_compressed(cell/'sampled_distances.npz',ref_to_out=fwd,out_to_ref=rev)
                (cell/'metrics.json').write_text(json.dumps(row,indent=2),encoding='utf-8');rows.append(row)
                print('ANALYZED',name,run['mode'],'F',row['F'],'CD2',row['chamfer_squared_averaged'],flush=True)
            except Exception as exc:
                failures.append({'model':name,'mode':run['mode'],'error':repr(exc)})
                rows.append({**base,'analysis_error':repr(exc)})
                print('ANALYSIS FAILED',failures[-1],flush=True)
        pd.DataFrame(rows).to_csv(root/'summary.csv',index=False)
        pd.DataFrame(local_rows).to_csv(root/'local_density.csv',index=False)
        (root/'field_summary.json').write_text(json.dumps(fields,indent=2),encoding='utf-8')
    local=pd.DataFrame(local_rows)
    for name in fields:
        sub=local[local.model==name]
        if sub.empty: continue
        table=sub.pivot(index='mode',columns='region',values='density_per_input_area')
        ratio=table/table.loc['cgal-adaptive-pmap']
        ratio.to_csv(root/name/'density_ratio_to_pmap.csv')
        ax=ratio.plot.bar(figsize=(10,5),color=[COLORS[REGIONS.index(r)] for r in ratio.columns])
        ax.axhline(1,color='black',lw=.8);ax.set_ylabel('Regional density / Pmap density');ax.set_title(name)
        plt.xticks(rotation=20,ha='right');plt.tight_layout();plt.savefig(root/name/'local_density.png',dpi=180);plt.close()
    versions={m:__import__(m).__version__ for m in ['numpy','pandas','trimesh','matplotlib']}
    versions.update(python=sys.version,vtk=vtk.vtkVersion.GetVTKVersion(),failures=failures)
    (root/'analysis_environment.json').write_text(json.dumps(versions,indent=2),encoding='utf-8')
    if failures: raise SystemExit('Analysis failures recorded; see analysis_environment.json')

if __name__=='__main__': main()
