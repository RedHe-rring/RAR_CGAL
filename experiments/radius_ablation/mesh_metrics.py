"""Quality and VTK locator functions reused from tools/compare_wine.py."""
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

def poly(m):
 p=vtk.vtkPolyData(); pts=vtk.vtkPoints();pts.SetData(numpy_to_vtk(np.asarray(m.vertices,dtype=np.float64),deep=True));p.SetPoints(pts)
 cells=vtk.vtkCellArray(); nc=len(m.faces);cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0,3*nc+1,3,dtype=np.int64),deep=True),numpy_to_vtkIdTypeArray(np.asarray(m.faces,dtype=np.int64).ravel(),deep=True));p.SetPolys(cells)
 return p

def locator(m):
 p=poly(m);loc=vtk.vtkStaticCellLocator();loc.SetDataSet(p);loc.BuildLocator();return loc,p

def quality(mesh):
 v=mesh.vertices;f=mesh.faces
 a=v[f[:,0]];b=v[f[:,1]];c=v[f[:,2]]
 l=np.stack([np.linalg.norm(a-b,axis=1),np.linalg.norm(b-c,axis=1),np.linalg.norm(c-a,axis=1)],axis=1)
 lsort=np.sort(l,axis=1)
 lmax=lsort[:,2];lmin=lsort[:,0];
 area=0.5*np.linalg.norm(np.cross(b-a,c-a),axis=1)
 angles=[]
 for u,vv,w in [(l[:,0],l[:,2],l[:,1]),(l[:,0],l[:,1],l[:,2]),(l[:,1],l[:,2],l[:,0])]:
  angles.append(np.degrees(np.arccos(np.clip((u*u+vv*vv-w*w)/(2*u*vv),-1.,1.))))
 ang=np.stack(angles,axis=1)
 min_ang=np.min(ang,axis=1);max_ang=np.max(ang,axis=1)
 tri_q=4*np.sqrt(3)*area/(l*l).sum(axis=1)
 ar=np.sqrt(3)*lmax*lmax/(4*area)
 es=mesh.edges_unique_length
 return {'V':len(v),'F':len(f),'area':float(mesh.area),'vol':float(mesh.volume),'minAngP05':float(np.percentile(min_ang,5)), 'minAngP50':float(np.median(min_ang)),'minAngMin':float(min_ang.min()),'maxAngP95':float(np.percentile(max_ang,95)),'fracMinAngUnder30':float(np.mean(min_ang<30)),'fracMaxAngOver90':float(np.mean(max_ang>90)),'triQ05':float(np.percentile(tri_q,5)),'AR95':float(np.percentile(ar,95)), 'EdgeCV':float(np.std(es)/np.mean(es)),'edgeMean':float(np.mean(es)),'edgeP05':float(np.percentile(es,5)),'edgeP95':float(np.percentile(es,95))}
