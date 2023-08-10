close all; clear

% add eztfem to path

eztfempath = "~/Desktop/eztfem/"
addpath(append(eztfempath))
addpath(append(eztfempath,"src"))
addpath(append(eztfempath,"addons/plotlib"))
addpath(append(eztfempath,"addons/meshes"))
addpath(append(eztfempath,"addons/poisson"))
addpath(append(eztfempath,"addons/stokes"))


% filename to which to write the python test file

fn = "~/Desktop/pytfem/dotest.py";


% write some header stuff

writelines("# run with: python -m unittest dotest.py",fn);
writelines("import numpy as np",fn,WriteMode="append");
writelines("import unittest",fn,WriteMode="append");
writelines("from distribute_elements import distribute_elements",fn,WriteMode="append");
writelines("from quadrilateral2d import quadrilateral2d",fn,WriteMode="append");
writelines("from mesh_defs import Mesh, Geometry",fn,WriteMode="append");

writelines("class Testsol_pytfem(unittest.TestCase):",fn,WriteMode="append");


% test for distribute_elements

cmd = "distribute_elements(8,1,2)";
sol_eztfem = eval(cmd);

writelines("  def test_distribute_elements(self):",fn,WriteMode="append");
writelines(append("    sol_pytfem = ",cmd),fn,WriteMode="append");
writelines("    sol_eztfem = np.array([",fn,WriteMode="append");
for i=1:size(sol_eztfem,1)
  writelines(append("     ",sprintf('%25.16e',sol_eztfem(i)),","),fn,WriteMode="append");
end
writelines("    ])",fn,WriteMode="append");
writelines("    self.assertTrue(np.allclose(sol_pytfem,sol_eztfem,atol=1e-15,rtol=0),'distribute_elements failed test!' )",fn,WriteMode="append");


% test for quadrilateral2d

cmd = "quadrilateral2d([1,2],'quad9')";
sol_eztfem = eval(cmd);

writelines("  def test_quadrilaterial2d(self):",fn,WriteMode="append");
writelines(append("    sol_pytfem = ",cmd),fn,WriteMode="append");
writelines("    sol_eztfem = Mesh()",fn,WriteMode="append");
writelines(append("    sol_eztfem.ndim = ",string(sol_eztfem.ndim)),fn,WriteMode="append");
writelines(append("    sol_eztfem.nnodes = ",string(sol_eztfem.nnodes)),fn,WriteMode="append");
writelines(append("    sol_eztfem.elshape = ",string(sol_eztfem.elshape)),fn,WriteMode="append");
writelines(append("    sol_eztfem.nelem = ",string(sol_eztfem.nelem)),fn,WriteMode="append");
writelines(append("    sol_eztfem.elnumnod = ",string(sol_eztfem.elnumnod)),fn,WriteMode="append");
writelines(append("    sol_eztfem.npoints = ",string(sol_eztfem.npoints)),fn,WriteMode="append");
writelines(append("    sol_eztfem.ncurves = ",string(sol_eztfem.ncurves)),fn,WriteMode="append");

writelines("    sol_eztfem.topology = np.array([",fn,WriteMode="append");
for i=1:size(sol_eztfem.topology,1)
  writelines(append("     [",sprintf('%12i,',sol_eztfem.topology(i,:)),"],"),fn,WriteMode="append");
end
writelines("    ],dtype=int)",fn,WriteMode="append");

writelines("    sol_eztfem.coor = np.array([",fn,WriteMode="append");
for i=1:size(sol_eztfem.coor,1)
  writelines(append("     [",sprintf('%25.16e,',sol_eztfem.coor(i,:)),"],"),fn,WriteMode="append");
end
writelines("    ])",fn,WriteMode="append");

writelines("    sol_eztfem.points = np.array([",fn,WriteMode="append");
for i=1:size(sol_eztfem.points,1)
  writelines(append("     ",sprintf('%12i,',sol_eztfem.points(i)),""),fn,WriteMode="append");
end
writelines("    ],dtype=int)",fn,WriteMode="append");

%         mesh.curves[0]= Geometry(elshape=2,ndim=2,elnumnod=3,nnodes=3,nelem=1)
%         mesh.curves[0].nodes = np.array([1,2,3],dtype=int)
%         mesh.curves[0].topology = np.array([[[1,2,3],[1,2,3]])
%  
%         mesh.curves[1]= Geometry(elshape=2,ndim=2,elnumnod=3,nnodes=5,nelem=2)
%         mesh.curves[1].nodes = np.array([3,6,9,12,15],dtype=int)
% 
%         mesh.curves[2]= Geometry(elshape=2,ndim=2,elnumnod=3,nnodes=3,nelem=1)
%         mesh.curves[2].nodes = np.array([15,14,13],dtype=int)
% 
%         mesh.curves[3]= Geometry(elshape=2,ndim=2,elnumnod=3,nnodes=5,nelem=2)
%         mesh.curves[3].nodes = np.array([13,10,7,4,1],dtype=int)

writelines("    # compensate for zero-based indexing",fn,WriteMode="append");
writelines("    sol_eztfem.topology = sol_eztfem.topology - 1",fn,WriteMode="append");
writelines("    sol_eztfem.points = sol_eztfem.points - 1",fn,WriteMode="append");

for i=1:sol_eztfem.ncurves



end




writelines("    self.assertTrue(sol_pytfem==sol_eztfem,'quadrilateral2d failed test!' )",fn,WriteMode="append");