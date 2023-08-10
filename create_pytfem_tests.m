close all; clear

% add eztfem to path

eztfempath = "~/Desktop/eztfem/";
addpath_eztfem(eztfempath)


% filename to which to write the python test file

global fn 
fn = "~/Desktop/pytfem/dotest.py";


% write some header stuff

writelines("# run with: python -m unittest dotest.py",fn);
mywritelines("import numpy as np")
mywritelines("import unittest");
mywritelines("from distribute_elements import distribute_elements");
mywritelines("from quadrilateral2d import quadrilateral2d");
mywritelines("from mesh_defs import Mesh, Geometry");

mywritelines("class Testsol_pytfem(unittest.TestCase):");


% test for distribute_elements

cmd = "distribute_elements(8,1,2)";
sol_eztfem = eval(cmd);

mywritelines("  def test_distribute_elements(self):");
mywritelines("    sol_pytfem = "+cmd);
write1Darr_r("    ",sol_eztfem,"sol_eztfem")
mywritelines("    self.assertTrue(np.allclose(sol_pytfem,sol_eztfem,atol=1e-15,rtol=0),'distribute_elements failed test!' )");


% test for quadrilateral2d

cmd = "quadrilateral2d([3,2],'quad9')";
sol_eztfem = eval(cmd);



mywritelines("  def test_quadrilaterial2d(self):");
mywritelines("    sol_pytfem = "+cmd);

mywritelines("    sol_eztfem = Mesh()");

% fieldnames = fieldnames(sol_eztfem);
% for k=1:numel(fieldnames)
%     tmp = sol_eztfem.(fieldnames{k}));
%     if isnumeric(tmp)
%         if iscalar(tmp)
%             mywritelines("    sol_eztfem.ndim = "+string(sol_eztfem.ndim));
% 
% 
%         end
%     end
% end



mywritelines("    sol_eztfem.ndim = "+string(sol_eztfem.ndim));
mywritelines("    sol_eztfem.nnodes = "+string(sol_eztfem.nnodes));
mywritelines("    sol_eztfem.elshape = "+string(sol_eztfem.elshape));
mywritelines("    sol_eztfem.nelem = "+string(sol_eztfem.nelem));
mywritelines("    sol_eztfem.elnumnod = "+string(sol_eztfem.elnumnod));
mywritelines("    sol_eztfem.npoints = "+string(sol_eztfem.npoints));
mywritelines("    sol_eztfem.ncurves = "+string(sol_eztfem.ncurves));
write2Darr_i("    ",sol_eztfem.topology,"sol_eztfem.topology")
write2Darr_r("    ",sol_eztfem.coor,"sol_eztfem.coor")
write1Darr_i("    ",sol_eztfem.points,"sol_eztfem.points")




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

mywritelines("    # compensate for zero-based indexing");
mywritelines("    sol_eztfem.topology = sol_eztfem.topology - 1");
mywritelines("    sol_eztfem.points = sol_eztfem.points - 1");

% for i=1:sol_eztfem.ncurves
% 
% 
% 
% end



mywritelines("    self.assertTrue(sol_pytfem==sol_eztfem,'quadrilateral2d failed test!' )");




function addpath_eztfem(eztfempath)

    addpath(eztfempath);
    addpath(append(eztfempath,"src"))
    addpath(append(eztfempath,"addons/plotlib"))
    addpath(append(eztfempath,"addons/meshes"))
    addpath(append(eztfempath,"addons/poisson"))
    addpath(append(eztfempath,"addons/stokes"))

end

function mywritelines(str)

    global fn

    writelines(str,fn,WriteMode="append");

end

function write1Darr_r(skip,arr,name)

    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines(skip+sprintf('%25.16e',arr(i))+",");
    end
    mywritelines("    ])");
end

function write1Darr_i(skip,arr,name)

    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines(skip+sprintf('%12i,',arr(i)));
    end
    mywritelines("    ],dtype=int)");

end


function write2Darr_r(skip,arr,name)

    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines(skip+"["+sprintf('%25.16e,',arr(i,:))+"],");
    end
    mywritelines("    ])");
end

function write2Darr_i(skip,arr,name)

    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines(skip+"["+sprintf('%12i,',arr(i,:))+"],");
    end
    mywritelines("    ],dtype=int)");

end
