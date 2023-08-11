close all; clear

%% Matlab file to generate testing code using eztfem for pytfem 
% NOTE: arrays are written "as is" (looping of the first index, then the
% second etc.). It is assumed that in Python they are also read like this
% 
% NOTE2: Python has real 1D arrays, whereas Matlab 1D arrays are N x 1
% (i.e., ndims([1,2]) gives as output 2)
% This can lead to problems when comparing arrays, since these are not
% equivalent. 
% 
% This is handled on the Matlab side by:
%  - write scalars as scalars (using isscalar)
%  - write 2D arrays (so can be 1D or "true" 2D) as 1D when second dim == 1
%  - write 3D arrays always as a full 3D array
%
% And on the Python side:
%  - in Python: use numpy.squeeze on all 2D arrays to get rid of dim == 1


%% add eztfem to path

eztfempath = "~/Desktop/eztfem/";
addpath_eztfem(eztfempath)

global fn 
fn = "~/Desktop/pytfem/dotest.py"; % filename for python test file


%% write some header stuff

writelines("# run with: python -m unittest dotest.py",fn);
mywritelines("import numpy as np")
mywritelines("import unittest");
mywritelines("from src.distribute_elements import distribute_elements");
mywritelines("from src_test.quadrilateral2d import quadrilateral2d");
mywritelines("from src_test.mesh_defs import Mesh, Geometry");
% mywritelines("from problem_definition import Problem");
% mywritelines("from gauss_legendre import gauss_legendre");
% mywritelines("from basis_function import basis_function");
% mywritelines("from user import User");

mywritelines("class TestPytfem(unittest.TestCase):");


%% test for distribute_elements

cmd_grid = "distribute_elements(8,1,3)";
grid_ez = eval(cmd_grid);

mywritelines("  def test_distribute_elements(self):");
mywritelines("    grid_py = "+cmd_grid);
write1Darr_r("    ",grid_ez,"grid_ez")
mywritelines("    self.assertTrue(np.allclose(grid_py,grid_ez," + ...
    "atol=1e-15,rtol=0),'distribute_elements failed test!' )");


%% run the problem in Matlab

% define the commands to run in Matlab and python
cmd_mesh_ez = "quadrilateral2d([1,2],'quad9')";
cmd_mesh_py = "quadrilateral2d([1,2],'quad9')";

% elementdof=[1,1,1,1,1,1,1,1,1;
%             2,2,2,2,2,2,2,2,2]' ;
% cmd_problem_ez = "    problem_definition(mesh_ez,elementdof,'nphysq',1);";
% cmd_problem_py = "    Problem(mesh_py,elementdof_py,nphysq=1);";
% 
% cmd_gauss_ez = "gauss_legendre('quad','n', 1 ) ;";
% cmd_gauss_py = "gauss_legendre('quad',n=1 )";
% 
% cmd_basis_ez = "basis_function('quad','Q2', user_ez.xr ) ;";
% cmd_basis_py = "basis_function('quad','Q2', user_py.xr )";


% run the Matlab code
mesh_ez = eval(cmd_mesh_ez);
% problem_ez = eval(cmd_problem_ez);
% [user_ez.xr,user_ez.wg] = eval(cmd_gauss_ez);
% [user_ez.phi,user_ez.dphi] = eval(cmd_basis_ez);

%% test for quadrilateral2d

mywritelines("  def test_quadrilaterial2d(self):");

% generate the mesh in pytfem
mywritelines("    mesh_py = "+cmd_mesh_py);

% copy the mesh from eztfem to pytfem
mywritelines("    mesh_ez = Mesh()");
write_attrib("    ",mesh_ez,"mesh_ez")
for i=1:mesh_ez.ncurves
    mywritelines("    mesh_ez.curves.append(Geometry())");
    write_attrib("    ",mesh_ez.curves(i),"mesh_ez.curves["+string(i-1)+"]")
    mywritelines("    mesh_ez.curves["+string(i-1)+"].topology = mesh_ez.curves["...
        +string(i-1)+"].topology - 1 # Python indexing");
    mywritelines("    mesh_ez.curves["+string(i-1)+"].nodes = mesh_ez.curves["...
        +string(i-1)+"].nodes - 1 # Python indexing"); 
end

% compensate the zero-based indexing
mywritelines("    # compensate for zero-based indexing");
mywritelines("    mesh_ez.topology = mesh_ez.topology - 1 # Python indexing");
mywritelines("    mesh_ez.points = mesh_ez.points - 1 # Python indexing");

% check for equivalence
mywritelines("    self.assertTrue(mesh_py==mesh_ez,'quadrilateral2d failed test!' )");


% %% test for problem_definition
% 
% mywritelines("  def test_problem_definition(self):");
% mywritelines("    mesh_py = "+cmd_mesh_py);
% write2Darr_i("    ",elementdof,"elementdof_py")
% mywritelines("    problem_py = "+cmd_problem_py);
% mywritelines("    problem_ez = Problem(mesh_py,elementdof_py)");
% write_attrib("    ",problem_ez,"problem_ez")
% 
% mywritelines("    self.assertTrue(problem_py==problem_ez,'problem_definition failed test!' )");
% 
% 
% %% test for gauss_legendre and basis function
% 
% mywritelines("  def test_gauss_legendre(self):");
% 
% mywritelines("    user_ez = User()");
% write_attrib("    ",user_ez,"user_ez")
% 
% mywritelines("    user_py = User()");
% mywritelines("    user_py.xr, user_py.wg = "+cmd_gauss_py);
% mywritelines("    user_py.phi, user_py.dphi = "+cmd_basis_py);
% 
% mywritelines("    self.assertTrue(user_py==user_ez,'problem_definition failed test!' )");


%% helper functions %%%%%%%%%%%%%%%%%%

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

function write3Darr_i(skip,arr,name)

    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines("    [")
        for j=1:size(arr,2)
            mywritelines(skip+"["+sprintf('%12i,',arr(i,j,:))+"],");
        end
        mywritelines("    ],")
    end
    mywritelines("    ],dtype=int)");

end

function write3Darr_r(skip,arr,name)

    mywritelines("    "+name+" = np.array([");
    for i=1:size(arr,1)
        mywritelines("    [")
        for j=1:size(arr,2)
            mywritelines(skip+"["+sprintf('%25.16e,',arr(i,j,:))+"],");
        end
        mywritelines("    ],")
    end
    mywritelines("    ],dtype=int)");

end

function write_attrib(skip,struct,name)

    fns = fieldnames(struct);
    for k=1:numel(fns)
        tmp = struct.(fns{k});
        if isnumeric(tmp)
            if ndims(tmp) > 3
                error("error write_attrib: dims("+fns{k}+") > 3")
            end
            if isscalar(tmp)
                mywritelines(skip+name+"."+fns{k}+" = "+string(tmp));
            elseif all(mod(tmp,1)==0,'all') % array of integers
                if ndims(tmp) == 2 
                    if size(tmp,2)==1
                        write1Darr_i(skip,tmp,name+"."+fns{k});
                    else 
                        write2Darr_i(skip,tmp,name+"."+fns{k});
                    end
                else
                    write3Darr_i(skip,tmp,name+"."+fns{k});
                end
            else 
                if ndims(tmp) == 2 
                    if size(tmp,2)==1
                        write1Darr_r(skip,tmp,name+"."+fns{k}); 
                    else
                        write2Darr_r(skip,tmp,name+"."+fns{k});
                    end
                else
                    write3Darr_r(skip,tmp,name+"."+fns{k});
                end
            end
        end
    end

end