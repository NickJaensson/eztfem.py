close all; clear

addpath('subs/')

eztfempath = "~/Desktop/eztfem/";
addpath(eztfempath);
addpath(append(eztfempath,"src"))
addpath(append(eztfempath,"addons/plotlib"))
addpath(append(eztfempath,"addons/meshes"))
addpath(append(eztfempath,"addons/poisson"))

addpath(append(eztfempath,"examples/poisson"))  % use func.m from poisson

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


%% filename for python test file
global fn 
fn = "~/Desktop/eztfem.py/dotest_poisson.py"; 


%% run the problem in eztfem

mesh_ez = quadrilateral2d([3,2],'quad9','origin',[1,1],'length',[4,3]);

elementdof_ez = [1,1,1,1,1,1,1,1,1; 2,2,2,2,2,2,2,2,2]' ;
problem_ez = problem_definition(mesh_ez,elementdof_ez,'nphysq',1);
[user_ez.xr,user_ez.wg] = gauss_legendre('quad','n', 3 );
[user_ez.phi,user_ez.dphi] = basis_function('quad','Q2', user_ez.xr );
user_ez.coorsys = 0 ;
user_ez.alpha = 1 ;
user_ez.funcnr = 4 ;
user_ez.func = @func ;
[A_ez,f_ez] = build_system ( mesh_ez, problem_ez, @poisson_elem, user_ez);
iess_ez = define_essential ( mesh_ez, problem_ez, 'curves', [1 2 3 4], 'degfd', 1 ) ;

uess_ez = fill_system_vector ( mesh_ez, problem_ez, 'curves', [1,2], @func, 'funcnr', 3 );
uess_ez = fill_system_vector ( mesh_ez, problem_ez, 'curves', [3,4], @func, 'funcnr', 3, 'fin', uess_ez );

[A_ez2, f_ez2] = apply_essential ( A_ez, f_ez, uess_ez, iess_ez );

u_ez = A_ez2\f_ez2 ;

user_ez2 = user_ez;
xr_ez = refcoor_nodal_points ( mesh_ez ) ;
[user_ez2.phi,user_ez2.dphi] = basis_function('quad','Q2', xr_ez ) ;
user_ez2.u = u_ez ;
gradu_ez = deriv_vector ( mesh_ez, problem_ez, @poisson_deriv, user_ez2 ) ;


%% define the same commands for pytfem

cmd_mesh_py =       "    mesh_py = quadrilateral2d([3,2],'quad9',origin=np.array([1,1]),length=np.array([4,3]))";

cmd_elementdof_py = "    elementdof_py = np.array([[1,1,1,1,1,1,1,1,1],[2,2,2,2,2,2,2,2,2]]).T"; 
cmd_problem_py =    "    problem_py = Problem(mesh_py,elementdof_py,nphysq=1)";
cmd_fill_user_py =  "    user_py = User();" + ...
                    "    user_py.coorsys = 0;"+...
                    "    user_py.alpha = 1;"+...
                    "    user_py.funcnr = 4; "+...
                    "    user_py.func = func";
cmd_gauss_py =      "    user_py.xr, user_py.wg = gauss_legendre('quad',n=3 )";
cmd_basis_py =      "    user_py.phi, user_py.dphi = basis_function('quad','Q2', user_py.xr );";

cmd_build_sys_py =  "    A_py,f_py = build_system ( mesh_py, problem_py, poisson_elem, user_py)";
cmd_define_ess_py = "    iess_py = define_essential ( mesh_py, problem_py,'curves', [0,1,2,3], degfd=0 );";

cmd_fill_sys_py =   "    uess_py = fill_system_vector ( mesh_py, problem_py, 'curves', [0,1], func, funcnr=3 );"+...
                    "    uess_py = fill_system_vector ( mesh_py, problem_py, 'curves', [2,3], func, funcnr=3, fin=uess_py )";

cmd_apply_ess_py =  "    A_py2, f_py2, _ = apply_essential ( A_py, f_py, uess_py, iess_py )";

cmd_solve_py     =  "    u_py = spsolve(A_py2.tocsr(), f_py2)";

cmd_deriv_vector =  "    user_py2 = user_py;" + ...
                    "    xr_py = refcoor_nodal_points ( mesh_py );"+...
                    "    user_py2.phi, user_py2.dphi = basis_function('quad','Q2', xr_py );"+...
                    "    user_py2.u = u_py;"+...
                    "    gradu_py = deriv_vector ( mesh_py, problem_py, poisson_deriv, user_py2 )";


%% write some header stuff

writelines("# run with: python -m unittest dotest.py",fn);
mywritelines("import numpy as np")
mywritelines("import unittest");
mywritelines("from src.distribute_elements import distribute_elements");
mywritelines("from src.quadrilateral2d import quadrilateral2d");
mywritelines("from src.mesh_class import Mesh, Geometry");
mywritelines("from src.problem_class import Problem");
mywritelines("from src.user_class import User");
mywritelines("from src.gauss_legendre import gauss_legendre");
mywritelines("from src.basis_function import basis_function");
mywritelines("from src.build_system import build_system");
mywritelines("from addons.poisson.poisson_elem import poisson_elem");
mywritelines("from addons.poisson.poisson_deriv import poisson_deriv");

mywritelines("from src.define_essential import define_essential");

mywritelines("from src.fill_system_vector import fill_system_vector");
mywritelines("from src.apply_essential import apply_essential");
mywritelines("from src.vector_class import Vector");
mywritelines("from src.deriv_vector import deriv_vector");
mywritelines("from src.refcoor_nodal_points import refcoor_nodal_points");

mywritelines("from scipy.sparse.linalg import spsolve")
mywritelines("from examples.poisson.func import func");

mywritelines("class TestPytfem(unittest.TestCase):");


%% test for quadrilateral2d

mywritelines("  def test_quadrilaterial2d(self):");
mywritelines(cmd_mesh_py);
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
mywritelines("    # compensate for zero-based indexing");
mywritelines("    mesh_ez.topology = mesh_ez.topology - 1 # Python indexing");
mywritelines("    mesh_ez.points = mesh_ez.points - 1 # Python indexing");
mywritelines("    self.assertTrue(mesh_py==mesh_ez,'quadrilateral2d failed test!' )");


%% test for problem_definition

mywritelines("  def test_problem_definition(self):");
mywritelines(cmd_mesh_py);
mywritelines(cmd_elementdof_py);
mywritelines(cmd_problem_py);
mywritelines("    problem_ez = Problem(mesh_py,elementdof_py)");
write_attrib("    ",problem_ez,"problem_ez")
mywritelines("    self.assertTrue(problem_py==problem_ez,'problem_definition failed test!' )");


%% test for user equivalence 
% NOTE: not all attributes checked in Python code, see user.py

mywritelines("  def test_user(self):");
mywritelines("    user_ez = User()");
write_attrib("    ",user_ez,"user_ez")
mywritelines(cmd_fill_user_py);
mywritelines(cmd_gauss_py);
mywritelines(cmd_basis_py);
mywritelines("    self.assertTrue(user_py==user_ez,'users failed test!' )");


%% test for gauss_legendre

mywritelines("  def test_gauss_legendre(self):");
mywritelines("    user_ez = User()");
write_attrib("    ",user_ez,"user_ez")
mywritelines("    user_py = User()");
mywritelines(cmd_gauss_py);
mywritelines("    check1=np.allclose(user_py.wg,user_ez.wg,atol=1e-15,rtol=0)")
mywritelines("    check2=np.allclose(user_py.wg,user_ez.wg,atol=1e-15,rtol=0)")
mywritelines("    self.assertTrue(check1 and check2,'gauss_legendre failed test!' )");


%% test for basis_function

mywritelines("  def test_basis_function(self):");
mywritelines("    user_ez = User()");
write_attrib("    ",user_ez,"user_ez")
mywritelines("    user_py = User()");
mywritelines(cmd_gauss_py);
mywritelines(cmd_basis_py);
mywritelines("    check1=np.allclose(user_py.phi,user_ez.phi,atol=1e-15,rtol=0)")
mywritelines("    check2=np.allclose(user_py.dphi,user_ez.dphi,atol=1e-15,rtol=0)")
mywritelines("    self.assertTrue(check1 and check2,'basis_functions failed test!' )");


%% test for build_system

mywritelines("  def test_build_system(self):");
mywritelines(cmd_mesh_py);
mywritelines(cmd_elementdof_py);
mywritelines(cmd_problem_py);
mywritelines("    problem_ez = Problem(mesh_py,elementdof_py)");
write_attrib("    ",problem_ez,"problem_ez")
mywritelines("    user_ez = User()");
write_attrib("    ",user_ez,"user_ez")
mywritelines(cmd_fill_user_py);
mywritelines(cmd_gauss_py);
mywritelines(cmd_basis_py);
write2Darr_r("    ",full(A_ez),"A_ez")
write1Darr_r("    ",f_ez,"f_ez")
mywritelines(cmd_build_sys_py);
mywritelines("    check1=np.allclose(A_py.toarray(),A_ez,atol=1e-12,rtol=0)")
mywritelines("    check2=np.allclose(f_py,f_ez,atol=1e-12,rtol=0)")
mywritelines("    self.assertTrue(check1 and check2,'build_system failed test!' )");


%% test for define_essential

mywritelines("  def test_define_essential(self):");
mywritelines(cmd_mesh_py);
mywritelines(cmd_elementdof_py);
mywritelines(cmd_problem_py);
write1Darr_i("    ",iess_ez,"iess_ez")
mywritelines(cmd_define_ess_py);
mywritelines("    self.assertTrue((iess_ez-1==iess_py).all(),'define_essential failed test!' )");


%% test for fill_sysvector

mywritelines("  def test_fill_system_vector(self):");
mywritelines(cmd_mesh_py);
mywritelines(cmd_elementdof_py);
mywritelines(cmd_problem_py);
write1Darr_r("    ",uess_ez,"uess_ez")
mywritelines(cmd_define_ess_py);
mywritelines(cmd_fill_sys_py);
mywritelines("    self.assertTrue(np.allclose(uess_py,uess_ez,atol=1e-15,rtol=0),'fill_system_vector failed test!' )");


%% test for apply_essential

mywritelines("  def test_apply_essential(self):");
mywritelines(cmd_mesh_py);
mywritelines(cmd_elementdof_py);
mywritelines(cmd_problem_py);
mywritelines(cmd_fill_user_py);
mywritelines(cmd_gauss_py);
mywritelines(cmd_basis_py);
mywritelines(cmd_build_sys_py);
mywritelines(cmd_define_ess_py);
mywritelines(cmd_fill_sys_py);
mywritelines(cmd_apply_ess_py);
write2Darr_r("    ",full(A_ez2),"A_ez2")
write1Darr_r("    ",f_ez2,"f_ez2")
mywritelines("    check1=np.allclose(A_py2.toarray(),A_ez2,atol=1e-12,rtol=0)")
mywritelines("    check2=np.allclose(f_py2,f_ez2,atol=1e-12,rtol=0)")
mywritelines("    self.assertTrue(check1 and check2,'apply_essential failed test!' )");


%% test for solve

mywritelines("  def test_solve(self):");
write1Darr_r("    ",u_ez,"u_ez")

mywritelines(cmd_mesh_py);
mywritelines(cmd_elementdof_py);
mywritelines(cmd_problem_py);
mywritelines(cmd_fill_user_py);
mywritelines(cmd_gauss_py);
mywritelines(cmd_basis_py);
mywritelines(cmd_build_sys_py);
mywritelines(cmd_define_ess_py);
mywritelines(cmd_fill_sys_py);
mywritelines(cmd_apply_ess_py);
mywritelines(cmd_solve_py);
mywritelines("    #print('max diff = ',(abs(u_ez-u_py)).max())");

mywritelines("    self.assertTrue(np.allclose(u_py,u_ez,atol=1e-12,rtol=0)," + ...
    "'solve failed test, max diff = '+str((abs(u_ez-u_py)).max()) )");


%% tests for deriv_vector

mywritelines("  def test_deriv_vector(self):");
mywritelines("    gradu_ez = Vector()");
write_attrib("    ",gradu_ez,"gradu_ez")
mywritelines("    gradu_ez.vec += -1 # compensate for Python indexing"); 

mywritelines(cmd_mesh_py);
mywritelines(cmd_elementdof_py);
mywritelines(cmd_problem_py);
mywritelines(cmd_fill_user_py);
mywritelines(cmd_gauss_py);
mywritelines(cmd_basis_py);
mywritelines(cmd_build_sys_py);
mywritelines(cmd_define_ess_py);
mywritelines(cmd_fill_sys_py);
mywritelines(cmd_apply_ess_py);
mywritelines(cmd_solve_py);
mywritelines(cmd_deriv_vector);
mywritelines("    self.assertTrue(gradu_ez==gradu_py," + ...
    "'deriv_vector failed test, max diff = '+str((abs(gradu_ez.u-gradu_py.u)).max()) )");


%% helper functions

function mywritelines(str)
    global fn
    writelines(str,fn,WriteMode="append");
end

function write1Darr_r(skip,arr,name)
    mywritelines("    "+name+" = np.array([");
    for i=1:length(arr)
        mywritelines(skip+sprintf('%25.16e',arr(i))+",");
    end
    mywritelines("    ])");
end

function write1Darr_i(skip,arr,name)
    mywritelines("    "+name+" = np.array([");
    for i=1:length(arr)
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
    mywritelines("    ])");
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
            elseif all(mod(tmp,1)<1e-15,'all') % array of integers
                if ndims(tmp) == 2 
                    if size(tmp,1)==1 || size(tmp,2)==1
                        write1Darr_i(skip,tmp,name+"."+fns{k});
                    else 
                        write2Darr_i(skip,tmp,name+"."+fns{k});
                    end
                else
                    write3Darr_i(skip,tmp,name+"."+fns{k});
                end
            else 
                if ndims(tmp) == 2 
                    if size(tmp,1)==1 || size(tmp,2)==1
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