close all; clear

addpath("subs/");

eztfempath = "~/Desktop/eztfem/";
addpath(eztfempath);
addpath(append(eztfempath,"src"))
addpath(append(eztfempath,"addons/plotlib"))
addpath(append(eztfempath,"addons/meshes"))
addpath(append(eztfempath,"addons/stokes"))

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
fn = "~/Desktop/eztfem.py/dotest_stokes.py"; 


%% run the problem in eztfem

problemtype = "stokes";

mesh_ez = quadrilateral2d([3,2],'quad9','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);

elementdof_ez = [2,2,2,2,2,2,2,2,2;
                 1,0,1,0,1,0,1,0,0;
                 1,1,1,1,1,1,1,1,1]' ;
problem_ez = problem_definition(mesh_ez,elementdof_ez,'nphysq',2);
[user_ez.xr,user_ez.wg] = gauss_legendre('quad','n', 3 );
[user_ez.phi,user_ez.dphi] = basis_function('quad','Q2', user_ez.xr );
[user_ez.psi,~] = basis_function('quad','Q1', user_ez.xr );
user_ez.coorsys = 0 ;
user_ez.mu = 1 ;
user_ez.funcnr = 0 ;
user_ez.func = @func ;
[A_ez,f_ez] = build_system ( mesh_ez, problem_ez, @stokes_elem, user_ez);
iess_ez = define_essential ( mesh_ez, problem_ez, 'curves', [1 2 3 4], 'degfd', 1 ) ;
iess_ez = define_essential ( mesh_ez, problem_ez, 'curves', [1 2 3 4], 'degfd', 2, 'iessp', iess_ez ) ;
iess_ez = define_essential ( mesh_ez, problem_ez, 'points', 1, 'physq', 2, 'iessp', iess_ez ) ;

uess_ez = fill_system_vector ( mesh_ez, problem_ez, 'curves', [1,2], @func, 'funcnr', 3 );
uess_ez = fill_system_vector ( mesh_ez, problem_ez, 'curves', [3,4], @func, 'funcnr', 3, 'fin', uess_ez );

[A_ez2, f_ez2] = apply_essential ( A_ez, f_ez, uess_ez, iess_ez );

u_ez = A_ez2\f_ez2 ;

user_ez2 = user_ez;
xr_ez = refcoor_nodal_points ( mesh_ez ) ;
[user_ez2.psi] = basis_function('quad','Q1', xr_ez ) ;
user_ez2.u = u_ez ;
pressure_ez = deriv_vector ( mesh_ez, problem_ez, @stokes_pressure, user_ez2 ) ;

[user_ez2.phi,user_ez2.dphi]=basis_function('quad','Q2', xr_ez ) ;
user_ez2.comp = 7 ; % divu, divergence of the velocity field
divu_ez = deriv_vector ( mesh_ez, problem_ez, @stokes_deriv, user_ez2 ) ;
user_ez2.comp = 8 ; % gammadot, effective strain rate = sqrt(2II_D) 
gammadot_ez = deriv_vector ( mesh_ez, problem_ez, @stokes_deriv, user_ez2 ) ;


%% define the same commands for pytfem

cmd_mesh_py =       "    mesh_py = quadrilateral2d([3,2],'quad9',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";

cmd_elementdof_py = "    elementdof_py = np.array([[2,2,2,2,2,2,2,2,2],[1,0,1,0,1,0,1,0,0],[1,1,1,1,1,1,1,1,1]]).T";
cmd_problem_py =    "    problem_py = Problem(mesh_py,elementdof_py,nphysq=2)";
cmd_fill_user_py =  "    user_py = User();" + ...
                    "    user_py.coorsys = 0;"+...
                    "    user_py.mu = 1;"+...
                    "    user_py.funcnr = 0; "+...
                    "    user_py.func = func";
cmd_gauss_py =      "    user_py.xr, user_py.wg = gauss_legendre('quad',n=3 )";
cmd_basis_py =      "    user_py.phi, user_py.dphi = basis_function('quad','Q2', user_py.xr );"+...
                    "    user_py.psi, _ = basis_function('quad','Q1', user_py.xr )";

cmd_build_sys_py =  "    A_py,f_py = build_system ( mesh_py, problem_py, stokes_elem, user_py)";
cmd_define_ess_py = "    iess_py = define_essential ( mesh_py, problem_py,'curves', [0,1,2,3], degfd=0 );"+...
                    "    iess_py = define_essential ( mesh_py, problem_py,'curves', [0,1,2,3], degfd=1, iessp=iess_py );"+...
                    "    iess_py = define_essential ( mesh_py, problem_py,'points', 0, physq=1, iessp=iess_py  )";

cmd_fill_sys_py =   "    uess_py = fill_system_vector ( mesh_py, problem_py, 'curves', [0,1], func, funcnr=3 );"+...
                    "    uess_py = fill_system_vector ( mesh_py, problem_py, 'curves', [2,3], func, funcnr=3, fin=uess_py )";

cmd_apply_ess_py =  "    A_py2, f_py2, _ = apply_essential ( A_py, f_py, uess_py, iess_py )";

cmd_solve_py     =  "    u_py = spsolve(A_py2.tocsr(), f_py2)";

cmd_deriv_vector =  "    user_py2 = user_py;" + ...
                    "    xr_py = refcoor_nodal_points ( mesh_py );"+...
                    "    user_py2.psi, _ = basis_function('quad','Q1', xr_py );"+...
                    "    user_py2.u = u_py;"+...
                    "    pressure_py = deriv_vector ( mesh_py, problem_py, stokes_pressure, user_py2 )";

cmd_deriv_vector2 =  "    user_py2 = user_py;" + ...
                    "    xr_py = refcoor_nodal_points ( mesh_py );"+... 
                    "    user_py2.phi, user_py2.dphi = basis_function('quad','Q2', xr_py );"+...
                    "    user_py2.u = u_py;"+...
                    "    user_py2.comp = 6;"+...
                    "    divu_py = deriv_vector ( mesh_py, problem_py, stokes_deriv, user_py2 );"+...
                    "    user_py2.comp = 7;"+...
                    "    gammadot_py = deriv_vector ( mesh_py, problem_py, stokes_deriv, user_py2 )";


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
mywritelines("from addons.stokes.stokes_elem import stokes_elem");
mywritelines("from addons.stokes.stokes_pressure import stokes_pressure");
mywritelines("from addons.stokes.stokes_deriv import stokes_deriv");

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
mywritelines("    pressure_ez = Vector()");
write_attrib("    ",pressure_ez,"pressure_ez")
mywritelines("    pressure_ez.vec += -1 # compensate for Python indexing"); 

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
mywritelines("    self.assertTrue(pressure_ez==pressure_py," + ...
    "'deriv_vector failed test, max diff = '+str((abs(pressure_ez.u-pressure_py.u)).max()) )");

mywritelines("  def test_deriv_vector2(self):");
mywritelines("    divu_ez = Vector()");
write_attrib("    ",divu_ez,"divu_ez")
mywritelines("    divu_ez.vec += -1 # compensate for Python indexing");
mywritelines("    gammadot_ez = Vector()");
write_attrib("    ",gammadot_ez,"gammadot_ez")
mywritelines("    gammadot_ez.vec += -1 # compensate for Python indexing"); 

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
mywritelines(cmd_deriv_vector2);
mywritelines("    self.assertTrue(divu_ez==divu_py and gammadot_ez==gammadot_py,"+...
    "'deriv_vector2 failed test, max diff = '+str((abs(divu_ez.u-divu_py.u)).max())"+...
    "+' and '+str((abs(gammadot_ez.u-gammadot_py.u)).max()))");