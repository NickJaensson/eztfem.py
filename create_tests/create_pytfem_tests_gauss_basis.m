close all; clear

addpath("subs/");

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
fn = "~/Desktop/eztfem.py/dotest_gauss_basis.py"; 


%% write some header stuff

writelines("# run with: python -m unittest dotest.py",fn);
mywritelines("import numpy as np")
mywritelines("import unittest");
mywritelines("from src.gauss_legendre import gauss_legendre");
mywritelines("from src.basis_function import basis_function");

mywritelines("class TestPytfem(unittest.TestCase):");


% write s number of tests for Guass and basis functions

tests = {'quad',2,'Q2';
         'quad',3,'Q2'};

for ii = 1:size(tests,1)

    shape = tests{ii,1};
    n = tests{ii,2};
    intpol = tests{ii,3};

    %% run the problem in eztfem
    
    [xr_ez,wg_ez] = gauss_legendre(shape,'n', n );
    [phi_ez,dphi_ez] = basis_function(shape, intpol, xr_ez );
    
    
    %% define the same commands for pytfem
    
    cmd_gauss_py = "    xr_py, wg_py = gauss_legendre('"+shape+"',n="+string(n)+" )";
    cmd_basis_py = "    phi_py, dphi_py = basis_function('"+shape+"','"+intpol+"', xr_py )";
    
    
    %% test for gauss_legendre
    
    mywritelines("  def test_gauss_legendre"+string(ii)+"(self):");
    write2Darr_r("     ",xr_ez,"xr_ez");
    write1Darr_r("     ",wg_ez,"wg_ez");
    mywritelines(cmd_gauss_py);
    mywritelines("    check1=np.allclose(xr_py,xr_ez,atol=1e-15,rtol=0)")
    mywritelines("    check2=np.allclose(wg_py,wg_ez,atol=1e-15,rtol=0)")
    mywritelines("    self.assertTrue(check1 and check2,'gauss_legendre failed test!' )");
    
    
    %% test for basis_function
    
    mywritelines("  def test_basis_function"+string(ii)+"(self):");
    write2Darr_r("     ",phi_ez,"phi_ez");
    write3Darr_r("     ",dphi_ez,"dphi_ez");
    mywritelines(cmd_gauss_py);
    mywritelines(cmd_basis_py);
    mywritelines("    check1=np.allclose(phi_py,phi_ez,atol=1e-15,rtol=0)")
    mywritelines("    check2=np.allclose(dphi_py,dphi_ez,atol=1e-15,rtol=0)")
    mywritelines("    self.assertTrue(check1 and check2,'basis_functions failed test!' )");

end