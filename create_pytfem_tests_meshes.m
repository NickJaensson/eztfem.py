close all; clear

eztfempath = "~/Desktop/eztfem/";
addpath(eztfempath);
addpath(append(eztfempath,"src"))
addpath(append(eztfempath,"addons/meshes"))

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
fn = "~/Desktop/eztfem.py/dotest_meshes.py"; 


%% run the problem in eztfem

mesh_ez1 = quadrilateral2d([2,3],'quad4','origin',[1,1],'length',[4,3]);
mesh_ez2 = quadrilateral2d([4,3],'quad4','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);
mesh_ez3 = quadrilateral2d([2,3],'quad9','origin',[1,1],'length',[4,3]);
mesh_ez4 = quadrilateral2d([4,3],'quad9','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);

%% define the same commands for pytfem

cmd_mesh_py1 =       "    mesh_py = quadrilateral2d([2,3],'quad4',origin=np.array([1,1]),length=np.array([4,3]))";
cmd_mesh_py2 =       "    mesh_py = quadrilateral2d([4,3],'quad4',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";
cmd_mesh_py3 =       "    mesh_py = quadrilateral2d([2,3],'quad9',origin=np.array([1,1]),length=np.array([4,3]))";
cmd_mesh_py4 =       "    mesh_py = quadrilateral2d([4,3],'quad9',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";


%% write some header stuff

writelines("# run with: python -m unittest dotest_meshes.py",fn);
mywritelines("import numpy as np")
mywritelines("import unittest");
mywritelines("from src.distribute_elements import distribute_elements");
mywritelines("from src.quadrilateral2d import quadrilateral2d");
mywritelines("from src.mesh_class import Mesh, Geometry");

mywritelines("class TestPytfem(unittest.TestCase):");


%% write the tests

write_test("test1",mesh_ez1,cmd_mesh_py1);
write_test("test2",mesh_ez2,cmd_mesh_py2);
write_test("test3",mesh_ez3,cmd_mesh_py3);
write_test("test4",mesh_ez4,cmd_mesh_py4);

%% test for quadrilateral2d

function write_test(test_name,mesh_ez_in,cmd_mesh_py_in)
    mywritelines("  def "+test_name+"_quadrilaterial2d(self):");
    mywritelines(cmd_mesh_py_in);
    mywritelines("    mesh_ez = Mesh()");
    write_attrib("    ",mesh_ez_in,"mesh_ez")
    for i=1:mesh_ez_in.ncurves
        mywritelines("    mesh_ez.curves.append(Geometry())");
        write_attrib("    ",mesh_ez_in.curves(i),"mesh_ez.curves["+string(i-1)+"]")
        mywritelines("    mesh_ez.curves["+string(i-1)+"].topology = mesh_ez.curves["...
            +string(i-1)+"].topology - 1 # Python indexing");
        mywritelines("    mesh_ez.curves["+string(i-1)+"].nodes = mesh_ez.curves["...
            +string(i-1)+"].nodes - 1 # Python indexing"); 
    end
    mywritelines("    # compensate for zero-based indexing");
    mywritelines("    mesh_ez.topology = mesh_ez.topology - 1 # Python indexing");
    mywritelines("    mesh_ez.points = mesh_ez.points - 1 # Python indexing");
    mywritelines("    self.assertTrue(mesh_py==mesh_ez,'quadrilateral2d failed test!' )");
end


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