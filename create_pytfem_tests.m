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

mywritelines("class TestPytfem(unittest.TestCase):");


% test for distribute_elements

cmd = "distribute_elements(8,1,2)";
grid_ez = eval(cmd);

mywritelines("  def test_distribute_elements(self):");
mywritelines("    grid_py = "+cmd);
write1Darr_r("    ",grid_ez,"grid_ez")
mywritelines("    self.assertTrue(np.allclose(grid_py,grid_ez," + ...
    "atol=1e-15,rtol=0),'distribute_elements failed test!' )");


% test for quadrilateral2d

cmd = "quadrilateral2d([3,2],'quad9')";
mesh_ez = eval(cmd);

mywritelines("  def test_quadrilaterial2d(self):");
mywritelines("    mesh_py = "+cmd);

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



% test for problem definition




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

% NOTE: arrays are written "as is" (looping of the first index, than the
% second etc.). It is assumed that in Python they are also read like this

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
                if size(tmp,2)==1
                    write1Darr_i(skip,tmp,name+"."+fns{k});
                elseif ndims(tmp) == 2
                    write2Darr_i(skip,tmp,name+"."+fns{k});
                else
                    write3Darr_i(skip,tmp,name+"."+fns{k});
                end
            else 
                if size(tmp,2)==1
                    write1Darr_r(skip,tmp,name+"."+fns{k}); 
                elseif ndims(tmp) == 2
                    write2Darr_r(skip,tmp,name+"."+fns{k});
                else
                    write3Darr_r(skip,tmp,name+"."+fns{k});
                end
            end
        end
    end

end

