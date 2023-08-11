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
mywritelines("    self.assertTrue(np.allclose(sol_pytfem,sol_eztfem," + ...
    "atol=1e-15,rtol=0),'distribute_elements failed test!' )");


% test for quadrilateral2d

cmd = "quadrilateral2d([3,2],'quad9')";
sol_eztfem = eval(cmd);

mywritelines("  def test_quadrilaterial2d(self):");
mywritelines("    sol_pytfem = "+cmd);

mywritelines("    sol_eztfem = Mesh()");
write_attrib("    ",sol_eztfem,"sol_eztfem")

for i=1:sol_eztfem.ncurves
    mywritelines("    sol_eztfem.curves.append(Geometry())");
    write_attrib("    ",sol_eztfem.curves(i),"sol_eztfem.curves["+string(i-1)+"]")
    mywritelines("    sol_eztfem.curves["+string(i-1)+"].topology = sol_eztfem.curves["...
        +string(i-1)+"].topology - 1 # Python indexing");
    mywritelines("    sol_eztfem.curves["+string(i-1)+"].nodes = sol_eztfem.curves["...
        +string(i-1)+"].nodes - 1 # Python indexing"); 
end

mywritelines("    # compensate for zero-based indexing");
mywritelines("    sol_eztfem.topology = sol_eztfem.topology - 1 # Python indexing");
mywritelines("    sol_eztfem.points = sol_eztfem.points - 1 # Python indexing");

mywritelines("    self.assertTrue(sol_pytfem==sol_eztfem,'quadrilateral2d failed test!' )");



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

