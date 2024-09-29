close all; clear

addpath("subs/");

eztfempath = "~/Desktop/eztfem/";
addpath(eztfempath);
addpath(append(eztfempath,"src/core"))
addpath(append(eztfempath,"src/addons/meshes"))

%% Matlab file to generate testing code using eztfem for pytfem 


%% filename for python test file
global fn
fn = "./test_meshes.py"; 


%% run the problem in eztfem

mesh_ez01 = quadrilateral2d([4,7],'quad4','origin',[1,1],'length',[7,5]);
mesh_ez02 = quadrilateral2d([7,5],'quad4','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);
mesh_ez03 = quadrilateral2d([4,7],'quad9','origin',[1,1],'length',[7,5]);
mesh_ez04 = quadrilateral2d([7,5],'quad9','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);
mesh_ez05 = quadrilateral2d([4,7],'quad5','origin',[1,1],'length',[7,5]);
mesh_ez06 = quadrilateral2d([7,5],'quad5','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);
mesh_ez07 = quadrilateral2d([4,7],'tria3','origin',[1,1],'length',[7,5]);
mesh_ez08 = quadrilateral2d([7,5],'tria3','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);
mesh_ez09 = quadrilateral2d([4,7],'tria4','origin',[1,1],'length',[7,5]);
mesh_ez10 = quadrilateral2d([7,5],'tria4','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);
mesh_ez11 = quadrilateral2d([4,7],'tria6','origin',[1,1],'length',[7,5]);
mesh_ez12 = quadrilateral2d([7,5],'tria6','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);
mesh_ez13 = quadrilateral2d([4,7],'tria7','origin',[1,1],'length',[7,5]);
mesh_ez14 = quadrilateral2d([7,5],'tria7','vertices',[1,1;2,2;2,4;1,4],'ratio',[1,2,3,4],'factor',[1.2,1.3,1.4,1.5]);

mesh_ez15 = line1d(5,'line2');
mesh_ez16 = line1d(5,'line2','ratio',1,'factor',1.2,'length',5,'origin',2);
mesh_ez17 = line1d(5,'line3');
mesh_ez18 = line1d(5,'line3','ratio',1,'factor',1.2,'length',5,'origin',2);

mesh_ez19 = l_shape2d([4,5,6,7],'quad4','factor',[10,10,10,10],'origin',[-1,2],'length',[1,2,3,4]);
mesh_ez20 = two_blocks2d([4,5,6],'quad4','factor',[10,10,10],'origin',[-1,2],'length',[1,2,3]);

mesh_tmp1 = quadrilateral2d([2, 3], 'quad4', 'origin', [-2.0, 0.0], 'length', [2.0, 3.0]);
mesh_tmp2 = quadrilateral2d([4, 3], 'quad4', 'origin', [0.0, 0.0], 'length',[4.0, 3.0]);
mesh_ez21 = mesh_merge(mesh_tmp1, mesh_tmp2, 'curves1', 2, 'curves2', -4, 'deletecurves1', 2, 'deletepoints1',1);


%% define the same commands for pytfem

cmd_mesh_py01 = "    mesh_py = ezt.quadrilateral2d([4,7],'quad4',origin=np.array([1,1]),length=np.array([7,5]))";
cmd_mesh_py02 = "    mesh_py = ezt.quadrilateral2d([7,5],'quad4',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";
cmd_mesh_py03 = "    mesh_py = ezt.quadrilateral2d([4,7],'quad9',origin=np.array([1,1]),length=np.array([7,5]))";
cmd_mesh_py04 = "    mesh_py = ezt.quadrilateral2d([7,5],'quad9',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";
cmd_mesh_py05 = "    mesh_py = ezt.quadrilateral2d([4,7],'quad5',origin=np.array([1,1]),length=np.array([7,5]))";
cmd_mesh_py06 = "    mesh_py = ezt.quadrilateral2d([7,5],'quad5',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";
cmd_mesh_py07 = "    mesh_py = ezt.quadrilateral2d([4,7],'tria3',origin=np.array([1,1]),length=np.array([7,5]))";
cmd_mesh_py08 = "    mesh_py = ezt.quadrilateral2d([7,5],'tria3',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";
cmd_mesh_py09 = "    mesh_py = ezt.quadrilateral2d([4,7],'tria4',origin=np.array([1,1]),length=np.array([7,5]))";
cmd_mesh_py10 = "    mesh_py = ezt.quadrilateral2d([7,5],'tria4',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";
cmd_mesh_py11 = "    mesh_py = ezt.quadrilateral2d([4,7],'tria6',origin=np.array([1,1]),length=np.array([7,5]))";
cmd_mesh_py12 = "    mesh_py = ezt.quadrilateral2d([7,5],'tria6',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";
cmd_mesh_py13 = "    mesh_py = ezt.quadrilateral2d([4,7],'tria7',origin=np.array([1,1]),length=np.array([7,5]))";
cmd_mesh_py14 = "    mesh_py = ezt.quadrilateral2d([7,5],'tria7',vertices=np.array([[1,1],[2,2],[2,4],[1,4]]),ratio=np.array([1,2,3,4]),factor=np.array([1.2,1.3,1.4,1.5]))";

cmd_mesh_py15 = "    mesh_py = ezt.line1d(5,'line2')";
cmd_mesh_py16 = "    mesh_py = ezt.line1d(5,'line2',ratio=1,factor=1.2,length=5.0,origin=2.0)";
cmd_mesh_py17 = "    mesh_py = ezt.line1d(5,'line3')";
cmd_mesh_py18 = "    mesh_py = ezt.line1d(5,'line3',ratio=1,factor=1.2,length=5.0,origin=2.0)";

cmd_mesh_py19 = "    mesh_py = ezt.l_shape2d([4,5,6,7],'quad4',factor=[10,10,10,10],origin=[-1,2],length=[1,2,3,4])";
cmd_mesh_py20 = "    mesh_py = ezt.two_blocks2d([4,5,6],'quad4',factor=[10,10,10],origin=[-1,2],length=[1,2,3])";

cmd_mesh_py21 = "    mesh_tmp1 = ezt.quadrilateral2d([2, 3], 'quad4', origin=[-2.0, 0.0], length=[2.0, 3.0]);"+...
                "    mesh_tmp2 = ezt.quadrilateral2d([4, 3], 'quad4', origin=[0.0, 0.0], length=[4.0, 3.0]);"+ ...
                "    mesh_py = ezt.mesh_merge(mesh_tmp1, mesh_tmp2, curves1=[1], curves2=[3], dir_curves2=[-1], deletecurves1=[1], deletepoints1=[0]);";


%% write some header stuff

writelines("# this test was automatically generated using create_pytfem_tests_meshes.m",fn);
mywritelines("# run with: python -m unittest test_meshes.py");
mywritelines("import numpy as np")
mywritelines("import unittest");
mywritelines("import sys");
mywritelines("sys.path.append('..')");
mywritelines("import eztfem as ezt");

mywritelines("class TestPytfem(unittest.TestCase):");


%% write the tests

write_test("test01",mesh_ez01,cmd_mesh_py01);
write_test("test02",mesh_ez02,cmd_mesh_py02);
write_test("test03",mesh_ez03,cmd_mesh_py03);
write_test("test04",mesh_ez04,cmd_mesh_py04);
write_test("test05",mesh_ez05,cmd_mesh_py05);
write_test("test06",mesh_ez06,cmd_mesh_py06);
write_test("test07",mesh_ez07,cmd_mesh_py07);
write_test("test08",mesh_ez08,cmd_mesh_py08);
write_test("test09",mesh_ez09,cmd_mesh_py09);
write_test("test10",mesh_ez10,cmd_mesh_py10);
write_test("test11",mesh_ez11,cmd_mesh_py11);
write_test("test12",mesh_ez12,cmd_mesh_py12);
write_test("test13",mesh_ez13,cmd_mesh_py13);
write_test("test14",mesh_ez14,cmd_mesh_py14);

write_test("test15",mesh_ez15,cmd_mesh_py15);
write_test("test16",mesh_ez16,cmd_mesh_py16);
write_test("test17",mesh_ez17,cmd_mesh_py17);
write_test("test18",mesh_ez18,cmd_mesh_py18);

write_test("test19",mesh_ez19,cmd_mesh_py19);
write_test("test20",mesh_ez20,cmd_mesh_py20);

write_test("test21",mesh_ez21,cmd_mesh_py21);


%% test for quadrilateral2d

function write_test(test_name,mesh_ez_in,cmd_mesh_py_in)
    mywritelines("  def "+test_name+"_quadrilaterial2d(self):");
    mywritelines(cmd_mesh_py_in);
    mywritelines("    mesh_ez = ezt.Mesh()");
    write_attrib("    ",mesh_ez_in,"mesh_ez")
    for i=1:mesh_ez_in.ncurves
        mywritelines("    mesh_ez.curves.append(ezt.Geometry())");
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