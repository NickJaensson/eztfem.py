"""Regression tests for a number of fixed bugs.

Each test fails on the code before the corresponding fix.
"""
import copy
import math
import unittest

import numpy as np
from scipy.sparse.linalg import spsolve

import eztfem as ezt


def _taylor_hood(num_el):
    """Q2/Q1 Stokes problem on the unit square."""
    mesh = ezt.quadrilateral2d(num_el, 'quad9')
    elementdof = np.array([[2, 2, 2, 2, 2, 2, 2, 2, 2],
                           [1, 0, 1, 0, 1, 0, 1, 0, 0]], dtype=int).T
    return mesh, ezt.Problem(mesh, elementdof)


def _velocity(x):
    """Divergence-free manufactured velocity field."""
    return np.array([math.sin(math.pi*x[0]) * math.cos(math.pi*x[1]),
                     -math.cos(math.pi*x[0]) * math.sin(math.pi*x[1])])


def _body_force(_nr, x):
    """f = -mu nabla^2 u + grad p, with mu = 1 and p = x*y."""
    return 2 * math.pi**2 * _velocity(x) + np.array([x[1], x[0]])


def _velocity_comp(nr, x):
    return _velocity(x)[nr]


class TestStokesBodyForce(unittest.TestCase):
    """stokes_elem with user.funcnr > 0 (body force)."""

    def _solve(self, num_int_points):
        mesh, problem = _taylor_hood([8, 8])
        curves = [0, 1, 2, 3]

        user = ezt.User()
        user.xr, user.wg = ezt.gauss_legendre(
            'quad', num_int_points=num_int_points)
        user.phi, user.dphi = ezt.basis_function('quad', 'Q2', user.xr)
        user.psi, _ = ezt.basis_function('quad', 'Q1', user.xr)
        user.coorsys = 0
        user.mu = 1
        user.funcnr = 1
        user.func = _body_force

        mat, rhs = ezt.build_system(mesh, problem, ezt.stokes_elem, user)

        iess = ezt.define_essential(mesh, problem, 'curves', curves, degfd=0)
        iess = ezt.define_essential(mesh, problem, 'curves', curves, degfd=1,
                                    iessp=iess)
        iess = ezt.define_essential(mesh, problem, 'points', [0], physq=1,
                                    iessp=iess)
        uess = ezt.fill_system_vector(mesh, problem, 'curves', curves,
                                      _velocity_comp, funcnr=0, degfd=0)
        ezt.fill_system_vector(mesh, problem, 'curves', curves,
                               _velocity_comp, funcnr=1, degfd=1, f=uess)
        ezt.apply_essential(mat, rhs, uess, iess)
        sol = spsolve(mat.tocsr(), rhs)

        nodes = np.arange(mesh.nnodes)
        exact = ezt.fill_system_vector(mesh, problem, 'nodes', nodes,
                                       _velocity_comp, funcnr=0, degfd=0)
        ezt.fill_system_vector(mesh, problem, 'nodes', nodes,
                               _velocity_comp, funcnr=1, degfd=1, f=exact)
        pos, _ = ezt.pos_array(problem, nodes, physq=0)
        return np.abs(sol[pos[0]] - exact[pos[0]]).max()

    def test_manufactured_solution(self):
        self.assertLess(self._solve(3), 1e-3)

    def test_num_int_points_differs_from_num_nodes(self):
        # 16 integration points, 9 velocity nodes
        self.assertLess(self._solve(4), 1e-3)


class TestAxisymmetricStreamfunction(unittest.TestCase):
    """streamfunction_elem with coorsys = 1."""

    def test_pipe_poiseuille(self):
        mesh = ezt.quadrilateral2d([4, 4], 'quad9')
        radius = mesh.coor[:, 1]
        elementdof = np.array([[1]*9, [2]*9], dtype=int).T
        problem = ezt.Problem(mesh, elementdof, nphysq=1)

        user = ezt.User()
        user.xr, user.wg = ezt.gauss_legendre('quad', num_int_points=3)
        user.phi, user.dphi = ezt.basis_function('quad', 'Q2', user.xr)
        user.coorsys = 1
        user.v = np.zeros(2*mesh.nnodes)
        user.v[0::2] = 1 - radius**2  # axial velocity

        mat, rhs = ezt.build_system(mesh, problem, ezt.streamfunction_elem,
                                    user, posvectors=True)

        xr, user.wg = ezt.gauss_legendre('line', num_int_points=3)
        user.phi, user.dphi = ezt.basis_function('line', 'P2', xr)
        for curve in range(4):
            ezt.add_boundary_elements(mesh, problem, rhs,
                                      ezt.streamfunction_natboun_curve, user,
                                      posvectors=True, curve=curve)

        iess = ezt.define_essential(mesh, problem, 'points', [0])
        ezt.apply_essential(mat, rhs, np.zeros(problem.numdegfd), iess)
        psi = spsolve(mat.tocsr(), rhs)

        exact = 2 * math.pi * (radius**2 / 2 - radius**4 / 4)
        self.assertLess(np.abs(psi - exact).max(), 1e-3)


class TestMissingDegreesOfFreedom(unittest.TestCase):
    """Nodes without the requested degree of freedom must be skipped."""

    def setUp(self):
        self.mesh, self.problem = _taylor_hood([2, 2])

    def test_define_essential_pressure_on_curve(self):
        iess = ezt.define_essential(self.mesh, self.problem, 'curves', [0],
                                    physq=1)
        # pressure is only defined in the 3 vertices on the curve
        self.assertEqual(len(iess), 3)

    def test_define_essential_degfd_out_of_range(self):
        iess = ezt.define_essential(self.mesh, self.problem, 'curves', [0],
                                    degfd=2)
        self.assertEqual(len(iess), 0)

    def test_fill_system_vector_pressure_on_curve(self):
        vec = ezt.fill_system_vector(self.mesh, self.problem, 'curves', [0],
                                     lambda nr, x: 1.0, physq=1)
        self.assertEqual(vec.sum(), 3.0)


class TestMiscellaneous(unittest.TestCase):
    """Smaller fixes."""

    def test_integrate_boundary_elements_curve_zero(self):
        mesh, problem = _taylor_hood([2, 2])
        user = ezt.User()
        xr, user.wg = ezt.gauss_legendre('line', num_int_points=3)
        user.phi, user.dphi = ezt.basis_function('line', 'P2', xr)
        user.coorsys = 0
        pos, _ = ezt.pos_array(problem, np.arange(mesh.nnodes), physq=0,
                               order='ND')
        user.u = np.zeros(problem.numdegfd)
        user.u[pos[0][1::2]] = -1.0  # v = -1: outflow through bottom curve
        flowrate = ezt.integrate_boundary_elements(
            mesh, problem, ezt.stokes_flowrate_curve, user, curve=0)
        self.assertAlmostEqual(flowrate, 1.0, places=12)

    def test_integrate_boundary_elements_curve_required(self):
        mesh, problem = _taylor_hood([2, 2])
        with self.assertRaises(ValueError):
            ezt.integrate_boundary_elements(
                mesh, problem, ezt.stokes_flowrate_curve, ezt.User())

    def test_problem_nphysq_too_large(self):
        mesh = ezt.quadrilateral2d([2, 2], 'quad9')
        elementdof = np.ones((9, 2), dtype=int)
        with self.assertRaises(ValueError):
            ezt.Problem(mesh, elementdof, nphysq=5)

    def test_mesh_merge_leaves_input_meshes_unchanged(self):
        mesh1 = ezt.quadrilateral2d([2, 2], 'quad4')
        mesh2 = ezt.quadrilateral2d([2, 2], 'quad4', origin=[1, 0])
        nodes1 = [copy.deepcopy(crv.nodes) for crv in mesh1.curves]
        nodes2 = [copy.deepcopy(crv.nodes) for crv in mesh2.curves]
        topo2 = [copy.deepcopy(crv.topology) for crv in mesh2.curves]

        ezt.mesh_merge(mesh1, mesh2, curves1=[1], curves2=[3],
                       dir_curves2=[-1], deletecurves1=[1])

        for crv, nodes in zip(mesh1.curves, nodes1):
            np.testing.assert_array_equal(crv.nodes, nodes)
        for crv, nodes, topo in zip(mesh2.curves, nodes2, topo2):
            np.testing.assert_array_equal(crv.nodes, nodes)
            np.testing.assert_array_equal(crv.topology, topo)


if __name__ == '__main__':
    unittest.main()
