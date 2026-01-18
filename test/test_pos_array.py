"""Unit tests for pos_array and pos_array_vec functions with different orderings."""
import unittest
import numpy as np
import eztfem as ezt
from eztfem.core.pos_array import pos_array, pos_array_vec
from eztfem.core.problem import Problem


class TestPosArrayOrdering(unittest.TestCase):
    """Test pos_array and pos_array_vec with DN and ND orderings."""

    def setUp(self):
        """Set up a simple mesh and problem for testing."""
        self.mesh = ezt.quadrilateral2d([3, 3], 'quad4')
        elementdof = np.ones((self.mesh.elnumnod, 1), dtype=int)
        self.problem = ezt.Problem(self.mesh, elementdof, nphysq=1)
        self.nodes = np.array([0, 1, 2, 3, 4])

    def test_pos_array_dn_order(self):
        """Test pos_array with DN (dof-then-nodes) ordering."""
        pos, ndof = pos_array(self.problem, self.nodes, order='DN')
        self.assertIsInstance(pos, list)
        self.assertIsInstance(ndof, np.ndarray)
        self.assertEqual(len(pos), self.problem.nphysq)
        self.assertEqual(len(ndof), self.problem.nphysq)
        self.assertGreater(sum(ndof), 0)

    def test_pos_array_nd_order(self):
        """Test pos_array with ND (nodes-then-dof) ordering."""
        pos, ndof = pos_array(self.problem, self.nodes, order='ND')
        self.assertIsInstance(pos, list)
        self.assertIsInstance(ndof, np.ndarray)
        self.assertEqual(len(pos), self.problem.nphysq)
        self.assertEqual(len(ndof), self.problem.nphysq)
        self.assertGreater(sum(ndof), 0)

    def test_pos_array_orderings_consistent(self):
        """Test that DN and ND orderings produce same total DOFs."""
        pos_dn, ndof_dn = pos_array(self.problem, self.nodes, order='DN')
        pos_nd, ndof_nd = pos_array(self.problem, self.nodes, order='ND')

        # Same total DOFs per physical quantity
        np.testing.assert_array_equal(ndof_dn, ndof_nd)

        # Same set of DOF indices (though in different order)
        for i in range(self.problem.nphysq):
            set_dn = set(pos_dn[i])
            set_nd = set(pos_nd[i])
            self.assertEqual(set_dn, set_nd,
                           f"DOF sets differ for physq {i}")

    def test_pos_array_orderings_different_order(self):
        """Test that DN and ND orderings may produce different sequences.
        """
        pos_dn, ndof_dn = pos_array(self.problem, self.nodes, order='DN')
        pos_nd, ndof_nd = pos_array(self.problem, self.nodes, order='ND')

        # Both orderings produce valid results
        self.assertEqual(len(pos_dn), len(pos_nd))
        for i in range(len(pos_dn)):
            self.assertIsInstance(pos_dn[i], list)
            self.assertIsInstance(pos_nd[i], list)

    def test_pos_array_vec_dn_order(self):
        """Test pos_array_vec with DN ordering."""
        pos, ndof = pos_array_vec(self.problem, self.nodes, order='DN')
        self.assertIsInstance(pos, list)
        self.assertIsInstance(ndof, np.ndarray)
        self.assertEqual(len(pos), self.problem.nvec)
        self.assertEqual(len(ndof), self.problem.nvec)

    def test_pos_array_vec_nd_order(self):
        """Test pos_array_vec with ND ordering."""
        pos, ndof = pos_array_vec(self.problem, self.nodes, order='ND')
        self.assertIsInstance(pos, list)
        self.assertIsInstance(ndof, np.ndarray)
        self.assertEqual(len(pos), self.problem.nvec)
        self.assertEqual(len(ndof), self.problem.nvec)

    def test_pos_array_vec_orderings_consistent(self):
        """Test that DN and ND orderings produce same total DOFs for vectors."""
        pos_dn, ndof_dn = pos_array_vec(self.problem, self.nodes, order='DN')
        pos_nd, ndof_nd = pos_array_vec(self.problem, self.nodes, order='ND')

        # Same total DOFs per vector
        np.testing.assert_array_equal(ndof_dn, ndof_nd)

        # Same set of DOF indices (though in different order)
        for i in range(self.problem.nvec):
            set_dn = set(pos_dn[i])
            set_nd = set(pos_nd[i])
            self.assertEqual(set_dn, set_nd,
                           f"DOF sets differ for vector {i}")

    def test_pos_array_vec_orderings_different_order(self):
        """Test that DN and ND orderings may produce different sequences.
        """
        pos_dn, ndof_dn = pos_array_vec(self.problem, self.nodes, order='DN')
        pos_nd, ndof_nd = pos_array_vec(self.problem, self.nodes, order='ND')

        # Both orderings produce valid results
        self.assertEqual(len(pos_dn), len(pos_nd))
        for i in range(len(pos_dn)):
            self.assertIsInstance(pos_dn[i], list)
            self.assertIsInstance(pos_nd[i], list)

    def test_pos_array_single_node(self):
        """Test pos_array with a single node."""
        node = 0
        pos_dn, ndof_dn = pos_array(self.problem, node, order='DN')
        pos_nd, ndof_nd = pos_array(self.problem, node, order='ND')

        np.testing.assert_array_equal(ndof_dn, ndof_nd)

    def test_pos_array_vec_single_node(self):
        """Test pos_array_vec with a single node."""
        node = 0
        pos_dn, ndof_dn = pos_array_vec(self.problem, node, order='DN')
        pos_nd, ndof_nd = pos_array_vec(self.problem, node, order='ND')

        np.testing.assert_array_equal(ndof_dn, ndof_nd)

    def test_pos_array_invalid_order(self):
        """Test pos_array raises error for invalid order."""
        with self.assertRaises(ValueError):
            pos_array(self.problem, self.nodes, order='INVALID')

    def test_pos_array_vec_invalid_order(self):
        """Test pos_array_vec raises error for invalid order."""
        with self.assertRaises(ValueError):
            pos_array_vec(self.problem, self.nodes, order='INVALID')


if __name__ == '__main__':
    unittest.main()
