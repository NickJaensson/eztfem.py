"""Hardcoded regression tests for pos_array and pos_array_vec.

These tests verify the correctness of DN and ND
orderings using a problem with multiple physical quantities.
"""
import unittest
import numpy as np
import eztfem as ezt
from eztfem.core.pos_array import pos_array, pos_array_vec


class TestPosArrayRegression(unittest.TestCase):
    """Hardcoded regression tests for pos_array and pos_array_vec.
    """

    def setUp(self):
        """Set up mesh and problem with multiple physical quantities."""
        self.mesh = ezt.quadrilateral2d([2, 2], 'quad4')
        # Create elementdof with 2 DOFs per node and 2 physical quantities
        # This allows DN and ND orderings to produce different sequences
        elementdof = 2*np.ones((self.mesh.elnumnod, 2), dtype=int)
        self.problem = ezt.Problem(self.mesh, elementdof, nphysq=2)
        self.nodes_array = np.array([0, 1, 2])
        self.nodes_list = [0, 1, 2]
        self.node_single = 0

    def test_pos_array_dn_values(self):
        """Test pos_array DN order produces expected hardcoded values."""
        pos, ndof = pos_array(self.problem, self.nodes_array, order='DN')
        
        # Hardcoded expected values (for 2 physical quantities)
        expected_ndof = np.array([6, 6])
        expected_pos = [[0, 4, 8, 1, 5, 9], [2, 6, 10, 3, 7, 11]]
        
        np.testing.assert_array_equal(ndof, expected_ndof)
        self.assertEqual(pos, expected_pos)

    def test_pos_array_nd_values(self):
        """Test pos_array ND order produces expected hardcoded values."""
        pos, ndof = pos_array(self.problem, self.nodes_array, order='ND')
        
        # Hardcoded expected values (for 2 physical quantities)
        expected_ndof = np.array([6, 6])
        expected_pos = [[0, 1, 4, 5, 8, 9], [2, 3, 6, 7, 10, 11]]
        
        np.testing.assert_array_equal(ndof, expected_ndof)
        self.assertEqual(pos, expected_pos)

    def test_pos_array_list_input_values(self):
        """Test pos_array with list input produces expected hardcoded values."""
        pos, ndof = pos_array(self.problem, self.nodes_list, order='DN')
        
        # Hardcoded expected values (for 2 physical quantities)
        expected_ndof = np.array([6, 6])
        expected_pos = [[0, 4, 8, 1, 5, 9], [2, 6, 10, 3, 7, 11]]
        
        np.testing.assert_array_equal(ndof, expected_ndof)
        self.assertEqual(pos, expected_pos)

    def test_pos_array_single_node_values(self):
        """Test pos_array with single node produces expected hardcoded values."""
        pos, ndof = pos_array(self.problem, self.node_single, order='ND')
        
        # Hardcoded expected values (for 2 physical quantities)
        expected_ndof = np.array([2, 2])
        expected_pos = [[0, 1], [2, 3]]
        
        np.testing.assert_array_equal(ndof, expected_ndof)
        self.assertEqual(pos, expected_pos)

    def test_pos_array_vec_dn_values(self):
        """Test pos_array_vec DN order produces expected hardcoded values."""
        pos, ndof = pos_array_vec(self.problem, self.nodes_array, order='DN')
        
        # Hardcoded expected values (for 2 vectors)
        expected_ndof = np.array([6, 6])
        expected_pos = [[0, 2, 4, 1, 3, 5], [0, 2, 4, 1, 3, 5]]
        
        np.testing.assert_array_equal(ndof, expected_ndof)
        self.assertEqual(pos, expected_pos)

    def test_pos_array_vec_nd_values(self):
        """Test pos_array_vec ND order produces expected hardcoded values."""
        pos, ndof = pos_array_vec(self.problem, self.nodes_array, order='ND')
        
        # Hardcoded expected values (for 2 vectors)
        expected_ndof = np.array([6, 6])
        expected_pos = [[0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]]
        
        np.testing.assert_array_equal(ndof, expected_ndof)
        self.assertEqual(pos, expected_pos)

    def test_pos_array_vec_single_node_values(self):
        """Test pos_array_vec with single node produces expected hardcoded values."""
        pos, ndof = pos_array_vec(self.problem, self.node_single, order='ND')
        
        # Hardcoded expected values (for 2 vectors)
        expected_ndof = np.array([2, 2])
        expected_pos = [[0, 1], [0, 1]]
        
        np.testing.assert_array_equal(ndof, expected_ndof)
        self.assertEqual(pos, expected_pos)


if __name__ == '__main__':
    unittest.main()
