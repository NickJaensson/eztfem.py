# run with: python -m unittest test_examples.py
from math import isclose
import unittest
import sys

# ensure test can be run from eztfem.py/test folder
sys.path.append('../examples/poisson')
sys.path.append('../examples/stokes')
# sys.path.append('../examples/project2')

# ensure test can be run from eztfem.py folder
sys.path.append('examples/poisson')
sys.path.append('examples/stokes')
# sys.path.append('examples/project2')

import poisson1, poisson2, poisson3, poisson4  # type: ignore
import stokes1, stokes2, streamfunction1  # type: ignore
# import gn_channel  # type: ignore

eps = 1e-15  # maximal difference for comparing floats


class TestPytfem(unittest.TestCase):

    # value from eztfem.m (Matlab): 4.071378150172222e-07
    def test_poisson1(self):
        result = poisson1.main()
        self.assertTrue(isclose(result, 4.071378136849546e-07, abs_tol=eps),
                        'poisson1 failed test!')

    # value from eztfem.m (Matlab): 0.294830659717740
    def test_poisson2(self):
        result = poisson2.main()
        self.assertTrue(isclose(result,  0.29483065971774913, abs_tol=eps),
                        'poisson2 failed test!')

    # value from eztfem.m (Matlab): 0.149343721153131 
    def test_poisson3(self):
        result = poisson3.main()
        self.assertTrue(isclose(result, 0.14934372115313188, abs_tol=eps),
                        'poisson3 failed test!')

    # value from eztfem.m (Matlab): 6.522070279402215e-07
    def test_poisson4(self):
        result = poisson4.main()
        self.assertTrue(isclose(result, 6.522070278291991e-07, abs_tol=eps),
                        'poisson4 failed test!')

    # value from eztfem.m (Matlab): 25.865650243176205
    def test_stokes1(self):
        result, _, _, _ = stokes1.main()
        self.assertTrue(isclose(result, 25.865650243176205, abs_tol=eps),
                        'stokes1 failed test!')

    # value from eztfem.m (Matlab): 0.083333333333333
    def test_stokes2(self):
        result, _, _, _ = stokes2.main()
        self.assertTrue(isclose(result, 0.0833333333333339, abs_tol=eps),
                        'stokes2 failed test!')

    # value from eztfem.m (Matlab): 0.008522786557203416
    def test_streamfunction1(self):
        _, mesh, problem, u = stokes1.main()
        result = streamfunction1.main(mesh, problem, u)
        self.assertTrue(isclose(result, 0.008522786557203416, abs_tol=eps),
                        'streamfunction1 failed test!')
        
    # # value from eztfem.m (Matlab): 0.031248765345553
    # def test_gn_channel(self):
    #     result, _, _, _ = gn_channel.main()
    #     self.assertTrue(isclose(result, 0.031248765345550343, abs_tol=eps),
    #                     'gn_channel failed test!')