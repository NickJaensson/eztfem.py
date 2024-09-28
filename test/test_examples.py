# run with: python -m unittest test_examples.py
from math import isclose
import unittest
import sys
sys.path.append('..')

import examples.poisson.poisson1
import examples.poisson.poisson2
import examples.poisson.poisson3
import examples.poisson.poisson4
import examples.stokes.stokes1
import examples.stokes.stokes2
import examples.stokes.streamfunction1

eps = 1e-15  # maximal difference for comparing floats


class TestPytfem(unittest.TestCase):

    # value from Matlab: 4.071378150172222e-07
    def test_poisson1(self):
        result = examples.poisson.poisson1.main()
        self.assertTrue(isclose(result, 4.071378136849546e-07, abs_tol=eps),
                        'poisson1 failed test!')

    # value from Matlab: 0.294830659717740
    def test_poisson2(self):
        result = examples.poisson.poisson2.main()
        self.assertTrue(isclose(result,  0.29483065971774913, abs_tol=eps),
                        'poisson2 failed test!')

    # value from Matlab: 0.149343721153131 
    def test_poisson3(self):
        result = examples.poisson.poisson3.main()
        self.assertTrue(isclose(result, 0.14934372115313188, abs_tol=eps),
                        'poisson3 failed test!')

    # value from Matlab: 6.522070279402215e-07
    def test_poisson4(self):
        result = examples.poisson.poisson4.main()
        self.assertTrue(isclose(result, 6.522070278291991e-07, abs_tol=eps),
                        'poisson4 failed test!')

    # value from Matlab: 25.865650243176205
    def test_stokes1(self):
        result, _, _, _ = examples.stokes.stokes1.main()
        self.assertTrue(isclose(result, 25.865650243176205, abs_tol=eps),
                        'stokes1 failed test!')

    # value from Matlab: 0.083333333333333
    def test_stokes2(self):
        result, _, _, _ = examples.stokes.stokes2.main()
        self.assertTrue(isclose(result, 0.0833333333333339, abs_tol=eps),
                        'stokes2 failed test!')

    # value from Matlab: 0.008522786557203416
    def test_streamfunction1(self):
        _, mesh, problem, u = examples.stokes.stokes1.main()
        result = examples.stokes.streamfunction1.main(mesh, problem, u)
        self.assertTrue(isclose(result, 0.008522786557203416, abs_tol=eps),
                        'streamfunction1 failed test!')