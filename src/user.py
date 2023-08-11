import numpy as np

# class definition for User objects
class User:

    # initialize attributes (will be filled elsewhere)
    def __init__(self):
        self.wg = np.array([])
        self.xr = np.array([])
        self.phi = np.array([])
        self.dphi = np.array([])

    # equivalence check for testing against Matlab code
    # NOTE: see the matlab file create_pytfem_tests.m for the use of np.squeeze
    def __eq__(self, other):
        check = [np.allclose(self.wg,other.wg,atol=1e-15,rtol=0),
                 np.allclose(np.squeeze(self.xr),other.xr,atol=1e-15,rtol=0),
                 np.allclose(np.squeeze(self.phi),other.phi,atol=1e-15,rtol=0),
                 np.allclose(self.dphi,other.dphi,atol=1e-15,rtol=0)]
        
       # print a warning when not equivalent (for debugging purposes)
        if not all(check):
           print("WARNING: Users not equivalent:")
           print(check)
           
        return all(check)