import numpy as np


# class definition for User objects
class User:

    # initialize attributes (will be filled elsewhere)
    def __init__(self):
        self.wg = np.array([])
        self.xr = np.array([])
        self.phi = np.array([])
        self.dphi = np.array([])
        self.psi = np.array([])

    # equivalence check for testing against Matlab code
    # NOTE: see the file NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze
    # NOTE2; only a selection of possible attributes is checked
    def __eq__(self, other):
        # get list of names of attributes
        attributes_self = []
        for attr, _ in self.__dict__.items():
            attributes_self.append(attr)
        attributes_other = []
        for attr, _ in other.__dict__.items():
            attributes_other.append(attr)

        # check if attributes exist in both Users
        check0 = []
        check0.append(('wg' in attributes_self) is ('wg' in attributes_other))
        check0.append(('xr' in attributes_self) is ('xr' in attributes_other))
        check0.append(('phi' in attributes_self)
                      is ('phi' in attributes_other))
        check0.append(('psi' in attributes_self)
                      is ('psi' in attributes_other))
        check0.append(('dphi' in attributes_self)
                      is ('dphi' in attributes_other))
        check0.append(('coorsys' in attributes_self)
                      is ('coorsys' in attributes_other))
        check0.append(('alpha' in attributes_self)
                      is ('alpha' in attributes_other))
        check0.append(('mu' in attributes_self)
                      is ('mu' in attributes_other))
        check0.append(('funcnr' in attributes_self)
                      is ('funcnr' in attributes_other))

        if not all(check0):
            print("WARNING: Users do not have equal attributes:")
            print(check0)
            return False

        # check if existing attributes have same values
        check1 = [True]  # avoid empty list
        if 'wg' in attributes_self:
            check1.append(np.allclose(self.wg, other.wg, atol=1e-15, rtol=0))
        if 'xr' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.xr), other.xr,
                                      atol=1e-15, rtol=0))
        if 'phi' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.phi), other.phi,
                                      atol=1e-15, rtol=0))
        if 'psi' in attributes_self:
            check1.append(np.allclose(np.squeeze(self.psi), other.psi,
                                      atol=1e-15, rtol=0))
        if 'dphi' in attributes_self:
            check1.append(np.allclose(self.dphi, other.dphi,
                                      atol=1e-15, rtol=0))
        if 'coorsys' in attributes_self:
            check1.append(self.coorsys == other.coorsys)
        if 'alpha' in attributes_self:
            check1.append(self.alpha == other.alpha)
        if 'mu' in attributes_self:
            check1.append(self.mu == other.mu)
        if 'funcnr' in attributes_self:
            check1.append(self.funcnr == other.funcnr)

        # print a warning when not equivalent (for debugging purposes)
        if not all(check1):
            print("WARNING: User attributes not equivalent:")
            print(check1)

        return all(check1)
