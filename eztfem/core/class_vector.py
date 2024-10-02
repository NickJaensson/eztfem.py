import numpy as np


class Vector:
    """
    A class used to represent a Vector.

    Attributes
    ----------
    vec : int
        The vector value.
    u : np.ndarray
        An array to store vector data.

    Methods
    -------
    __init__(self, vec=0):
        Initializes the Vector object with the given attributes.
    __eq__(self, other):
        Checks equivalence of two Vector objects (overloads == sign).

    """

    def __init__(self, vec=0):
        """
        Initializes the Vector object with the given attributes.

        Parameters
        ----------
        vec : int or float, optional
            The vector value (default is 0).

        """
        self.vec = vec
        self.u = np.array([])

    def __eq__(self, other):
        """
        Checks equivalence of two Vector objects (overloads == sign).

        Parameters
        ----------
        other : Vector
            The other Vector object to compare with.

        Returns
        -------
        bool
            True if the Vector objects are equivalent, False otherwise.

        Notes
        -----
        See NOTE_ON_COMPARING_ARRAYS.md for the use of np.squeeze.

        """
        check = [self.vec == other.vec,
                 np.allclose(np.squeeze(self.u), other.u, atol=1e-12, rtol=0)]
        return all(check)
