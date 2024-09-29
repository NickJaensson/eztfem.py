

def stokes_pressure(elem, coor, user, pos):
    """
    Compute the pressure for post-processing in Stokes flow.

    Parameters
    ----------
    elem : int
        Element number.
    coor : ndarray
        Coordinates of the nodes of the element, shape (n_points, n_dim).
    user : User
        User object containing shape function, Gauss points, parameters
    pos : list of ndarray
        Positions of the degrees of freedom of each physical quantity.

    Returns
    -------
    elemvec : ndarray
        The pressure at all nodes.

    Notes
    -----
    This function is intended for use with `deriv_vector` or `plot_sol`.

    """

    elemvec = user.psi @ user.u[pos[1]]

    return elemvec
