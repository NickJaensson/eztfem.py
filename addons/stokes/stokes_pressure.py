import numpy as np

def stokes_pressure(elem, coor, user, pos):
  """"
  STOKES_PRESSURE  Element routines (deriv) for pressure post-processing
  [ elemvec ] = STOKES_PRESSURE ( elem, coor, user, pos )
  input:
    elem: element number
    coor: coordinates of the nodes of the element
          coor(i,j) with i the point in space and j the direction in space
    user: used for transferring data to the element routine:
            user.psi the basis functions for the pressure
            user.u the solution vector
    pos: cell array of the positions of the degrees of freedom of each physq
  output:
    elemvec: the pressure in all nodes.
  NOTE: this element is for use with deriv_vector or plot_sol.
  """
  elemvec = user.psi @ user.u [ pos[1] ]
  return elemvec