import numpy as np
from src.isoparametric_deformation import isoparametric_deformation

def poisson_deriv(elem, coor, user, pos):
  """"
  POISSON_DERIV  Element routines (deriv) for post-processing the velocity 
    [ elemvec ] = POISSON_DERIV ( elem, coor, user, pos )
    input:
      elem: element number
      coor: coordinates of the nodes of the element
            coor(i,j) with i the point in space and j the direction in space
      user: used for transferring data to the element routine:
              user.phi the basis functions
              user.dphi the derivatives of basis functions
              user.coorsys =0 for Cartesian, =1 for axisymmetric 
              user.u the solution vector
            NOTE: all must have a value.
      pos: cell array of the positions of the degrees of freedom of each physq
    output:
      elemvec: the gradient (dudx,dudy) in all nodes.
    NOTE: this element is for use with deriv_vector or plot_sol.
  """
  # Set some values
  ndim = coor.shape[1]
  nodalp = user.phi.shape[0]
  ndf = user.phi.shape[1]
  # Compute mapping of reference to real element
  F, Finv, detF = isoparametric_deformation(coor, user.dphi)
  # Compute derivative of the basis functions with respect to the real coordinates
  dphidx = np.zeros((nodalp, ndf, ndim))
  ndr = Finv.shape[1]
  for ip in range(nodalp):
      dphidx[ip, :, :] = user.dphi[ip, :, :].dot(Finv[ip, :, :])
  # compute gradient vector
  gradu = np.zeros((nodalp,ndim))
  for j in range(ndim):
    gradu[:,j] = dphidx[:,:,j] @ user.u[pos[0]]

  elemvec = np.reshape(gradu, [ndim * nodalp],order='F')

  return elemvec