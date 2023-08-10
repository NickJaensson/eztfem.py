import numpy as np
from math import isclose

def distribute_elements(n, ratio, factor):

# DISTRIBUTE_ELEMENTS  Generate noequidistant n elements on the interval [0,1] 
#   x = DISTRIBUTE_ELEMENTS ( n, ratio, factor )
#   input:
#     n: number of elements
#     ratio, factor: 
#       ratio = 0: equidistant mesh
#       ratio = 1: the size of the last element is factor times the first
#       ratio = 2: the size of an element is factor times the previous one
#       ratio = 3: the size of the last element is 1/factor times the first
#       ratio = 4: the size of an element is 1/factor times the previous one
#   output:
#     x : coordinates of n+1 points [numpy array]
# 
#   The interval [0:1] is divided into n elements:
# 
#      dx_{i+1} = g dx_{i-1}
#                                                1 - g^n
#      1 = (1+g+g^2+g^3+....+g^{n-1}) dx_1   (= -------- dx_1)
#                                                1 - g
#    with fac = 1+g+g^2+g^3+....+g^{n-1} we have
# 
#      dx_1 = 1 / fac        (=(1-g)/(1-g^n))
#      dx_n = g^{n-1} dx_1     

    if factor < 0:
        raise ValueError("negative factor")

    if n <= 1:
        raise ValueError("number of elements must be at least 2")

    if ratio == 0:
        g = 1
    elif ratio == 1:
        g = np.exp(np.log(factor) / (n - 1))
    elif ratio == 2:
        g = factor
    elif ratio == 3:
        g = np.exp(-np.log(factor) / (n - 1))
    elif ratio == 4:
        g = 1 / factor
    else:
        raise ValueError("Invalid value for ratio: {}".format(ratio))

    # generate mesh
    fac = 1.0
    for i in range(1, n):
        fac = fac + g ** i

    # size of first element
    dx_1 = 1.0 / fac

    # generate all elements
    x = np.zeros(n+1)
    x[0] = 0
    dx = dx_1
    for i in range(n):
        x[i+1] = x[i] + dx
        dx = g * dx

    # test whether x(n+1) = 1
    if not isclose(x[n],1.0,abs_tol=1e-10):
        raise ValueError("end value x(n+1) /= 1: {:0.5e}".format(x[n]))
    else:
        x[n] = 1.0

    return x