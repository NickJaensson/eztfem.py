import numpy as np

def distribute_elements(n, ratio, factor):
    """
    DISTRIBUTE_ELEMENTS  Generate noequidistant n elements on the interval [0,1] 
      [ x ] = DISTRIBUTE_ELEMENTS ( n, ratio, factor )
      input:
        n: number of elements
        ratio, factor: 
          ratio = 0: equidistant mesh
          ratio = 1: the size of the last element is factor times the first
          ratio = 2: the size of an element is factor times the previous one
          ratio = 3: the size of the last element is 1/factor times the first
          ratio = 4: the size of an element is 1/factor times the previous one
      output:
        x : coordinates of n+1 points
    
      The interval [0:1] is divided into n elements:
    
         dx_{i+1} = g dx_{i-1}
                                                   1 - g^n
         1 = (1+g+g^2+g^3+....+g^{n-1}) dx_1   (= -------- dx_1)
                                                   1 - g
       with fac = 1+g+g^2+g^3+....+g^{n-1} we have
    
         dx_1 = 1 / fac        (=(1-g)/(1-g^n))
         dx_n = g^{n-1} dx_1
         
    """
    
    if factor < 0:
        raise ValueError("Negative factor")
    
    if n <= 1:
        raise ValueError("Number of elements must be at least 2")

    match ratio:
        case 0:
            g = 1
        case 1:
            g = np.exp(np.log(factor) / (n - 1))
        case 2:
            g = factor
        case 3:
            g = np.exp(-np.log(factor) / (n - 1))
        case 4:
            g = 1 / factor
        case _:
            raise ValueError(f"Invalid value for ratio: {ratio}")

#   generate mesh

    fac = 1.0 + sum(g**i for i in range(1, n))

#   size of first element

    dx_1 = 1.0 / fac

#   generate all elements 

    dx_values = [dx_1 * g ** i for i in range(n)]
    
    x = np.zeros(n + 1)
    x[1:] = np.cumsum(np.array(dx_values))

#   test whether x(n+1) = 1

    if abs(x[-1] - 1) > 1e-10:
        raise ValueError(f"End value x(n+1) != 1: {x[-1]:.5e}")
    else:
        x[-1] = 1.0

    return x