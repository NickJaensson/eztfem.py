import numpy as np

def distribute_elements(n, ratio, factor):
    """
    Generate non-equidistant n elements on the interval [0,1].
    
    Parameters:
    - n: number of elements
    - ratio, factor:
        * ratio = 0: equidistant mesh
        * ratio = 1: the size of the last element is factor times the first
        * ratio = 2: the size of an element is factor times the previous one
        * ratio = 3: the size of the last element is 1/factor times the first
        * ratio = 4: the size of an element is 1/factor times the previous one
          
    Returns:
    - x: coordinates of n+1 points

    The interval [0:1] is divided into n elements:

    dx_{i+1} = g dx_{i-1}

                                                1 - g^n
    1 = (1+g+g^2+g^3+....+g^{n-1}) dx_1     ( = -------- dx1)
                                                 1 - g
    with fac = 1+g+g^2+g^3+....+g^{n-1} we have

    dx_1 = 1 / fac        (=(1-g)/(1-g^n))
    dx_n = g^{n-1} dx_1
    """
    
    if factor < 0:
        raise ValueError("Negative factor")
    
    if n <= 1:
        raise ValueError("Number of elements must be at least 2")

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
        raise ValueError(f"Invalid value for ratio: {ratio}")

    fac = 1

    for i in range(1, n):
        fac += g ** i

    dx_1 = 1 / fac

    x = np.zeros(n + 1)

    x[0] = 0
    dx = dx_1
    for i in range(1, n + 1):
        x[i] = x[i - 1] + dx
        dx *= g

    if abs(x[-1] - 1) > 1e-10:
        raise ValueError(f"End value x(n+1) != 1: {x[-1]:.5e}")
    else:
        x[-1] = 1

    return x
