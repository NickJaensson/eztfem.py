import numpy as np

def gauss_legendre(shape, **kwargs):
    """
    Gauss-Legendre integration points and weights
    
    Args:
        shape: shape of the domain:
            'line', interval [-1,1]
            'quad', domain [-1,1]x[-1,1]
            'triangle', reference triangle, left-lower half of [0,1]x[0,1]
        **kwargs: Optional keyword arguments to set optional parameters:
                  'n': the number of integration points in one direction. This only applies
                       to shape='line' and shape='quad', where the number of integration
                       points will be n and n^2, respectively.
                  'p': the order of the integration rule (=order polynomial integrated exact).
                       This only applies to shape='triangle'.
    Returns:
        x, w: coordinates and weights of the integration scheme.
    """
    n = kwargs.get('n', -1)
    p = kwargs.get('p', -1)

    if shape == 'line':
        if n < 0:
            raise ValueError('n must be specified for the integration rule')
        x, w = gauss_legendre_line(n)
    elif shape == 'quad':
        if n < 0:
            raise ValueError('n must be specified for the integration rule')
        x, w = gauss_legendre_quad(n)
    elif shape == 'triangle':
        if p < 0:
            raise ValueError('p must be specified for the integration rule')
        x, w = gauss_legendre_triangle(p)
    else:
        if isinstance(shape, str):
            ch = f'Invalid shape: {shape}'
        else:
            ch = 'shape must be a string'
        raise ValueError(ch)
    
    return x, w


def gauss_legendre_line(n):
    """
    Gauss Legendre in 1D defined on the interval [-1,1]
    
    Args:
        n: number of integration points
    Returns:
        x, w: the reference coordinates and weights of the integration points
    """
    x = np.zeros(n)
    w = np.zeros(n)

    if n == 1:
        x[0] = 0
        w[0] = 2

    elif n == 2:
        x[0] = 2 * 0.2113248654051871177454256097490212 - 1
        x[1] = -x[0]
        w[0:2] = 1

    elif n == 3:
        x[0] = 2 * 0.1127016653792583114820734600217600 - 1
        x[1] = 0
        x[2] = -x[0]
        w[0:2] = [5/9, 8/9]
        w[2] = w[0]

    else:
        raise ValueError(f"Invalid number of points, n = {n}")

    return x, w

def gauss_legendre_quad(n):
    """
    Gauss Legendre in 2D defined on the region [-1,1]x[-1,1]
    
    Parameters:
    n (int): number of integration points in one dimension
    
    Returns:
    x (numpy.ndarray): the reference coordinates of the integration points (shape: (n^2, 2))
    w (numpy.ndarray): the weights of the integration points (shape: (n^2,))
    """
    x1, w1 = gauss_legendre_line(n)

    m = n**2
    x = np.zeros((m, 2))
    w = np.zeros(m)

    for i in range(n):
        for j in range(n):
            ip = i + n * (j)
            x[ip, 0] = x1[i]
            x[ip, 1] = x1[j]
            w[ip] = w1[i] * w1[j]

    return x, w

def gauss_legendre_triangle(n):
    # Implementation of gauss_legendre_triangle
    raise NotImplementedError