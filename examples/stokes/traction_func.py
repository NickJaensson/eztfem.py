import numpy as np
import math

def traction_func(nr, x):

    y = np.zeros(2)
     
    if nr == 1:
      y[0] = 1.0
      y[1] = 0.0
    else:
        raise ValueError("Invalid value for 'nr'")
    
    return y
    