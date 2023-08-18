import math

def func(nr, x):
    if nr == 1:
        return (x[1] + 1) ** 2
    elif nr == 3:
        return 1 + math.cos(math.pi * x[0]) * math.cos(math.pi * x[1])
    elif nr == 4:
        return 2 * math.pi ** 2 * math.cos(math.pi * x[0]) * math.cos(math.pi * x[1])
    elif nr == 5:
        return 1
    elif nr == 6:
        return math.cos(math.pi * x[0]) * math.cos(math.pi * x[1]) + (x[0] * x[1]) ** 3
    elif nr == 7:
        return 2 * math.pi ** 2 * math.cos(math.pi * x[0]) * math.cos(math.pi * x[1]) - 6 * (x[0] * x[1] ** 3 + x[0] ** 3 * x[1])
    elif nr == 8:
        return math.pi * math.sin(math.pi * x[0]) * math.cos(math.pi * x[1]) - 3 * x[0] ** 2 * x[1] ** 3
    else:
        raise ValueError("Invalid value for 'nr'")
