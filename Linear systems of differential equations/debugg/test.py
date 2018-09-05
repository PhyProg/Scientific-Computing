import numpy as np
import scipy as sp
import scipy.integrate as int
from scipy.integrate import*

def integrate(M, t):
    x = len(M(0)[0])
    y = len(M(0))
    temp = [0 for tt in t]

    m = np.array([[1.0 for i in range(x)] for i in range(y)])

    for i in range(y):
        for j in range(x):
            mm = [M(tt) for tt in t]

    for i in range(y):
        for j in range(x):
            for tt in range(len(t)):
                temp[tt] = mm[tt][i][j]
            m[i][j] = float(trapz(temp, x = t, dx = t[1]-t[0]))

    return m


M = lambda t: np.array([[1, t, 2*t], [np.exp(t), np.sin(t), t**3], [t**2, 3*t, t**10]])

t = np.linspace(1.234, 10.123, 666)

m = integrate(M, t)

v = np.array([[1],[2],[3]])

print(np.matmul(m, v))
