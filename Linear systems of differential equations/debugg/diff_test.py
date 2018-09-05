import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
import scipy.integrate as int

y0 = [1, 2, 3]
t0 = 2

def f(t, y):
    return [y[0] + 3*y[1] + 6*y[2], y[1] - 3*y[2], y[0] + y[1]]

y = []

y1 = np.array([1+i/1000 for i in range(1000)])

sol = sp.integrate.solve_ivp(f, [2,3], y0, t_eval = [2 + i/1000 for i in range(1000)])
print(sol.t, sol.y)
plt.plot(sol.t, sol.y[0], label = '0')
plt.plot(sol.t, sol.y[1], label = '1')
plt.plot(sol.t, sol.y[2], label = '2')
plt.legend()
#plt.show()

M = np.array([[1,1],[2,3]])
y = np.array([[1],[2]])
print(np.matmul(M,y))
