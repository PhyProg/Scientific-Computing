import numpy as np
import scipy as sp
import scipy.linalg as lin
import scipy.integrate as int
from matplotlib import pyplot as plt

"""
Linear ODE system solver using matrix exponent.
Here arguments are:
y0 - Initial value;
ti, tf, dt - Initial and final time, and time step;
M - Matrix of linear system
"""

def sol_lin_system_ode_c(y0, ti, tf, dt, M):
    t = np.arange(0, tf-ti, dt)
    y = np.array([[] for i in range(len(y0))])
    for tt in t:
        x = np.matmul(lin.expm(M*tt), y0)
        y = np.c_[y, x]

    for i in range(len(t)):
        t[i] += ti

    return y, t

"""
Linear ODE system solver using Runge-Kutta algorithm (set as default in scipy.integrate.solve_ivp()).
Here arguments are:
y0 - Initial value;
ti, tf, dt - Initial and final time, and time step;
M - Matrix of linear system
"""

def lin_system_ivp_c(y0, ti, tf, dt, M):
    def f(y, M):
        return np.matmul(M, y)
    def g(t, y):
        return f(y, M)
    sol = int.solve_ivp(g, [ti, tf], y0, t_eval = np.arange(ti, tf, dt))
    return sol.y, sol.t

"""
Function for integrating matrix M (lambda t: np.array([[],[],...])) through points t (np.array([])) using
trapezoidal rule.
"""

def integrate_matrix_l(M, t):
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
            m[i][j] = float(int.trapz(temp, x = t, dx = t[1]-t[0]))

    return m

"""
Linear ODE system solver using matrix exponent.
Here arguments are:
y0 - Initial value;
ti, tf, dt - Initial and final time, and time step;
M - Matrix of linear system, defined as: lambda t: np.array().
"""

def sol_lin_system_ode_l(y0, ti, tf, dt, M):
    t = np.arange(ti, tf, dt)
    y = np.array([[] for i in range(len(y0))])
    t1 = [ti]
    for tt in t:
        if tt == ti:
            y = np.c_[y, y0]
            continue

        t1.append(tt)
        m = integrate_matrix_l(M, t1)
        x = np.matmul(lin.expm(m), y0)
        y = np.c_[y, x]

    for i in range(len(t)):
        t[i] += ti

    return y, t

"""
Linear ODE system solver using Runge-Kutta algorithm (set as default in scipy.integrate.solve_ivp()).
Here arguments are:
y0 - Initial value;
ti, tf, dt - Initial and final time, and time step;
M - Matrix of linear system, defined as: lambda t: np.array().
"""

def lin_system_ivp_l(y0, ti, tf, dt, M):
    def f(y, t, M):
        return np.matmul(M(t), y)
    def g(t, y):
        return f(y, t, M)
    sol = int.solve_ivp(g, [ti, tf], y0, t_eval = np.arange(ti, tf, dt))
    return sol.y, sol.t

"""
Function for integrating matrix M (np.array([[[]],[[]],...])) through points t (np.array([])) using
trapezoidal rule. Matrix M at each index [i][j] has array of values at the points t
"""

def integrate_matrix_a(M, t):
    x = len(M(0)[0])
    y = len(M(0))
    temp = [0 for tt in t]

    m = np.array([[1.0 for i in range(x)] for i in range(y)])

    for i in range(y):
        for j in range(x):
            for tt in range(len(t)):
                temp[tt] = M[i][j][tt]
            m[i][j] = float(int.trapz(temp, x = t, dx = t[1]-t[0]))

    return m



"""
Linear ODE system solver using matrix exponent.
Here arguments are:
y0 - Initial value;
ti, tf, dt - Initial and final time, and time step;
M - Matrix of linear system, defined as matrix of arrays representing its values at each point.
"""

def sol_lin_system_ode_a(y0, ti, tf, dt, M):
    t = np.arange(ti, tf, dt)
    if len(M[0][0]) != len(t):
        print('Lenght error')
        return 0, 0
        
    y = np.array([[] for i in range(len(y0))])
    t1 = [ti]
    for tt in t:
        if tt == ti:
            y = np.c_[y, y0]
            continue

        t1.append(tt)
        m = integrate_matrix_a(M, t1)
        x = np.matmul(lin.expm(m), y0)
        y = np.c_[y, x]

    for i in range(len(t)):
        t[i] += ti

    return y, t

"""
Linear ODE system solver using Runge-Kutta algorithm (set as default in scipy.integrate.solve_ivp()).
Here arguments are:
y0 - Initial value;
ti, tf, dt - Initial and final time, and time step;
M - Matrix of linear system, defined as matrix of arrays representing its values at each point.
"""

def lin_system_ivp_a(y0, ti, tf, dt, M):
    t = np.arrange(ti, tf, dt)
    if len(M[0][0]) != len(t):
        print('Lenght error')
        return 0, 0

    def f(y, i, M):

        x = len(M(0)[0])
        y = len(M(0))

        for k in range(y):
            for j in range(x):
                mm[k][j] = M[k][j][i]

        return np.matmul(mm, y)

    def g(t, y):
        i = int((t-ti)/dt)
        return f(y, i, M)

    sol = int.solve_ivp(g, [ti, tf], y0, t_eval = np.arange(ti, tf, dt))
    return sol.y, sol.t
