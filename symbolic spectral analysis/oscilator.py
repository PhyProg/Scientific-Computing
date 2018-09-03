import sympy as sy
import scipy as sp
import numpy as np
from sympy import *

#Solves forced oscilator EOM with force F (beautiful enough such that Fourier transform can be applied)

def forced_oscilator(F, t, w, w_0, b):
    Fx, Fy, Fz = F
    x, y, z = Function('x')(t), Function('y')(t), Function('z')(t)

    fx = fourier_transform(Fx, t, w)
    fy = fourier_transform(Fy, t, w)
    fz = fourier_transform(Fz, t, w)

    xw = fx/(w**2 - w_0**2 + 2*I*b*w)
    yw = fy/(w**2 - w_0**2 + 2*I*b*w)
    zw = fz/(w**2 - w_0**2 + 2*I*b*w)

    x = inverse_fourier_transform(xw, w, t)
    y = inverse_fourier_transform(yw, w, t)
    z = inverse_fourier_transform(zw, w, t)

    return x, y, z

def oscilator(x0, v0, w, t, w0, b):
    return x0 * exp(-2bt) * cos(w*t) + v0/w0 * exp(-2bt) * sin(w*t)
