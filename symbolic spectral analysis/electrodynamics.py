import sympy as sy
import scipy as sp
import numpy as np
from sympy import *

#Returns Fourier amplitudes of 3D vector

def fourier_amplitudes_3(F, r, k, t, w):
    Fx, Fy, Fz = F
    x, y, z = r
    kx, ky, kz = k

    if diff(Fx, x) != 0:
        fx = fourier_transform(Fx, x, kx)
    if diff(fx, y) != 0:
        fx = fourier_transform(fx, y, ky)
    if diff(fx, z) != 0:
        fx = fourier_transform(fx, z, kz)

    if diff(Fy, x) != 0:
        fy = fourier_transform(Fy, x, kx)
    if diff(fy, y) != 0:
        fy = fourier_transform(fy, y, ky)
    if diff(fy, z) != 0:
        fy = fourier_transform(fy, z, kz)

    if diff(Fz, x) != 0:
        fz = fourier_transform(Fz, x, kx)
    if diff(fz, y) != 0:
        fz = fourier_transform(fz, y, ky)
    if diff(fz, z) != 0:
        fz = fourier_transform(fz, z, kz)

    fx = fourier_transform(fx, t, w)
    fy = fourier_transform(fy, t, w)
    fz = fourier_transform(fz, t, w)

    return fx, fy, fz

def inverse_fourier_amplitudes_3(f, r, k, t, w):
    fx, fy, fz = f
    x, y, z = r
    kx, ky, kz = k

    Fx = inverse_fourier_amplitudes_1(fx, r, k, t, w)
    Fy = inverse_fourier_amplitudes_1(fy, r, k, t, w)
    Fz = inverse_fourier_amplitudes_1(fz, r, k, t, w)

    return Fx, Fy, Fz

#Returns Fourier amplitudes of 1D vector

def fourier_amplitudes_1(F, r, k, t, w):
    x, y, z = r
    kx, ky, kz = k

    if diff(F, x) != 0:
        f = fourier_transform(F, x, kx)
    if diff(f, y) != 0:
        f = fourier_transform(f, y, ky)
    if diff(f, z) != 0:
        f = fourier_transform(f, z, kz)

    f = fourier_transform(f, t, w)
    return f

def inverse_fourier_amplitudes_1(f, r, k, t, w):
    x, y, z = r
    kx, ky, kz = k

    if diff(f, kx) != 0:
        F = inverse_fourier_transform(f, kx, x)
    if diff(F, ky) != 0:
        F = inverdr_fourier_transform(F, ky, y)
    if diff(F, kz) != 0:
        F = fourier_transform(F, kz, z)

    F = inverse_fourier_transform(F, w, t)
    return F

#Calculates electric field fourier components from electric potential V

def electric_field_f(V, r, k):
    x, y, z = r
    kx, ky, kz = k
    if diff(V, x) != 0:
        v = fourier_transform(V, x, kx)
    if diff(v, y) != 0:
        v = fourier_transform(v, y, ky)
    if diff(v, z) != 0:
        v = fourier_transform(v, z, kz)

    E = I * v * k
    return E

#Calculates electric potential with time-space charge density distribution;
#e0 = dielectric permitivity of vacuum, c = speed of light in vacuum

def electric_potential_from_rho(rho, r, t, e0, c):

    k = (symbols('kx'), symbols('ky'), symbols('kz'))
    w = symbols('w')

    Rho = fourier_amplitudes_1(rho, r, k, t, w)

    v = Rho/(e0 * (k[0]**2 + k[1]**2 + k[2]**2 - w**2/c**2))

    V = inverse_fourier_amplitudes_1(v, r, k, t, v)

    return V

#calculates electric field from potential

def electric_field_from_potential(V, r):
    Ex = diff(V, x)
    Ey = diff(V, y)
    Ez = diff(V, z)

    return Ex, Ey, Ez

#calculates magnetic field fourier components from magnetic potential A

def magnetic_field_f(A, r, k):
    Ax, Ay, Az = A
    x, y, z = r
    kx, ky, kz = k

    if diff(Ax, x) != 0:
        ax = fourier_transform(Ax, x, kx)
    if diff(ax, y) != 0:
        ax = fourier_transform(ax, y, ky)
    if diff(ax, z) != 0:
        ax = fourier_transform(ax, z, kz)

    if diff(Ay, x) != 0:
        ay = fourier_transform(Ay, x, kx)
    if diff(ay, y) != 0:
        ay = fourier_transform(ay, y, ky)
    if diff(ay, z) != 0:
        ay = fourier_transform(ay, z, kz)

    if diff(Az, x) != 0:
        az = fourier_transform(Az, x, kx)
    if diff(az, y) != 0:
        az = fourier_transform(az, y, ky)
    if diff(az, z) != 0:
        az = fourier_transform(az, z, kz)

    bx = I * (ky*az - kz*ay)
    by = I * (kz*ax - kx*az)
    bz = I * (kx*ay - ky*ax)

    return bx, by, bz

#calculates magnetic potential from time-space current distribution

def magnetic_potential_from_j(j, r, t, e0, c):

    jx, jy, jz = j
    k = (symbols('kx'), symbols('ky'), symbols('kz'))
    w = symbols('w')

    Jx, Jy, Jz = fourier_amplitudes_3(j, r, k, t, w)

    ax = Jx/(e0**2 * c**2 * (k[0]**2 + k[1]**2 + k[2]**2 - w**2/c**2))
    ay = Jy/(e0**2 * c**2 * (k[0]**2 + k[1]**2 + k[2]**2 - w**2/c**2))
    az = Jz/(e0**2 * c**2 * (k[0]**2 + k[1]**2 + k[2]**2 - w**2/c**2))

    A = inverse_fourier_amplitudes_3((ax,ay,az), r, k, t, v)

    return A

#Calculates magnetic field from vector potential

def magnetic_field_from potential(A, r):
    Ax, Ay, Az = A
    x, y, z = r

    Bx = diff(Ay, z) - diff(Az, y)
    By = diff(Az, x) - diff(Ax, z)
    Bz = diff(Ax, y) - diff(Ay, x)

    return Bx, By, Bz
