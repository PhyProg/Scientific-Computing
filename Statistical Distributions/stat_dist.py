import sympy as sy
import scipy as sp
import numpy as np
import random as rd
import plotly as py
from ipywidgets import interact
from matplotlib import pyplot

x = []
y = []

def plot(x, y, x_min, x_max, y_min, y_max):
    pyplot.ylim(y_min, y_max)
    pyplot.xlim(x_min, x_max)

    pyplot.scatter(x, y, s = 0.01)
    pyplot.show()

def stat(x1, low, high, dx, n):
    hl = high - low
    y = [0 for i in range(int(hl/dx))]
    x = [i*dx - abs(low) for i in range(int(hl/dx))]
    ymax = 0

    for i in range(len(x1)):
        j = int((x1[i]-low)/dx)
        y[j]+=1

    for i in range(int(hl/dx)):
        if y[i]>ymax:
            ymax = float(y[i])

    plot(x, y, low-1, high+1, 0, ymax+5.0)

def gauss(low, high, mu, sigma, n):
    for j in range(n):
        x1 = rd.gauss(mu, sigma)
        if x1>low and x1<high:
            x.append(x1)
        else: j-=1

    return x

def dist_gauss(N, low, high, mu, sigma):
    N = int(input())
    low, high, mu, sigma = map(float, input().split())
    x = gauss(low, high, mu, sigma, N)
    stat(x, low, high, 0.001, N)

def uniform(low, high, n):
    for j in range(n):
        x1 = rd.uniform(low, high)
        if x1 > low and x1 < high:
            x.append(x1)
        else:
            j-=1
    return x

def dist_uniform():
    N = int(input())
    low, high = map(float, input().split())
    x = uniform(low, high, N)
    stat(x, low, high, 0.001, N)

def triangular(low, high, mid, n):
    for j in range(n):
        x.append(rd.triangular(low, high, mid))

    return x

def dist_triangular():
    N = int(input())
    low, high, mid = map(float, input().split())
    x = triangular(low, high, mid, N)
    stat(x, low, high, 0.001, N)

def beta(alfa, bet, n):
    for j in range(n):
        x.append(rd.betavariate(alfa, bet))

    return x

def dist_beta():
    N = int(input())
    alfa, bet = map(float, input().split())
    x = beta(alfa, bet, N)
    stat(x, 0, 1, 0.0001, N)

def gamma(high, alfa, beta, n):
    for j in range(n):
        x1 = rd.gammavariate(alfa, beta)
        if x1<high:
            x.append(x1)
        else: j-=1

    return x

def dist_gamma():
    N = int(input())
    high, alfa, bet = map(float, input().split())
    x = gamma(high, alfa, bet, N)
    stat(x, 0, high, 0.0001, N)

def exponential(low, high, lambd, n):
    for j in range(n):
        x1 = rd.expovariate(lambd)
        if x1>low and x1<high:
            x.append(x1)
        else: j-=1

    return x

def dist_exp():
    N = int(input())
    low, high, lambd = map(float, input().split())
    x = exponential(low, high, lambd, N)
    stat(x, low, high, 0.0001, N)

def lognorm(low, high, mu, sigma, n):
    for j in range(n):
        x1 = rd.lognormvariate(mu, sigma)
        if x1>low and x1<high:
            x.append(x1)
        else: j-=1

    return x

def dist_log():
    N = int(input())
    low, high, mu, sigma = map(float, input().split())
    x = lognorm(low, high, mu, sigma, N)
    stat(x, low, high, 0.0001, N)

def vonmises(mu, kappa, n):
    for j in range(n):
        x1 = rd.vonmisesvariate(mu, kappa)
        x.append(x1)

    return x

def dist_vonmises():
    N = int(input())
    mu, kappa = map(float, input().split())
    x = vonmises(mu, kappa, N)
    stat(x, 0, 6.3, 0.0001, N)

def pareto(low, high, alfa, n):
    for j in range(n):
        x1 = rd.paretovariate(alfa)
        if x1>low and x1<high:
            x.append(x1)
        else: j-=1

    return x

def dist_pareto():
    N = int(input())
    low, high, alfa = map(float, input().split())
    x = pareto(low, high, alfa, N)
    stat(x, low, high, 0.0001, N)