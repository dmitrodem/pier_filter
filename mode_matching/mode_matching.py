#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
import sys
from scipy.integrate import quad

# assuming a < b

a = np.sqrt(0.5)
b = 1.0
f = 1.0 # f/c0 actually, f_c = 1/(2a)

assert (a < b)

M = 1
N = M

DTYPE = np.complex

def eigenfunction(a, n, x):
    """
    собственная функция f_a^n(x) прямоугольного волновода, нормированная на 1
    учитываем только
    1. TE-волны
    2. зависимость от y отсутствует (m = 0)
    3. только нечетные моды (n = 2 n' + 1)
    """
    n = 2*n + 1
    return np.sqrt(2.0/a) * np.cos(np.pi*n*x/a)

def _V(m, n):
    """
    матричный элемент V_{m, n} = \int_\limits_{-a/2}^{a/2} f_a^m(x) f_b^n(x) dx
    """
    global a
    global b
    _m = m
    m = DTYPE(2*m + 1)
    n = DTYPE(2*n + 1)
    if m/a == n/b:
        result = np.sqrt(a/b)
    else:
        result = 4./np.pi * a**(-1.5) * b**(-0.5) * m * (-1)**(_m) * np.cos(0.5 * np.pi * a/b * n) / ((m/a)**2 - (n/b)**2)
    return result

def _Y(n, a):
    global f
    n = DTYPE(2*n + 1)
    return np.sqrt(1.0 - (0.5*n/(a*f))**2)

def _Ya(m, n):
    global a
    return _Y(n, a)

def _Yb(m, n):
    global b
    return _Y(n, b)

V  = np.vectorize(_V)
Ya = np.vectorize(_Ya)
Yb = np.vectorize(_Yb)

I   = np.identity(M, dtype = DTYPE)
Vmn = np.fromfunction(V, (M, N), dtype = DTYPE)
Vnm = np.transpose(Vmn)

Yamn = np.fromfunction(Ya, (M, N), dtype = DTYPE)
Ybmn = np.fromfunction(Yb, (M, N), dtype = DTYPE)

P = np.zeros((2*M, 2*M), dtype = DTYPE)
P[0:M,   0:M]   = I
P[0:M, M:2*M]   = I
P[M:2*M, 0:M]   = -Vnm * Yamn
P[M:2*M, M:2*M] =  Vnm * Yamn

Q = np.zeros((2*M, 2*M), dtype = DTYPE)
Q[0:M, 0:M]     = Vmn
Q[0:M, M:2*M]   = Vmn
Q[M:2*M, 0:M]   = I*Ybmn
Q[M:2*M, M:2*M] = -I*Ybmn

R = np.dot(lin.inv(P), Q)

def conv_to_S(R):
    """
    convert R matrix to S matrix as follows:
    R = [[X Y] [Y X]]
    S = [[XY' Y-XY'X] [Y' -Y'X]], where Y' = inv(Y)
    """
    global M
    X = R[0:M, 0:M]
    Y = R[0:M, M:2*M]
    Yi = lin.inv(Y)

    result = np.zeros(R.shape)
    result[0:M, 0:M] = np.dot(X, Yi)
    result[0:M, M:2*M] = Y - np.dot(X, np.dot(Yi, X))
    result[M:2*M, 0:M] = Yi
    result[M:2*M, M:2*M] = -np.dot(Yi, X)

    return result

S = conv_to_S(R)

A = np.zeros(2*M)
A[0] = 1.0
B = np.dot(S, A)

xa = np.linspace(-a/2, a/2, num = 100)
Ea = np.zeros(xa.shape, dtype = DTYPE)
for i in xrange(M):
    Ea += (A[i] + B[i]) * eigenfunction(a, i, xa)

xb = np.linspace(-b/2, b/2, num = 100)
Eb = np.zeros(xb.shape, dtype = DTYPE)
for i in xrange(M):
    Eb += (A[M+i] + B[M+i]) * eigenfunction(b, i, xb)



testfunc = lambda x: eigenfunction(a, 11, x) * eigenfunction(b, 17, x)
print "Q = ", quad(testfunc, -a/2, a/2)[0]
print "P = ", _V(11, 17)
