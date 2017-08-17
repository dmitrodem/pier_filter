#!/usr/bin/env python2
# -*- coding: utf-8 -*- 

import numpy as np
from matplotlib import pyplot as plt

def butterworth_prototype(order):
    """
    butterworth_prototype(order):
    
    Return g[n] prototype values for maximally flat lowpass filter.
    Values are normalized, i.e. Z[load] = 1, omega[c] = 1
    """
    g = np.fromfunction(lambda k: 2*np.sin((2*k-1)*np.pi/(2*order)), (order+2, ))
    g[0] = 1.0
    g[-1] = 1.0
    return g

def tchebysheff_prototype(order, ripple):
    """
    tchebysheff_prototype(order, ripple):
    
    Return g[n] prototype values for Tchebysheff lowpass filter of
    order `order' and maximum ripple in passband `ripple' dB.
    Values are normalized, i.e. Z[load] = 1, omega[c] = 1
    """
    n = order
    g = np.zeros(n+2)
    beta = -np.log(np.tanh(ripple/17.37))
    gamma = np.sinh(beta/(2*n))

    g[0] = 1.0
    g[1] = 2.0/gamma * np.sin(np.pi/(2*n))
    for i in xrange(2, n+1):
        g[i] = 1/g[i-1] * 4*np.sin((2*i-1)*np.pi/(2*n)) * np.sin((2*i-3)*np.pi/(2*n)) / \
        (gamma**2 + np.sin((i-1)*np.pi/n)**2)

    if n & 1 == 1:
        # odd order
        g[n+1] = 1.0
    else:
        # even order
        g[n+1] = 1.0/np.tanh(beta/4)**2

    return g
        

def coupled_microstrip_elements(g, fbw):
    """
    g -- prototype elements,
    fbw -- fractional bandwidth : (f_max - f_min)/f_c
    """
    n = g.shape[0] - 2
    c = np.zeros(n + 1)
    c[0] = np.sqrt(np.pi/2 * fbw / g[0] / g[1])
    for j in xrange(1, n):
        c[j] = np.pi/2 * fbw * 1/np.sqrt(g[j]*g[j+1])
    c[n] = np.sqrt(np.pi/2 * fbw / g[n] / g[n+1])
    return c

def ABCD_tline(theta):
    return np.matrix([[np.cos(theta), 1j*np.sin(theta)], [1j*np.sin(theta), np.cos(theta)]])

def ABCD_series(Z):
    return np.matrix([[1, Z], [0, 1]])

def ABCD_shunt(Y):
    return np.matrix([[1, 0], [Y, 1]])


def run():
    n = 6
    ripple = 0.5
    fbw = 0.052
    wc = 0.5766490461777093 # wc/w0 -- нормированная критическая частота
    
    # прототип фильтра
    g = tchebysheff_prototype(n, ripple)

    # коэффициенты трансформации (Z0 = 1)
    k = -1.0/coupled_microstrip_elements(g, fbw)

    # импедансы индуктивностей на частоте omega = 1
    L = k/(1-k**2)

    print u"индуктивности катушек: ", L

    #углы
    theta = np.arctan(2*L)
    print theta

    ang = np.zeros(n)
    for i in xrange(n):
        ang[i] = np.pi - 0.5*(theta[i] + theta[i+1])

    print u"длины резонаторов в радианах:", ang

    Omega = np.linspace(1 - 2*fbw, 1 + 2*fbw, num = 1001)[1:]
    Gamma = np.zeros(Omega.shape, dtype = complex)

    for index, omega in enumerate(Omega):
        m = np.matrix([[1, 0], [0, 1]], dtype = complex)
        for i in xrange(n):
            m = m*ABCD_shunt(-1j/(L[i] * omega))
            m = m* ABCD_tline(ang[i] * np.sqrt(omega**2-wc**2)/np.sqrt(1.0 - wc**2)) # учитываем дисперсию волновода
        m = m*ABCD_shunt(-1j/(L[n] * omega))
        zin = (m[0,0] + m[0,1])/(m[1,0] + m[1,1])
        gamma = (zin-1)/(zin+1)
        Gamma[index] = gamma

    reflection = np.abs(Gamma)**2
    transmission = 1.0 - reflection
    plt.plot(Omega, 10*np.log10(transmission))
    plt.grid()
    plt.ylim([-120, 10])
    plt.show()

run()
