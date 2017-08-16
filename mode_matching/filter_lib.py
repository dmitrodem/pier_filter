#!/usr/bin/env python2

import numpy as np

def butterworth_prototype(order):
    """
    butterworth_prototype(order):
    
    Return g[n] prototype values for maximally flat lowpass filter.
    Values are normalized, i.e. Z[load] = 1, omega[c] = 1
    """
    g = np.fromfunction(lambda k: 2*np.sin((2*k-1)*np.pi/(2*order)), (order+2, ))
    g[0] = 1.0
    g[-1] = 1.0
    return g[1:]

def tchebysheff_prototype(order, ripple):
    """
    tchebysheff_prototype(order, ripple):
    
    Return g[n] prototype values for Tchebysheff lowpass filter of
    order `order' and maximum ripple in passband `ripple' dB.
    Values are normalized, i.e. Z[load] = 1, omega[c] = 1
    """
    beta = np.log(1.0/np.tanh(ripple/17.37))
    gamma = np.sinh(beta/(2*order))
    a = np.fromfunction(lambda k: np.sin((2*k+1)*np.pi/(2*order)), (order, ))
    b = np.fromfunction(lambda k: gamma**2 + np.sin((k+1)*np.pi/order)**2, (order,))
    g = np.zeros(order+1)
    g[0] = 2*a[0]/gamma
    for i in xrange(1, order):
        g[i] = 4*a[i-1]*a[i]/(b[i-1]*g[i-1])
    if order % 2 == 0:
        g[order] = 1.0/np.tanh(beta/4)**2
    else:
        g[order] = 1.0
    return g

def coupled_microstrip_elements(g, delta):
    """
    g -- prototype elements,
    delta -- fractional bandwidth : (f_max - f_min)/f_c
    """
    c = np.zeros(g.shape)
    for i in xrange(1, g.shape[0]-1):
        c[i] = np.pi*delta/(2*np.sqrt(g[i-1]*g[i]))

    c[0] = np.sqrt(np.pi*delta/(2*g[0]))
    c[-1] = np.sqrt(np.pi*delta/(2*g[-2]*g[-1]))
    return c

g = tchebysheff_prototype(3, 0.5)
K = coupled_microstrip_elements(g, 0.1)
X = K/(1.0 - K**2)
print X
