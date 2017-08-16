#!/usr/bin/env pyton2

import numpy as np
from matplotlib import pyplot as plt

def tline_matrix(theta, Z0):
    """
    ABCD matrix for T-line
    """
    return np.matrix([[np.cos(theta), 1j*Z0*np.sin(theta)], [1j/Z0*np.sin(theta), np.cos(theta)]])

def L_matrix(X):
    """
    ABCD matrix of shunt inductance
    """
    return np.matrix([[1, 0], [1/(1j*X), 1]])

omega_range = np.linspace(0, 2, num = 101)[1:]
gamma_range = np.zeros(omega_range.shape)

for L in [1, 10]:
    for i, omega in enumerate(omega_range):
        theta = omega*2
        X = omega * 1

        m_tline = tline_matrix(theta, 1)
        m_l = L_matrix(X)
        m = m_l*m_tline*m_l*m_tline*m_l*m_tline*m_l*m_tline*m_l
        zin = (m[0,0]*1 + m[0,1])/(m[1,0]*1 + m[1,1])
        gamma = (zin-1)/(zin+1)
        gamma_range[i] = 10*np.log10(np.abs(gamma)**2)

    plt.plot(omega_range, gamma_range, label = 'L = %f' % L)

plt.legend()
plt.show()

