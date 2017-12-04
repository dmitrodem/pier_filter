#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from WaveguideElement import WaveguideElement as WG
from WaveguideJunction import WaveguideJunction as WJ
from GenericElement import GenericElement as GE
from ThickIris import ThickIris
from Utils import toSI as SI

import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import c as c0
from scipy.optimize import newton_krylov
from scipy.interpolate import interp1d

def thick_iris_approximaiton(fc, thickness, b_range = (0.05, 0.75), b_num = 100):
    """
    Рассчитываем зависимость нормированного шунтирующего импеданса X/Z0
    от нормированной ширины щели b/a для фиксированного значения 
    частоты среза fc = fcutoff/fcenter

    -[====]--+--[====]--
     theta/2 |  theta/2
             <
             < j*X             Gamma = -exp(-j*theta)/(1+2*j*X)
             <
             |
            ---
    Result: b_interp     -- нормированная на a ширина щели,
            theta_interp -- длина линии
    """
    
    minmodes = 50
    fcenter = SI("1GHz") # от балды -- центральную частоту
    fcutoff = fc*fcenter # частота среза в герцах
 
    a = c0/2.0/fcutoff # ширина секции в метрах
    l = thickness*a    # длина перегородки
    bmin, bmax = b_range
    bs_norm = np.linspace(bmin, bmax, num = b_num)
    g     = np.zeros(bs_norm.shape, dtype = complex)
    t     = np.zeros(bs_norm.shape, dtype = complex)
    x     = np.zeros(bs_norm.shape, dtype = complex)
    theta = np.zeros(bs_norm.shape, dtype = complex)
    
    for i, b_norm in enumerate(bs_norm):
        numModesA = int(minmodes/b_norm)
        numModesB = minmodes
        mp = ThickIris(a, a*b_norm, l, numModesA, numModesB, fcenter)
        g[i] = mp.s11[0,0]
        t[i] = mp.s21[0,0]

    x     = 0.5j*t/g
    theta = np.angle(1.0/(t - g))
    
    b_interp     = interp1d(x.real, bs_norm)
    theta_interp = interp1d(x.real, theta)

    assert np.allclose(np.abs(t - g), 1.0) # |-S11 + S21| = 1 для данного приближения
    assert np.allclose(x.imag, 0.0)
    
    return b_interp, theta_interp

def thick_iris_approximaiton1(fc, thickness, b_range = (0.05, 0.75), b_num = 100):
    """
    Рассчитываем зависимость нормированного шунтирующего импеданса X/Z0
    от нормированной ширины щели b/a для фиксированного значения 
    частоты среза fc = fcutoff/fcenter

    -[====]--+--[====]--
     theta/2 |  theta/2
             <
             < j*X             Gamma = -exp(-j*theta)/(1+2*j*X)
             <
             |
            ---
    Result: b_interp     -- нормированная на a ширина щели,
            theta_interp -- длина линии
    """
    
    minmodes = 50
    fcenter = SI("1GHz") # от балды -- центральную частоту
    fcutoff = fc*fcenter # частота среза в герцах
 
    a = c0/2.0/fcutoff # ширина секции в метрах
    l = thickness*a    # длина перегородки
    bmin, bmax = b_range
    bs_norm = np.linspace(bmin, bmax, num = b_num)
    g     = np.zeros(bs_norm.shape, dtype = complex)
    t     = np.zeros(bs_norm.shape, dtype = complex)
    x     = np.zeros(bs_norm.shape, dtype = complex)
    theta = np.zeros(bs_norm.shape, dtype = complex)
    
    for i, b_norm in enumerate(bs_norm):
        numModesA = int(minmodes/b_norm)
        numModesB = minmodes
        mp = ThickIris(a, a*b_norm, l, numModesA, numModesB, fcenter)
        g[i] = mp.s11[0,0]
        t[i] = mp.s21[0,0]

    x     = 0.5j*t/g
    theta = np.angle(1.0/(t - g))
    
    return bs_norm, x.real, theta

def run():

    import pylab
    from matplotlib import rc
    from cycler import cycler

    monochrome = (cycler('color', ['k']) * cycler('marker', ['']) *
              cycler('linestyle', ['-', '-.', '--', (0, (2, 2)), ':']))

    rc('font',**{'family':'serif'})
    rc('text', usetex=True)
    rc('text.latex',unicode=True)
    rc('text.latex',preamble='\usepackage[utf8]{inputenc}')
    rc('text.latex',preamble='\usepackage[russian]{babel}')
    rc('axes', prop_cycle=monochrome)
    
    fcs = np.arange(0.5, 1.0, 0.1)
    b_num = 14

    z     = np.zeros((fcs.shape[0], b_num))
    theta = np.zeros((fcs.shape[0], b_num))
    
    for i, fc in enumerate(fcs):
        print fc
        bs, z[i, :], theta[i, :] = thick_iris_approximaiton1(fc, 2.0 / 7.112, b_range = (0.4, 0.705), b_num = b_num)

    for i, fc in enumerate(fcs):
        plt.plot(bs, z[i, :], label = ('$fc = %4.1f$' % fc))

    plt.legend()
    plt.grid()
    plt.xlim(0.4, 0.7)
    plt.ylim(0, 1)
    plt.xlabel('$b/a$')
    plt.ylabel('$X/Z_0$')
    plt.gcf().set_size_inches(np.array([8, 10])/2.54)
    plt.tight_layout()
    plt.savefig('../../art/images/impedance.pdf')
    plt.show()

    plt.close()
    for i, fc in enumerate(fcs):
        plt.plot(bs, np.rad2deg(theta[i, :]), label = ('$fc = %4.1f$' % fc))
    plt.grid()
    plt.legend()
    plt.xlim(0.4, 0.7)
    plt.xlabel('$b/a$')
    plt.ylabel('$\\theta$, deg')
    plt.gcf().set_size_inches(np.array([8, 10])/2.54)
    plt.tight_layout()
    plt.savefig('../../art/images/theta.pdf')
    plt.show()
    
if __name__ == '__main__':
    run()
