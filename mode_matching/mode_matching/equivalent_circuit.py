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

def run():
    fd = open("/tmp/results.txt", "w")
    bs = np.linspace(0.05, 0.5, num = 100)
    for fc in np.linspace(0.4, 0.9, num = 21):
        print fc
        z, theta = thick_iris_approximaiton(fc, 0.1)
        for b in bs:
            fd.write("{}\t{}\t{}\n".format(fc, b, theta(b)))
        plt.plot(bs, z(bs), label = '$fc/f0 = %f$' % fc)

    fd.close()
    plt.xlabel("$b/a$")
    plt.ylabel("$X/Z_0$")
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == '__main__':
    run()
