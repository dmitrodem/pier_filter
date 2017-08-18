#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from WaveguideElement import WaveguideElement as WG
from WaveguideJunction import WaveguideJunction as WJ
from GenericElement import GenericElement as GE
from Utils import toSI as SI

import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import c as c0
from scipy.optimize import newton_krylov
from scipy.interpolate import interp1d

def shunt_impedance_approximation(fc):
    """
    Рассчитываем зависимость нормированного шунтирующего импеданса X/Z0
    от нормированной ширины щели b/a для фиксированного значения 
    частоты среза fc = fcutoff/fcenter

    Result: функция, интерполирующая рассчитанные значения
    """
    minmodes = 50
    fcenter = SI("1GHz") # от балды -- центральную частоту
    fcutoff = fc*fcenter # частота среза в герцах
 
    a = c0/2.0/fcutoff # ширина секции в метрах
    bs_norm = np.linspace(0.05, 0.55, num = 100)
    s11 = np.zeros(bs_norm.shape, dtype = complex)

    for i, b_norm in enumerate(bs_norm):
        numModesA = int(minmodes/b_norm)
        numModesB = minmodes
        m1 = WJ(a, a*b_norm, numModesA, numModesB, fcenter)
        m1p = GE(m1.s22, m1.s21, m1.s12, m1.s11)
        m = m1*m1p
        s11[i] = m.s11[0,0]

    x = 0.5j*(1 + 1/s11)
    z = interp1d(bs_norm, x.real)
    zi = interp1d(x.real, bs_norm)
    return z, zi

def run():
    for fc in np.linspace(0.5, 0.9, num = 11):
        print fc
        z, _ = shunt_impedance_approximation(fc)
        bs = np.linspace(0.05, 0.5, num = 100)
        plt.plot(bs, z(bs), label = '$fc/f0 = %f$' % fc)

    plt.xlabel("$b/a$")
    plt.ylabel("$X/Z_0$")
    plt.grid()
    plt.legend()
    plt.show()

if __name__ == '__main__':
    run()
