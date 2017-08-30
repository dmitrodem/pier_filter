#!/usr/bin/env python2
# -*- coding: utf-8 -*- 

import numpy as np
from matplotlib import pyplot as plt
from equivalent_circuit import thick_iris_approximaiton
from Utils import toSI as SI
from scipy.constants import c as c0
from WaveguideJunction import WaveguideJunction
from WaveguideElement import WaveguideElement
from GenericElement import GenericElement
from ThickIris import ThickIris

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

def build_tchebysheff_filter(n, ripple, fbw):
    """
    Расчет фильтра Чебышева порядка `n' c рябью в полосе пропускания
    `ripple' и относительной полосой частот `fbw' = (fmax - fmin)/sqrt(fmax*fmin)

    Return: X/Z0 -- нормированные импедансы эквивалентных шунтирующих индуктивностей
            Ang  -- длины резонаторов в радианах (на центральной частоте)
    """
    # прототип фильтра
    g = tchebysheff_prototype(n, ripple)

    # коэффициенты трансформации (Z0 = 1)
    k = -1.0/coupled_microstrip_elements(g, fbw)

    # импедансы индуктивностей на частоте omega = 1
    L = k/(1-k**2)

    #углы
    theta = np.arctan(2*L)

    ang = np.zeros(n)
    for i in xrange(n):
        ang[i] = np.pi - 0.5*(theta[i] + theta[i+1])

    return L, ang

def calculate_filter_response(L, ang, fbw, wc):
    """
    Рассчитываем отклик фильтра
    L -- шунтирующие импедансы, 
    ang -- длины резонирующих секций,
    fbw -- относительная полоса частот,
    wc -- относительная критическая частота (fcutoff/fcenter)

    Return: 
        reflection -- коэффициент отражения
        transmission -- коэффициент прохождения через фильтр
    """
    n = ang.shape[0]
    
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

    return Omega, reflection, transmission

def run():
    print 
    n = 6
    ripple = 0.1
    fbw = 0.1 # TODO: понять, как связана реальная полоса фильтра
              #       с этим параметром

    a = SI("7.112mm")
    h = SI("2.0mm")
    ha = h/a # нормированная толщина диафрагмы
    fcenter = SI("30.1GHz")
    fcutoff = c0/2.0/a

    wc = fcutoff / fcenter # нормированная критическая частота
    
    L, ang = build_tchebysheff_filter(n, ripple, fbw)

    b_interp, theta_interp = thick_iris_approximaiton(wc, ha)

    b       = a*np.array(map(b_interp, L))    # ширины щелей в метрах
    thetas  = np.array(map(theta_interp, L))
    delta_ang = 0.5*(thetas[:-1] + thetas[1:]) # учитываем конечную толщину диафрагм
    angc = ang - delta_ang

    print u"Входные параметры:"
    print u"\tШирина волновода a = %4.2f мм" % (a/SI("1mm"))
    print u"\tТолщина диафрагм h = %4.2f мм" % (h/SI("1mm"))
    print u"\tЦентральная частота Fc = %4.2f ГГц ± %4.2f%%" % (fcenter/SI("1GHz"), 100*fbw)
    print u"\tЧастота среза Fcutoff = %4.2f * Fc" % wc
    print u"\tПорядок фильтра n = %i" % n
    print u"\tНеравномерность в полосе пропускания %4.2f дБ" % ripple
    
    
    print u"Ширины щелей:"
    for i, l in enumerate(L):
        print "\tb[%i] = %4.2f мм" % (i, a*b_interp(l)/SI("1mm"))

    norm_lengths = angc/(2.0*np.pi)/np.sqrt(1-wc**2)
    lengths = c0/fcenter * norm_lengths # длины резонаторов в метрах

    print u"Длины резонирующих секций:"
    for i, length in enumerate(lengths):
        print "\tl[%i] = %4.2f мм" % (i, length/SI("1mm"))
    
    print "Полная длина фильтра = %4.2f мм" % ((np.sum(lengths) + n*h)/SI("1mm"))
    
    nummodes = 10
    minB = np.min(b)

    numModesA = int(a/minB*nummodes)

    frequencies = np.linspace((1-2*fbw) * fcenter, (1+2*fbw) * fcenter, num = 1000)
    s11 = np.zeros(frequencies.shape, dtype = complex)
    
    for i, f in enumerate(frequencies):
        m = WaveguideElement(a, 0, numModesA, f) # затычка нулевой длины
        for j, _ in enumerate(lengths):
            numModesB = int(b[j]/minB*nummodes)
            n1 = ThickIris(a, b[j], h, numModesA, numModesB, f)
            n2 = WaveguideElement(a, lengths[j], numModesA, f)
            m = m*(n1*n2)

        n1 = ThickIris(a, b[-1], h, numModesA, numModesB, f)
        m = m * n1
        s11[i] = m.s21[0,0]

    plt.plot(frequencies/fcenter, 20*np.log10(np.abs(s11)),
             label = 'Exact solution (mode matching)')
    
    Omega, reflection, transmission = calculate_filter_response(L, ang, fbw, wc)
    plt.plot(Omega, 10*np.log10(transmission),
             label = 'Approximate solution (equivalent circuit)')
    plt.grid()
    plt.ylim([-150, 10])
    plt.legend()
    plt.show()

run()
