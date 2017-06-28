#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from pylab import *
from scipy.constants import c as c0
from pint import UnitRegistry
from numpy import linalg as LA

def eigenfunction(n, l, x):
    return sqrt(2.0/l) * cos(pi*(2*n+1)*x/a)

def Imn(m, n, a, b):
    return 4*(-1)**(n)*a**1.5*b**0.5*cos(0.5*pi*b*(2*m+1)/a)*(2*n+1)/pi/(-b**2*(2*m+1)**2+a**2*(2*n+1)**2)

def gamma(n, l, kv):
    return pi*sqrt(((2*n+1)/l)**2 - kv**2)

def __compute_QP(a, b, Mmax, Nmax, kv):
    Q = fromfunction(lambda m, n: Imn(m, n, a, b), (Mmax, Nmax))
    P = fromfunction(lambda m, n: gamma(n, a, kv)/gamma(m, b, kv)*Imn(n, m, a, b), (Nmax, Mmax))
    return (Q, P)

def __compute_S(a, b, Mmax, Nmax, kv):
    Q, P = __compute_QP(a, b, Mmax, Nmax, kv)
    QP = dot(Q, P)
    PQ = dot(P, Q)
    S11 = dot(inv(eye(Mmax) + QP), -eye(Mmax) + QP)
    S12 = 2 * dot(inv(eye(Mmax) + QP), Q)
    S21 = 2 * dot(inv(eye(Nmax) + PQ), P)
    S22 = dot(inv(eye(Nmax) + PQ), eye(Nmax) - PQ)

    return (S11, S12, S21, S22)

def __adjust_units(a, b, freq, length_unit = 1.0e-3, freq_unit = 1e9, dtype = complex):
    # by default, lengths are in mm, freqs are in GHz
    a = a * length_unit
    b = b * length_unit
    freq = freq * freq_unit
    kv = 2 * freq / c0

    a = dtype(a)
    b = dtype(b)
    kv = dtype(kv)

    return (a, b, kv)

def compute_S(a, b, Mmax, Nmax, freq, **kwargs):
    a, b, kv = __adjust_units(a, b, freq, **kwargs)
    if abs(a) >= abs(b):
        (S11, S12, S21, S22) = __compute_S(a, b, Mmax, Nmax, kv)
    else:
        (S22, S21, S12, S11) = __compute_S(b, a, Nmax, Mmax, kv)

    return (S11, S12, S21, S22)

def __compute_T(a, b, Mmax, Nmax, kv):
    Q, P = __compute_QP(a, b, Mmax, Nmax, kv)
    T11 = 0.5 * (-P + pinv(Q)) # = T22
    T12 = 0.5 * ( P + pinv(Q)) # = T21

    return (T11, T12, T12, T11)

def compute_T(a, b, Mmax, Nmax, freq, **kwargs):
    a, b, kv = __adjust_units(a, b, freq, **kwargs)
    if abs(a) >= abs(b):
        (T11, T12, T21, T22) = __compute_T(a, b, Mmax, Nmax, kv)
    else:
        (T22, T21, T12, T11) = __compute_T(a, b, Mmax, Nmax, kv)

    return (T11, T12, T21, T22)

def compute_gamma(n, l, freq, **kwargs):
    l, _, kv = __adjust_units(l, 0, freq, **kwargs)
    return gamma(n, l, kv)

def convert_T_S(T11, T12, T21, T22):
    S11 = -dot(pinv(T12), T11)
    S12 = inv(T12)
    S21 = T21 - dot(dot(T22, pinv(T12)), T11)
    S22 = dot(T22, pinv(T12))

    return (S11, S12, S21, S22)

def cascade_S(A, B):
    """
    A = (A11, A12, A21, A22) -- S-матрица первого каскада,
    B = (B11, B12, B21, B22) -- S-матрица второго каскада
    Возвращается матрица кортеж C, соответствующий последовательному
    включению каскадов
    """ 
    A11, A12, A21, A22 = A
    B11, B12, B21, B22 = B

    assert(A22.shape == B11.shape)
    N = A22.shape[0]
    
    X = inv(eye(N) - dot(B11, A22))
    Y = inv(eye(N) - dot(A22, B11))

    C11 = A11 + dot(A12, dot(X, dot(B11, A21)))
    C12 = dot(A12, dot(X, B12))
    C21 = dot(B21, dot(Y, A21))
    C22 = B22 + dot(B21, dot(Y, dot(A22, B12)))

    return (C11, C12, C21, C22)

def db(x):
    return 20*log10(abs(x))

def phase(x):
    return rad2deg(angle(x))

def waveguide_S(a, l, M, freq, **kwargs):
    a, l, kv = __adjust_units(a, l, freq, **kwargs)
    gamma_vect = fromfunction(lambda m: gamma(m, a, kv), (M, ))
    S12 = diag(exp(-gamma_vect*l))
    S21 = diag(exp(-gamma_vect*l))
    S11 = zeros((M, M))
    S22 = zeros((M, M))
    return (S11, S12, S21, S22)

def verify_diaphragm():
    a = 1.0001e3 # mm
    b = 0.5001e3 # mm
    l = 20.0     # mm
    f = 0.449688687   # GHz

    M = 100 # in "a" section
    N = 50 # in "b" section
    
    (S11, S12, S21, S22) = compute_S(a, b, M, N, f, dtype = complex128)
    print "S11 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(S11[0,0])), rad2deg(angle(S11[0, 0])))
    print "S11 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(S11[1,1])), rad2deg(angle(S11[1, 1])))
    print "S11 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(S11[2,2])), rad2deg(angle(S11[2, 2])))
    print "S11 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(S11[2,0])), rad2deg(angle(S11[2, 0])))
    print "S22 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(S22[0,0])), rad2deg(angle(S22[0, 0])))

    norm11 = sqrt(fromfunction(lambda m, n: compute_gamma(m, a, f)/compute_gamma(n, a, f), (M, M)))
    norm12 = sqrt(fromfunction(lambda m, n: compute_gamma(m, a, f)/compute_gamma(n, b, f), (M, N)))
    norm21 = sqrt(fromfunction(lambda m, n: compute_gamma(m, b, f)/compute_gamma(n, a, f), (N, M)))
    norm22 = sqrt(fromfunction(lambda m, n: compute_gamma(m, b, f)/compute_gamma(n, b, f), (N, N)))    

    print S12.shape, norm12.shape
    
    normS11 = S11*norm11
    normS12 = S12*norm12
    normS21 = S21*norm21
    normS22 = S22*norm22

    print "S11 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS11[2,0])), rad2deg(angle(normS11[2, 0])))
    print "S12 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS12[0,0])), rad2deg(angle(S12[0, 0])))
    print "S12 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS12[0,2])), rad2deg(angle(S12[0, 2])))
    print "S12 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS12[2,0])), rad2deg(angle(S12[2, 0])))

    print
    print "S21 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS21[0,0])), rad2deg(angle(S21[0, 0])))
    print "S21 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS21[0,2])), rad2deg(angle(S21[0, 2])))
    print "S21 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS21[2,0])), rad2deg(angle(S21[2, 0])))

    print
    print "S22 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS22[0,0])), rad2deg(angle(S22[0, 0])))
    print "S22 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS22[0,2])), rad2deg(angle(S22[0, 2])))
    print "S22 = %4.2f dB, angle = %4.2f deg" % (20*log10(abs(normS22[2,0])), rad2deg(angle(S22[2, 0])))

    freq = linspace(0.25, 1, num = 201)
    s11 = zeros(freq.shape[0])

    for i in xrange(freq.shape[0]):
        f = freq[i]
        A = compute_S(a, b, M, N, f, dtype = complex)
        B = waveguide_S(b, l, N, f, dtype = complex)
        C = compute_S(b, a, N, M, f, dtype = complex)
        D = D11, D12, D21, D22 = cascade_S(cascade_S(A, B), C)
        s11[i] = db(D[0][0,0])

    #print db(D11[0,0]), phase(D11[0,0])

    plot(freq, s11)
    grid()

    d = loadtxt("diaphragm.csv", skiprows=1, delimiter=',')
    ref_freq = d[:, 0]
    ref_s11  = d[:, 1]
    plot(ref_freq, ref_s11)
    show()

def load_constants(filename = "vars.csv"):
    data = loadtxt(filename, delimiter = '\t', dtype = str)
    ureg = UnitRegistry(system = 'mks')
    result = {}
    for (k, v) in data:
        v = ureg.parse_expression(v)
        try:
            v.ito(ureg.millimeter)
            v = v.magnitude
        except:
            pass
        result.update({k : v})
    return result

def eval_filter6_30_a(filename = "vars.csv"):
    
    freq = 30.5
    d = load_constants(filename)

    freqs = linspace(26.0, 36.0, num = 1001)
    s11 = zeros(freqs.shape)


    for M in [5, 10, 20]:
        print M
        def getM(val):
            return int(M * d[val]/d['a'])
    
        for i, freq in enumerate(freqs):
            res = waveguide_S(d['a1'], d['ddd'], getM('a1'), freq, dtype = complex)
            res = cascade_S(res,
                            compute_S(d['a1'], d['c4'], getM('a1'), getM('c4'), freq, dtype = complex))
            res = cascade_S(res,
                            waveguide_S(d['c4'], d['dd'], getM('c4'), freq, dtype = complex))
            res = cascade_S(res,
                            compute_S(d['c4'], d['a'], getM('c4'), getM('a'), freq, dtype = complex))
            res = cascade_S(res,
                            waveguide_S(d['a'], d['d3'], getM('a'), freq, dtype = complex))
            res = cascade_S(res,
                             compute_S(d['a'], d['c3'], getM('a'), getM('c3'), freq, dtype = complex))
            res = cascade_S(res,
                             waveguide_S(d['c3'], d['dd'], getM('c3'), freq, dtype = complex))
            res = cascade_S(res,
                             compute_S(d['c3'], d['a'], getM('c3'), getM('a'), freq, dtype = complex))
            res = cascade_S(res,
                            waveguide_S(d['a'], d['d2'], getM('a'), freq, dtype = complex))
            res = cascade_S(res,
                            compute_S(d['a'], d['c2'], getM('a'), getM('c2'), freq, dtype = complex))
            res = cascade_S(res,
                            waveguide_S(d['c2'], d['dd'], getM('c2'), freq, dtype = complex))
            res = cascade_S(res,
                            compute_S(d['c2'], d['a'], getM('c2'), getM('a'), freq, dtype = complex))
            res = cascade_S(res,
                            waveguide_S(d['a'], d['d1'], getM('a'), freq, dtype = complex))
            res = cascade_S(res,
                            compute_S(d['a'], d['c1'], getM('a'), getM('c1'), freq, dtype = complex))
            res = cascade_S(res,
                            waveguide_S(d['c1'], 0.5*d['dd'], getM('c1'), freq, dtype = complex))

            res1 = res[3], res[2], res[1], res[0]

            res = cascade_S(res, res1)
            
            s11[i] = db(res[1][0,0])

        plot(freqs, s11, label = 'M = %i' % M)
        xlabel("Frequency [GHz]")
        ylabel("$S_{21}$ [dB]")
        title("Transmission coefficient (a.k.a. $S_{21}[0,0]$)")

    grid(True)
    legend()
    show()

def check_self_inverse(S):
    """
    check_self_inverse(S): 
    проверяем, выполняется ли соотношение S*S^-1 = E
    """
    s11, s12, s21, s22 = S
    M = s11.shape[0]
    N = s22.shape[0]

    X = zeros((M+N, M+N), dtype = s11.dtype)
    E = eye(M+N, dtype = s11.dtype)
    X[0:M, 0:M] = s11
    X[0:M, M:M+N] = s12
    X[M:M+N, 0:M] = s21
    X[M:M+N, M:M+N] = s22

    return allclose(E, dot(X, X))

def test_waveguide_S():
    M = 1
    d = load_constants("vars.csv")
    res = waveguide_S(8.636, 44.334, M, 40.0, dtype = complex)
    Sa = zeros((2, 2), dtype = complex)
    Sa[0, 0] = res[0][0,0]
    Sa[0, 1] = res[1][0,0]
    Sa[1, 0] = res[2][0,0]
    Sa[1, 1] = res[3][0,0]
    Sr = zeros((2, 2), dtype = complex)
    Sr[0, 0] = -3.507559781674718E-04 + 6.236349064019116E-04j
    Sr[0, 1] = -4.681279540461796E-01 - 8.836604023572072E-01j
    Sr[1, 0] = -4.681279540443931E-01 - 8.836604023538351E-01j
    Sr[1, 1] = -7.129774651205567E-04 - 6.011145426055048E-05j
    passed = abs(phase(Sr[0, 1]) - phase(Sa[0, 1])) < 1
    msg = "S12 : Reference phase = %4.2f deg, calculated phase = %4.2f deg" % (phase(Sr[0, 1]), phase(Sa[0, 1]))
    if  passed:
        msg = "[passed] " + msg
    else:
        msg = "[failed] " + msg
    print msg 

#test_waveguide_S()
eval_filter6_30_a()
