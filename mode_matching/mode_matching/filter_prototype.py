#!/usr/bin/env python2

import numpy as np
from scipy.signal import cheby1, butter
from matplotlib import pyplot as plt
import sys


z, p, k = cheby1(5, 0.1, 1, btype = 'lowpass', analog = True, output = 'zpk')

p += 5e-3

pB = np.poly1d(k*np.poly(z))
pA = np.poly1d(np.poly(p))

def minus_poly(p):
    order = p.order
    return np.poly1d(p.coeffs * np.fromfunction(lambda n: (-1)**(order-n), (order+1,)))

pAp = minus_poly(pA)
pBp = minus_poly(pB)

P = pB*pBp
Q = pA*pAp

error = P.deriv()*Q - P*Q.deriv()

extremae = np.abs(P(error.roots)/Q(error.roots))
gain = 1.0/np.max(extremae)

print 'gain = ', gain
print 'k = ', k

omega = np.linspace(0, 3, num = 1001)
H = gain*np.abs(pB(1j*omega)/pA(1j*omega))**2
print "Hmin = ", 10*np.log10(H[0])

plt.plot(omega, 10*np.log10(H))
plt.grid()
plt.show()

# renormalized B
pB  = np.poly1d(gain * pB.coeffs)
pBp = np.poly1d(gain * pBp.coeffs)


pP2 = pA*pAp - pB*pBp
r = pP2.roots

r_upper = r[np.where(r.real < 0.0)]
if 0+0j in r:
    r_upper = np.append(r_upper, 0.0)
    
pP = np.poly1d(np.poly(r_upper))
pPp = minus_poly(pP)
pP2p = pP*pPp

k = np.sqrt(pP2.coeffs[0]/pP2p.coeffs[0])

pP = np.poly1d(k*pP.coeffs)
pQ = pA

ZB = pQ + pP
ZA = pQ - pP

order = ZB.order
for i in xrange(ZB.order):
    q, r = ZB / ZA
    print 'g[{}] = {}'.format(i+1, q.coeffs[0])
    ZB = ZA
    ZA = r

print 'g[{}] = {}'.format(order+1, 1/q.coeffs[1])

