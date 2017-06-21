#!/usr/bin/env python2

from pylab import *
from sys import stdout
from time import time

Mmax = 10
Nmax = 10
a = 1.0001+0j
b = 0.5001+0j
kv = 2.5*a

def Imn(m, n, a, b):
    return 4*(-1)**(n)*a**1.5*b**0.5*cos(0.5*pi*b*(2*m+1)/a)*(2*n+1)/pi/(-b**2*(2*m+1)**2+a**2*(2*n+1)**2)

def beta(n, l, kv):
    return sqrt(1.0-((2*n+1)/l/kv)**2)

t0 = time()
s11 = zeros(1001)
i = 0
kv_space = linspace(0, 5/a, 1001)
kv_space = array([0j, 3 + 0j])
for kv in kv_space[1:]:
    Q = fromfunction(lambda m, n: Imn(m, n, a, b), (Mmax, Nmax))
    P = fromfunction(lambda m, n: beta(n, a, kv)/beta(m, b, kv)*Imn(n, m, a, b), (Nmax, Mmax))

    M = dot(inv(eye(Mmax, dtype=complex) + dot(Q, P)), (-eye(Mmax, dtype=complex) + dot(Q, P)))
    s11[i] = 20*log10(abs(M[0,0]))
    i += 1

# plot(kv_space*a, s11)
# grid()
# show()

ai = zeros(Mmax, dtype=complex); ai[0] = 1.0
ar = dot(M, ai)
br = dot(P, ai - ar)


xs = linspace(-0.5*a, 0.5*a, 1001)
lhsE = zeros(xs.shape, dtype=complex)
rhsE = zeros(xs.shape, dtype=complex)
lhsH = zeros(xs.shape, dtype=complex)
rhsH = zeros(xs.shape, dtype=complex)
lhsS = 0j
rhsS = 0j

for m in range(Mmax):
    lhsE +=  ai[m]*cos(pi/a*(2*m+1)*xs)*sqrt(2.0/a)
    lhsE +=  ar[m]*cos(pi/a*(2*m+1)*xs)*sqrt(2.0/a)
    lhsH +=  beta(m, a, kv) * ai[m]*cos(pi/a*(2*m+1)*xs)*sqrt(2.0/a)
    lhsH += -beta(m, a, kv) * ar[m]*cos(pi/a*(2*m+1)*xs)*sqrt(2.0/a)
    lhsS += beta(m, a, kv)*(ai[m]**2 - ar[m]**2)
    
for n in range(Nmax):
    rhsE += br[n]*cos(pi/b*(2*n+1)*xs)*sqrt(2.0/b)
    rhsH += beta(n, b, kv) * br[n]*cos(pi/b*(2*n+1)*xs)*sqrt(2.0/b)
    rhsS += beta(n, b, kv) * br[n]**2

subplot(211)
plot(xs, abs(lhsE))
plot(xs, abs(rhsE))
grid()
subplot(212)
plot(xs, abs(lhsH))
plot(xs, abs(rhsH))
grid()
show()

delta = (lhsS - rhsS) / (lhsS + rhsS) * 2
print "Energy imbalance = ", np.abs(delta)
# np.set_printoptions(precision=3, suppress=True, linewidth=250)
# print("Showing table P : ")
# print(abs(P))

with open("/tmp/Q.npy", "w") as fd:
    np.save(fd, Q)

print "S11 = ", 20*log10(abs(M[0, 0]))
print "ang(S11) = ", np.rad2deg(np.angle(M[0,0]))
print "S11 = ", 20*log10(abs(M[0, 1]))
print "ang(S11) = ", np.rad2deg(np.angle(M[0,1]))
