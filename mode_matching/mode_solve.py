#!/usr/bin/env python3

from pylab import *
from sys import stdout
from time import time

Mmax = 10
Nmax = 10
a = 1.0+0j
b = 0.5001+0j
kv = 2.5*pi*a

def Imn(m, n, a, b):
    return 2*(-1)**(1+n)*a**2*b*cos(0.5*pi*b*(2*m-1)/a)*(2*n-1)/pi/(-b**2*(2*m-1)**2+a**2*(2*n-1)**2)

def beta(n, l, kv):
    return sqrt(kv**2-pi**2*(2*n-1)**2/l**2)

t0 = time()
s11 = zeros(1001)
i = 0
kv_space = linspace(0, 5*pi/a, 1001)
for kv in kv_space:
    Q = fromfunction(lambda m, n: 2.0/a*Imn(m+1, n+1, a, b), (Mmax, Nmax))
    P = fromfunction(lambda m, n: 2.0/b*beta(n+1, a, kv)/beta(m+1, b, kv)*Imn(n+1, m+1, a, b), (Nmax, Mmax))

    M = dot(inv(eye(Mmax) + dot(Q, P)), (-eye(Mmax) + dot(Q, P)))
    s11[i] = 20*log10(abs(M[0,0]))
    i += 1

plot(kv_space/pi*a, s11)
grid()
show()

ai = zeros(Mmax, dtype=complex); ai[0] = 1.0
ar = dot(M, ai)
br = dot(P, ai - ar)

print(ar)
print(br)

xs = linspace(-0.5*a, 0.5*a, 1001)
lhs = zeros(xs.shape, dtype=complex)
rhs = zeros(xs.shape, dtype=complex)

for m in range(Mmax):
    lhs += ar[m]*cos(pi/a*(2*m+1)*xs)
    
for n in range(Nmax):
    rhs += br[n]*cos(pi/b*(2*n+1)*xs)

plot(xs, abs(cos(pi/a*xs)+lhs))
plot(xs, abs(rhs))
grid()
show()
