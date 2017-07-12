from pyswarm import pso
import numpy as np
import Utils as u
from Filter import Filter
from GenericElement import GenericElement
from matplotlib import pyplot as plt

a	= u.toSI("8.636mm")
c1	= u.toSI("3.375mm")
c2      = u.toSI("3.563mm")
dd	= u.toSI("2.5mm")
d1	= u.toSI("5.925mm")
d2      = u.toSI("4.722mm")
ddd	= u.toSI("6mm")

fmin = u.toSI("26GHz")
fmax = u.toSI("36GHz")
freqs = np.linspace(fmin, fmax, num = 101)
min_modes = 5

f_low  = u.toSI("29.3GHz")
f_high = u.toSI("30.9GHz")
f_band  = np.where(np.logical_and(freqs >= f_low, freqs <= f_high))

def cost(freqs, s12):
    s12[f_band] -= 1.0
    return np.trapz(np.abs(s12), freqs)

def filter_eval(x, return_s21 = False):
    d1 = x[0]
    c1 = x[1]
    d2 = x[2]
    c2 = x[3]
    descr = [
        (a, ddd),
        (c2, dd),
        (a, d2),
        (c1, dd),
        (a, 0.5*d1)
        ]

    s21 = np.zeros(freqs.shape, dtype = np.complex)

    for i, freq in enumerate(freqs):
        f = Filter(descr, min_modes, freq)
        f.calculate()
        fp = GenericElement(f.s22, f.s21, f.s12, f.s11)
        f1 = f*fp
        s21[i] = (f1.s21[0,0])

    if return_s21:
        return (20*np.log10(np.abs(s21)), np.unwrap(np.angle(s21)))
    else:
        return cost(freqs, np.abs(s21))

lb = [u.toSI("2mm"), u.toSI("1mm"), u.toSI("2mm"), u.toSI("1mm")]
ub = [u.toSI("7mm"), u.toSI("7mm"), u.toSI("7mm"), u.toSI("7mm")]

if False:
    xopt, fopt = pso(filter_eval, lb, ub, swarmsize = 100, maxiter = 100, debug = True)

    print xopt[0] * 1e3
    print fopt

freqs = np.linspace(fmin, fmax, num = 1001)

xopt = [ 0.00486402,  0.00355182,  0.00433982,  0.00444804]

s21, phase = filter_eval(xopt, return_s21 = True)
plt.subplot(211)
plt.plot(freqs/1e9, s21)
plt.grid(True)
plt.subplot(212)
plt.plot(freqs/1e9, phase/2/np.pi)
plt.grid(True)
plt.show()
