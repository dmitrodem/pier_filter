from pyswarm import pso
import numpy as np
import Utils as u
from Filter import Filter
from GenericElement import GenericElement
from matplotlib import pyplot as plt

a	= u.toSI("8.636mm")
b	= u.toSI("3.556mm")
y0	= u.toSI("20mm")
c1	= u.toSI("3.375mm")
c2	= u.toSI("3.563mm")
c3	= u.toSI("3.784mm")
c4	= u.toSI("4.924mm")
d1	= u.toSI("4.925mm")
d2	= u.toSI("4.722mm")
d3	= u.toSI("3.77mm")
dd	= u.toSI("2.5mm")
a1	= u.toSI("7.112mm") # const
ddd	= u.toSI("6mm") #const

fmin = u.toSI("26GHz")
fmax = u.toSI("36GHz")
freqs = np.linspace(fmin, fmax, num = 101)
min_modes = 5

f_low  = u.toSI("29.3GHz")
f_high = u.toSI("30.9GHz")
f_band  = np.where(np.logical_and(freqs >= f_low, freqs <= f_high))
f_outband = np.where(np.logical_or(freqs < f_low, freqs > f_high))

def cost(freqs, s12):
    s12[f_band] -= 1.0
    return np.trapz(np.abs(s12), freqs)

def filter_eval(x, return_s21 = False):

    c4, c3, c2, c1, d3, d2, d1, dd = x
    descr = [
        (a1, ddd),
        (c4, dd),
        (a, d3),
        (c3, dd),
        (a, d2),
        (c2, dd),
        (a, d1),
        (c1, 0.5*dd)
        ]

    s21 = np.zeros(freqs.shape, dtype = np.complex)

    f = Filter(descr, min_modes, freqs[0])
    for i, freq in enumerate(freqs):
        f.update(freq)
        fp = GenericElement(f.s22, f.s21, f.s12, f.s11)
        f1 = f*fp
        s21[i] = (f1.s21[0,0])

    if return_s21:
        return (20*np.log10(np.abs(s21)), np.unwrap(np.angle(s21)))
    else:
        return cost(freqs, np.abs(s21))

lb = np.array([c4, c3, c2, c1, d3, d2, d1, dd]) * 0.5
ub = np.array([c4, c3, c2, c1, d3, d2, d1, dd]) * 2

if True:
    xopt, fopt = pso(filter_eval, lb, ub, swarmsize = 1000, maxiter = 1000, debug = True)

    print xopt[0] * 1e3
    print fopt

freqs = np.linspace(fmin, fmax, num = 1001)

#xopt = [ 0.00486402,  0.00355182,  0.00433982,  0.00444804]

s21, phase = filter_eval(xopt, return_s21 = True)
plt.subplot(211)
plt.plot(freqs/1e9, s21)
plt.grid(True)
plt.subplot(212)
plt.plot(freqs/1e9, phase/2/np.pi)
plt.grid(True)
plt.show()
