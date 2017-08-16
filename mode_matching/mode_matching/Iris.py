from GenericElement import GenericElement
from WaveguideJunction import WaveguideJunction
from WaveguideElement import WaveguideElement
from Utils import toSI as SI
import numpy as np
from matplotlib import pyplot as plt
from scipy.constants import c as c0

a = SI("8.636mm")
l = SI("100.0mm")
dd = SI("1mm")
f_c = 0.5*c0/a

fmin = SI("17.5GHz")
fmax = SI("30GHz")
frange = np.linspace(fmin, fmax, num = 1001)
s12 = np.zeros(frange.shape, dtype = complex)

iris_widths = np.array([0.11, 0.21, 0.31, 0.41, 0.51, 0.61, 0.66, 0.71]) * a
impedances = np.zeros(iris_widths.shape)

irisModes = 10

plt.figure()

dd_lengths = np.linspace(0, SI("2.5mm"), num = 11)

for dd_idx, dd in enumerate(dd_lengths):
    for idx, iris in enumerate(iris_widths):
        waveguideModes = np.round(irisModes * a/iris)
        waveguide = WaveguideElement(a, l, waveguideModes, frange[0])
        junction = WaveguideJunction(a, iris, waveguideModes, irisModes, frange[0])
        waveguide1 = WaveguideElement(iris, 0.5*dd, irisModes, frange[0])
        for i, freq in enumerate(frange):
            waveguide.update(freq)
            junction.update(freq)
            waveguide1.update(freq)
            r1 = waveguide * junction * waveguide1
            r2 = GenericElement(r1.s22, r1.s21, r1.s12, r1.s11)
            r3 = r1*r2
            s12[i] = r3.s12[0,0]


        x = 1.0/(1.0 - (f_c/frange)**2)
        y = np.abs(s12)**(-2) - 1.0

        p = np.polyfit(x, y, 1)
        impedances[idx] = 0.5/np.sqrt(p[0])


    plt.plot(iris_widths/a, impedances, label = ("dd = %f" % dd))

plt.xlabel("iris width [waveguide widths]")
plt.ylabel("Impedance [Z0]")
plt.legend()
plt.grid()
plt.title("Iris load inductance vs. iris width")
plt.show()
    
