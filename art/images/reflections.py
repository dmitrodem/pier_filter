#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
from cycler import cycler
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib2tikz import save as tikz_save
from matplotlib.patches import Ellipse, Arrow
from matplotlib.collections import Collection
import pylab


monochrome = (cycler('color', ['k']) * cycler('marker', ['']) *
              cycler('linestyle', ['-', '--', ':'
              ]))

use_tex = True
if use_tex:
    from matplotlib import rc
    rc('font',**{'family':'serif'})
    rc('text', usetex=True)
    rc('text.latex',unicode=True)
    rc('text.latex',preamble='\usepackage[utf8]{inputenc}')
    rc('text.latex',preamble='\usepackage[russian]{babel}')
    rc('axes', prop_cycle=monochrome)

def run():
    plt.figure()
    d1 = np.loadtxt('../data/01_mode_matching.txt')
    d2 = np.loadtxt('../data/02_equivalent_circuit.txt')
    d3 = np.loadtxt('../data/03_hfss.txt')

    ax = plt.gca()
    
    ax.plot(d1[:, 0], d1[:, 1], label = u'Mode matching') 
    ax.plot(d2[:, 0], d2[:, 1], label = u'Equivalent circuit model')
    ax.plot(d3[:, 0], d3[:, 1], label = u'HFSS model')

    ann = []
    mark1 = (0.93, -75)
    ell = Ellipse(xy = mark1,
                  width = 0.005,
                  height = 30,
                  facecolor = 'none',
                  edgecolor = 'black',
                  linestyle = 'solid',
                  linewidth = 1)
    ax.add_artist(ell)
    ax.annotate('', xy = mark1, xytext = (0.905, mark1[1]), arrowprops=dict(arrowstyle="<-"))
    
    ax.set_xlim([0.9, 1.1])
    ax.set_ylim([-100, 5])
    plt.xticks(np.arange(0.9, 1.1, 0.05))
    plt.grid()
    plt.xlabel('$F/F_{center}')
    plt.ylabel('$S_{21}$, dB')
    ax2 = plt.gca().twinx()
    ax2.set_ylabel('$S_{11}$, dB')

    conv_s21 = lambda x: 10*np.log10(1-10**(0.1*x))
    ax2.plot(d1[:, 0], conv_s21(d1[:, 1]), label = u'Mode matching') 
    ax2.plot(d2[:, 0], conv_s21(d2[:, 1]), label = u'Equivalent circuit model')
    ax2.plot(d3[:, 0], conv_s21(d3[:, 1]), label = u'HFSS model')

    mark2 = (1.08, -3)
    ell = Ellipse(xy = mark2,
                  width = 0.005,
                  height = 10,
                  facecolor = 'none',
                  edgecolor = 'black',
                  linestyle = 'solid',
                  linewidth = 1)
    ax.add_artist(ell)
    ax.annotate('', xy = mark2, xytext = (1.095, mark2[1]), arrowprops=dict(arrowstyle="<-"))
    
    ax2.set_ylim([-100, 5])
    ax2.legend(loc = 'center right')
    
    axins = inset_axes(ax, '70%', '30%', loc = 'lower right')
    axins.set_xlim(0.97, 1.03)
    axins.set_ylim(-0.2, 0.1)
    axins.plot(d1[:, 0], d1[:, 1])
    axins.plot(d2[:, 0], d2[:, 1])
    axins.plot(d3[:, 0], d3[:, 1])
    axins.get_xaxis().set_visible(False)
    plt.gcf().set_size_inches(np.array([17, 15]) / 2.54)
    plt.tight_layout()
    plt.grid()
    plt.savefig('reflections.pdf')
    tikz_save('reflections_tikz.tex')
    plt.show()


run()
