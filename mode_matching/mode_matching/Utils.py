from pint import UnitRegistry
import numpy as np
from numpy.lib.scimath import sqrt # handles sqrt(-1) as 1j
from scipy.constants import c as c0

"""
Various utilities
"""

__ureg__ = UnitRegistry(system = 'mks')

def __convSI(x):
    return __ureg__.parse_expression(x).to_base_units().magnitude

def toSI(*args):
    """
    Convert input arguments to SI units (returns tuple of conversion results)
    """
    result = []
    for arg in args:
        try:
            result.append(__convSI(arg))
        except:
            result.append(arg)
    if len(result) == 1:
        return result[0]
    else:
        return tuple(result)

def EigenfunctionTEn0(n, l):
    """
    Eigenfunction for TE(2*n+1, 0) mode in rectangular waveguide of width ``l``
    """
    return lambda x: np.sqrt(2.0/l) * np.cos(np.pi*(2*n+1)*x/l)

def GammaTEn0(n, l, freq):
    """
    Propagation coefficients for TE(2*n+1, 0) mode in rectangular waveguide 
    of width ``l``
    """
    return np.pi * sqrt(((2*n+1)/l)**2 - (2 * freq / c0)**2)

def OverlappingTEn0(m, n, a, b):
    """
    Calculate matrix elements Imn = \int_{S} Eigenfunction(m, a) * Eigenfunction(n, b) dx
    NB: ensure that a > b before calling this function
    """
    return 4*(-1)**(n)*a**1.5*b**0.5*np.cos(0.5*np.pi*b*(2*m+1)/a)*(2*n+1)/np.pi/(-b**2*(2*m+1)**2+a**2*(2*n+1)**2)
