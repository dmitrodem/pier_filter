from GenericElement import GenericElement
import Utils
import numpy as np

class WaveguideElement(GenericElement):
    """
    WaveguideElement represents S matrix of a simple waveguide section
    """
    def __init__(self, width, length, numModes, frequency):
        """
        Initialize S matrix from know geometric properties (width and length of a waveguide),
        frequency and number of modes (numModes)
        """
        self.width, self.length, self.frequency = Utils.toSI(width, length, frequency)
        gamma_vect = np.fromfunction(lambda m: Utils.GammaTEn0(m, self.width, self.frequency), (numModes, ))
        reflection  = np.zeros((numModes, numModes), dtype = complex)
        propagation = np.diag(np.exp(-gamma_vect * self.length))
        self.s11 = reflection
        self.s12 = propagation
        self.s21 = propagation
        self.s22 = reflection
    
    def __repr__(self):
        fmt = \
"""
Rectangular waveguide 
    Width = {width}m, length = {length}m
    Frequency = {frequency}Hz
""".strip()
        return fmt.format(
            width = self.width,
            length = self.length,
            frequency = self.frequency)
