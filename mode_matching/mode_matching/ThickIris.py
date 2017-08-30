from GenericElement import GenericElement
from WaveguideElement import WaveguideElement
from WaveguideJunction import WaveguideJunction
from Utils import toSI as SI

class ThickIris(GenericElement):
    """
    ThickIris represents S matrix of a thick iris in a waveguide
    """

    def __init__(self, widthA, widthB, length, numModesA, numModesB, frequency):
        self.widthA, self.widthB = SI(widthA, widthB)
        self.length = SI(length)
        self.frequency = SI(frequency)
        self.numModesA = numModesA
        self.numModesB = numModesB

        self.wj = WaveguideJunction(self.widthA,
                                    self.widthB,
                                    self.numModesA,
                                    self.numModesB,
                                    self.frequency)
        self.wg = WaveguideElement(self.widthB,
                                   self.length,
                                   self.numModesB,
                                   self.frequency)
        self.update(frequency, first_time = True)

    def update(self, frequency, first_time = False):
        if (not first_time):
            self.frequency = SI(frequency)
            self.wj.update(frequency)
            self.wg.update(frequency)

            
        wjp = GenericElement(self.wj.s22,
                             self.wj.s21,
                             self.wj.s12,
                             self.wj.s11)
        bb = self.wj * self.wg * wjp
        self.s11, self.s12, self.s21, self.s22 = bb.s11, bb.s12, bb.s21, bb.s22

    def __repr__(self):
        fmt = \
"""
Rectangular waveguide junction, 
    Left: width = {widthA}m, num. modes = {modesA},
    Right: width = {widthB}m, num. modes = {modesB},
    Thickness = {thickness}m,
    Frequency = {frequency}Hz
""".strip()
        return fmt.format(widthA = self.widthA, widthB = self.widthB,
                          modesA = self.numModesA, modesB = self.numModesB,
                          thickness = self.length,
                          frequency = self.frequency)
