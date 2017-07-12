from GenericElement import GenericElement
import Utils
import numpy as np
import numpy.linalg as lin

class WaveguideJunction(GenericElement):
    """
    WaveguideJunction represents S matrix of a symmetric rectangular waveguide junction
    """
    def __init__(self, widthA, widthB, numModesA, numModesB, frequency):
        self.widthA, self.widthB, self.frequency = Utils.toSI(widthA, widthB, frequency)
        self.numModesA, self.numModesB = numModesA, numModesB

        if abs(self.widthA) > abs(self.widthB):
            (self.s11, self.s12, self.s21, self.s22) = \
            self.__compute_S__(self.widthA, self.widthB,
                               self.numModesA, self.numModesB,
                               self.frequency)
        else:
            (self.s22, self.s21, self.s12, self.s11) = \
            self.__compute_S__(self.widthB, self.widthA,
                               self.numModesB, self.numModesA,
                               self.frequency)
            

    def __compute_S__(self, a, b, Ma, Mb, frequency):
        """
        Compute S matrix elements assuming that a > b
        """
        Q = np.fromfunction(lambda m, n: Utils.OverlappingTEn0(m, n, a, b), (Ma, Mb))
        P = np.fromfunction(lambda m, n:
                            Utils.GammaTEn0(n, a, frequency) / Utils.GammaTEn0(m, b, frequency) * Utils.OverlappingTEn0(n, m, a, b),
                            (Mb, Ma))

        QP = np.dot(Q, P)
        PQ = np.dot(P, Q)
        X = lin.inv(np.eye(Ma) + QP)
        Y = lin.inv(np.eye(Mb) + PQ)
        s11 = np.dot(X, -np.eye(Ma) + QP)
        s12 = 2 * np.dot(X, Q)
        s21 = 2 * np.dot(Y, P)
        s22 = np.dot(Y, np.eye(Mb) - PQ)

        return (s11, s12, s21, s22)

    def __repr__(self):
        fmt = \
"""
Rectangular waveguide junction, 
    Left: width = {widthA}m, num. modes = {modesA},
    Right: width = {widthB}m, num. modes = {modesB},
    Frequency = {frequency}Hz
""".strip()
        return fmt.format(widthA = self.widthA, widthB = self.widthB,
                          modesA = self.numModesA, modesB = self.numModesB,
                          frequency = self.frequency)
