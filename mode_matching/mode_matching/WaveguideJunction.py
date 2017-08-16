from GenericElement import GenericElement
import Utils
import numpy as np
import numpy.linalg as lin

class WaveguideJunction(GenericElement):
    """
    WaveguideJunction represents S matrix of a symmetric rectangular waveguide junction
    """
    def __init__(self, widthA, widthB, numModesA, numModesB, frequency):
        self.widthA, self.widthB = Utils.toSI(widthA, widthB)
        self.numModesA, self.numModesB = numModesA, numModesB

        self.direct =  abs(self.widthA) > abs(self.widthB) 
        
        if self.direct:
            self.Q = np.fromfunction(lambda m, n: Utils.OverlappingTEn0(m, n, self.widthA, self.widthB), (self.numModesA, self.numModesB))
        else:
            self.Q = np.fromfunction(lambda m, n: Utils.OverlappingTEn0(m, n, self.widthB, self.widthA), (self.numModesB, self.numModesA))

        self.eye_modesA = np.eye(self.numModesA)
        self.eye_modesB = np.eye(self.numModesB)
        self.update(frequency)

    def update(self, frequency):
        self.frequency = Utils.toSI(frequency)
        if self.direct:
            self.P = np.fromfunction(lambda m, n: Utils.GammaTEn0(n, self.widthA, frequency) / Utils.GammaTEn0(m, self.widthB, frequency), (self.numModesB, self.numModesA)) * self.Q.T
            QP = np.dot(self.Q, self.P)
            PQ = np.dot(self.P, self.Q)
            X = lin.inv(self.eye_modesA + QP)
            Y = lin.inv(self.eye_modesB + PQ)
            self.s11 = np.dot(X, -self.eye_modesA + QP)
            self.s12 = 2 * np.dot(X, self.Q)
            self.s21 = 2 * np.dot(Y, self.P)
            self.s22 = np.dot(Y, self.eye_modesB - PQ)
        else:
            self.P = np.fromfunction(lambda m, n: Utils.GammaTEn0(n, self.widthB, frequency) / Utils.GammaTEn0(m, self.widthA, frequency), (self.numModesA, self.numModesB)) * self.Q.T
            QP = np.dot(self.Q, self.P)
            PQ = np.dot(self.P, self.Q)
            X = lin.inv(self.eye_modesB + QP)
            Y = lin.inv(self.eye_modesA + PQ)
            self.s22 = np.dot(X, -self.eye_modesB + QP)
            self.s21 = 2 * np.dot(X, self.Q)
            self.s12 = 2 * np.dot(Y, self.P)
            self.s11 = np.dot(Y, self.eye_modesA - PQ)

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
