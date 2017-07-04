import numpy as np
from numpy.linalg import inv

class GenericElement(object):
    """
    GenericElement represents basic class for mode matching
    technique. It has generic constructor which builds generalized
    S matrix and defines a way to "multiply" S matrices (i.e. to 
    cascade networks)
    """
    def __init__(self, s11, s12, s21, s22):
        """
        Initialize S matrix from its bloc elements
        """
        self.s11 = s11
        self.s12 = s12
        self.s21 = s21
        self.s22 = s22
        
    def __mul__(self, A):
        """
        Embed current element S matrix with element A S matrix
        """
        assert(self.s22.shape == A.s11.shape)
        N = self.s22.shape[0]
        X = inv(np.eye(N) - np.dot(A.s11, self.s22))
        Y = inv(np.eye(N) - np.dot(self.s22, A.s11))

        s11 = self.s11 + np.dot(self.s12, np.dot(X, np.dot(A.s11, self.s21)))
        s12 = np.dot(self.s12, np.dot(X, A.s12))
        s21 = np.dot(A.s21, np.dot(Y, self.s21))
        s22 = A.s22 + np.dot(A.s21, np.dot(Y, np.dot(self.s22, A.s12)))

        return GenericElement(s11, s12, s21, s22)
