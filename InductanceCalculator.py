''' This class calculates the off-diagonal components of the inductance matrix for a system of coils. '''

import numpy as np
from scipy.integrate import dblquad

class InductanceCalculator:
    def __init__(self, *, coil_system, initial_L):
        self.coil_system = coil_system
        self.L = np.copy(initial_L) 
        # this is assumed to have valid self-inductances but possibly wrong mutual inductances (e.g. from moving coils about)
        self.N = self.coil_system.NUM_COILS

    def compute_L(self, write_loc):
        for i in range(self.N):
            for j in range(i + 1, self.N):
                R1, Z1 = self.coil_system.coil_coords[i]
                R2, Z2 = self.coil_system.coil_coords[j]

                computed_Lij = self.__calc_Lij(R1, R2, Z1 - Z2)
                self.L[i, j] = computed_Lij
                self.L[j, i] = computed_Lij

        np.savetxt(write_loc, self.L, delimiter=',')
        return self.L

    def __calc_Lij(self, R1, R2, Z):
        # takes in meters, returns Henries

        MU_0 = 4 * np.pi * 1e-7 

        def integrand(theta1, theta2):
            numer = np.cos(theta1) * np.cos(theta2) + np.sin(theta1) * np.sin(theta2)
            denom = np.sqrt((R1 * np.cos(theta1) - R2 * np.cos(theta2)) ** 2 + \
                    (R1 * np.sin(theta1) - R2 * np.sin(theta2)) ** 2 + Z ** 2)
            return numer / denom
        
        integral = dblquad(
            integrand, 0, 2 * np.pi, 
            lambda theta1: 0, lambda theta1: 2 * np.pi
        )[0]

        return MU_0 / (4 * np.pi) * R1 * R2 * integral
