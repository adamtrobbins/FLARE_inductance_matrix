import numpy as np
from CoilSystem import CoilSystem
from InductanceSolver import InductanceSolver
from magnet_formulae import ring_Br, ring_Bz

class TransientResponse:
    def __init__(self, *, coil_system, inductance_solution):
        self.coil_system = coil_system
        self.sol = inductance_solution
        self.N = len(self.sol.t)

    def tindex(self, t):
        return np.abs(self.sol.t - t).argmin()

    def calculate_field(self, tindex, r, z):
        currents = self.sol.y[:, tindex]

        Br = 0
        Bz = 0

        for i in range(self.coil_system.NUM_COILS):
            [coil_r, coil_z] = self.coil_system.coil_coords[i]
            z_delta = z - coil_z

            Br += self.__calc_Br(currents[i], coil_r, r, z_delta)
            Bz += self.__calc_Bz(currents[i], coil_r, r, z_delta)

        return Br, Bz
    
    def calc_field_waveform(self, r, z):
        Brs = np.zeros(self.N)
        Bzs = np.zeros(self.N)

        for i in range(self.N):
            Brs[i], Bzs[i] = self.calculate_field(i, r, z)

        return Brs, Bzs

    def __calc_Br(self, I, R, r, z_delta):
        return ring_Br(I, R, z_delta, r)

    def __calc_Bz(self, I, R, r, z_delta):
        return ring_Bz(I, R, z_delta, r)   

