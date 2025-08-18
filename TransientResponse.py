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

    def calculate_field(self, tindex, r, z, activated_coils = None):
        currents = self.sol.y[:, tindex]

        Br = 0
        Bz = 0

        if activated_coils is None:
            activated_coils = range(self.coil_system.NUM_COILS)
        else:
            activated_coils = self.coil_system.get_cindices(activated_coils)

        for i in activated_coils:
            [coil_r, coil_z] = self.coil_system.coil_coords[i]
            z_delta = z - coil_z

            inc_Br, inc_Bz = self.__calc_B(currents[i], coil_r, r, z_delta, spread = self.coil_system.coil_spreads[i])

            Br += inc_Br
            Bz += inc_Bz

        return Br, Bz
    
    def calc_field_waveform(self, r, z):
        Brs = np.zeros(self.N)
        Bzs = np.zeros(self.N)

        for i in range(self.N):
            Brs[i], Bzs[i] = self.calculate_field(i, r, z)

        return Brs, Bzs

    def __calc_B(self, I, R, r, z_delta, spread = [np.nan, np.nan], N = 20):
        if (np.isnan(I)):
            return 0, 0

        # Spread [z1, z2] gives range over which to spread out the coil from the center position, N the # of discretization points
        zmin, zmax = spread
        if np.isnan(zmin) or np.isnan(zmax):
            return ring_Br(I, R, z_delta, r), ring_Bz(I, R, z_delta, r)
        else:
            Br = 0
            Bz = 0
            
            z_vals = np.linspace(zmin, zmax, N)

            for z in z_vals:
                Br += ring_Br(I / N, R, z_delta + z, r)
                Bz += ring_Bz(I / N, R, z_delta + z, r)

            return Br, Bz

