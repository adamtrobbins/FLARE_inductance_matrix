import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

class InductanceSolver:
    def __init__(self, *, coil_system, waveform_dir, feedback_supression = 1e4):
        self.coil_system = coil_system
        self.NUM_COILS = coil_system.NUM_COILS

        # Master PF waveform

        pf_waveform = np.loadtxt(f'{waveform_dir}/pf_waveform.csv', delimiter=',', skiprows=1)

        self.pf_times = pf_waveform[:, 0]
        self.pf_currents = pf_waveform[:, 1] * 1e3 # assumes csv data in kiloamps

        self.master_pf = interp1d(self.pf_times, self.pf_currents, kind = 'linear', bounds_error = True)

        self.current_dist = np.loadtxt(f'{waveform_dir}/current_dist.csv')
        self.current_dist = self.coil_system.prune_coil_array(self.current_dist)

        M = coil_system.L / coil_system.R[:, np.newaxis]
        driven_coils = np.where(self.current_dist != 0)[0]
        M[driven_coils, :] /= feedback_supression

        self.M_inv = np.linalg.inv(M)

    def solve_system(self, *, ntimesteps):
        I0 = lambda t: self.current_dist * self.master_pf(t)

        t_span = (0, self.pf_times[-1])
        I_init = np.zeros(self.NUM_COILS)

        def dI_dt(t, I):
            return self.M_inv @ (I0(t) - I)

        sol = solve_ivp(dI_dt, t_span, I_init, t_eval=np.linspace(*t_span, ntimesteps))

        return sol 

