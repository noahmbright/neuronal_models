import neuronal_models as nm
from scipy.integrate import solve_ivp
import numpy as np


class Izhikevic(nm.NeuronalModel):
    def __init__(self):
        super().__init__()
        self.v_threshold = 1
        self.a = 1
        self.b = 1
        self.c = 1
        self.d = 1

    @staticmethod
    def dALLdt(t, X, self, clamp):
        v, u = X
        if clamp:
            dvdt = 0
        else:
            dvdt = self.I0 + self.I_inj(t) + v**2.0 - u
        dudt = self.a * (self.b * self.v - u)
        return np.array([dvdt, dudt])

    def v_cross(t, X, self, clamp):
        return X[0] - self.v_threshold

    def integrate(self, ics, clamp=False):
        y0 = ics
        ys = np.empty((2, 0))
        ts = np.array([])
        t0 = self.t0
        spike_ts = np.array([])
        while len(ts) == 0 or ts[-1] < self.tf:
            sol = solve_ivp(
                self.dALLdt, (t0, self.tf), y0, events=self.v_cross, args=(self, clamp)
            )

            ts = np.concatenate((ts, sol.t))
            ys = np.concatenate((ys, sol.y), axis=1)
            spike_ts = np.concatenate((spike_ts, sol.t_events[0]))

            y0 = [self.c, sol.y[1][-1] + self.d]
            t0 = ts[-1]

        return ts, ys, spike_ts
