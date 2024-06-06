import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def ddt(t, y, I, a, b):
    v, u = y
    dvdt = I(t) + v**2.0 - u
    dudt = a * (b * v - u)
    return [dvdt, dudt]


def v_cross(t, y, I, a, b):
    return y[0] - v_threshold


v_cross.terminal = True


def I_step(t, I, t_on):
    return I * (t > t_on)


v_threshold = 1
a = 0.05
b = 0.4
c = -1
d = 5

I_amp = 5
t_on = 10
current = lambda t: I_step(t, I_amp, t_on)

t0 = 0
tf = 100
u0 = 0
v0 = c


def integrate_izhikevich():
    t = t0
    y0 = [v0, u0]
    ys = np.empty((2, 0))
    ts = np.array([])
    spike_ts = np.array([])
    i = 0
    while len(ts) == 0 or ts[-1] < tf:
        if i % 100 == 0 and i > 0:
            print(f"iteration {i}, latest time: {ts[-1]}")
        i += 1
        sol = solve_ivp(ddt, (t, tf), y0, events=v_cross, args=(current, a, b))

        ts = np.concatenate((ts, sol.t))
        ys = np.concatenate((ys, sol.y), axis=1)
        spike_ts = np.concatenate((spike_ts, sol.t_events[0]))

        y0 = [c, sol.y[1][-1] + d]
        t = ts[-1]

    return ts, ys, spike_ts


ts, ys, spike_ts = integrate_izhikevich()
vs = ys[0]
plt.plot(ts, vs)
plt.show()
