import numpy as np
from scipy.integrate import solve_ivp


def square_current_pulse(I_0, t_on=0, t_off=-1.0, width=-1.0):
    if t_off == -1.0 and width == -1.0:
        print("Either pulse width or time off must be set")
    if t_off != 1.0 and width != -1.0 and t_off - t_on != width:
        print(
            "Both pulse width and t_off set, but difference between t_on and t_off != width"
        )
    if t_off != -1.0 and t_off < t_on:
        print("Pulse off time before on time")

    if width != -1.0:
        t_off = t_on + width

    return lambda t: I_0 * (t > t_on and t < t_off)


def integrate_with_reset(dVdt, t0, tf, V_reset, V_threshold):
    threshold_crossed = lambda V: V[0] - V_threshold
    V_start = [V_reset]
    t_start = t0
    ts = np.array([])
    Vs = np.array([])
    spike_ts = np.array([])

    while len(ts) == 0 or ts[-1] < tf:
        sol = solve_ivp(
            dVdt,
            (t_start, tf),
            V_start,
            events=(threshold_crossed),
        )

        ts = np.concatenate((ts, sol.t))
        Vs = np.concatenate((Vs, sol.y[0]))
        spike_ts = np.concatenate((spike_ts, sol.t_events[0]))

        t_start = ts[-1]
        V_start = [V_reset]

    return ts, Vs, spike_ts
