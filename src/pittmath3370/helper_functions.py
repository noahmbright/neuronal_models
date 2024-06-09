import numpy as np
from scipy.integrate import solve_ivp
from math import floor


def square_current_pulse(I_0, t_on=0, t_off=-1.0, width=-1.0):
    if t_off == -1.0 and width == -1.0:
        raise Exception("Either pulse width or time off must be set")
    if t_off != 1.0 and width != -1.0 and t_off - t_on != width:
        raise Exception(
            "Both pulse width and t_off set, but difference between t_on and t_off != width"
        )
    if t_off != -1.0 and t_off < t_on:
        raise Exception("Pulse off time before on time")

    if width != -1.0:
        t_off = t_on + width

    return lambda t: I_0 * (t > t_on and t < t_off)


def integrate_with_reset(dVdt, t0, tf, V_reset, V_threshold):
    if t0 > tf:
        raise Exception("t0 is after tf")
    threshold_crossed = lambda V: V[0] - V_threshold
    threshold_crossed.terminal = True
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


def _square_wave(t, A, duty_cycle, T, t_end):
    if t > t_end:
        return 0

    quotient = floor(t / T)
    t_mod = t - quotient * T
    return A * (t_mod < (duty_cycle * T))


def square_wave(A, duty_cycle, T, t_end):
    return lambda t: _square_wave(t, A, duty_cycle, T, t_end)


def _infinity_func(V, Vth, sigma):
    return 1 / (1 + np.exp((V - Vth) / sigma))


def infinity_func(Vth, sigma):
    return lambda V: _infinity_func(V, Vth, sigma)


def gillespie(alpha, beta, t_max, n0, N):
    t = 0
    if np.random.uniform() < n0:
        n = N
    else:
        n = np.random.randint(0, N - 1)

    ns = [int(n == N)]
    ts = [0]
    while t < t_max:
        r = (N - n) * alpha + beta * n
        rho = np.random.uniform()
        t += -1 / r * np.log(rho)
        x = np.random.uniform(0, r)

        if x < alpha * (N - n):
            n = min(N, n + 1)
        else:
            n = max(0, n - 1)

        ts.append(t)
        ns.append(int(n == N))

    return np.array(ts), np.array(ns)
