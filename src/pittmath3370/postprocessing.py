import neuronal_models
import numpy as np
from scipy.signal import find_peaks


def _get_period_if_steady_state(Vs, ts, tolerance=1e-9, min_periods=10):
    peaks, _ = find_peaks(Vs)
    if len(peaks) < min_periods + 1:
        return None
    width = ts[peaks[-1]] - ts[peaks[-2]]

    for idx in peaks[-1:-min_periods]:
        if np.abs(ts[peaks[idx]] - ts[peaks[idx - 1]]) > tolerance:
            return None
    return 1 / width


def compute_FI_curve(model, Is):
    Fs = np.zeros(len(Is))
    for i in range(len(Is)):
        model.set_I0(Is[i])
        X = model.integrate()

        T = _get_period_if_steady_state()
        Fs[i] = 1 / T

    return Fs


def compute_peak_IV(model, Vs):
    n = len(Vs)
    # TODO implement a steady state function
    ic = lambda V: [V, hh.m_inf(V_eq), hh.h_inf(V_eq), hh.n_inf(V_eq)]
    peak_Is = np.zeros(n)
    for i in range(n):
        X = model.integrate(I_0=0, ics=ic(Vs[i]), clamp=True)
        I_Na = hh.I_Na(X[:, 0], X[:, 1], X[:, 2])
        peak_Is[i] = np.max(np.abs(I_Na))
