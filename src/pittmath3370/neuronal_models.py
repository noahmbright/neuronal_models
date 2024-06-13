import numpy as np
from scipy.integrate import solve_ivp
from abc import ABC, abstractmethod
# from typing import override

# import helper_functions


# stuff I can do easily
# TODO tweak infinity functions to take a hyperpolarization parameter
# and you would add that by subtracting off the hyperpolarization in both alpha and beta right?
# TODO allow clamping for infinity functions in case you're asked to make alpha inf big or something

# stuff I don't know how to do
# TODO make phi sophisticated
# TODO stability analysis? If yes, figure out how to do that
# TODO fix Morris Lecar parameters - what is right?

# stuff I could do if asked
# TODO add the Fitzhugh Nagumo ones from homework 4?
# TODO add regular Fitzhugh Nagumo?


class NeuronalModel(ABC):
    # TODO make phi do what phi is actually supposed to do
    phi = 1

    def __init__(self):
        self.dimensions = 0
        self.t0 = 0
        self.I0 = 0
        self.tf = 100
        self.I_inj = lambda: 0
        self.cm = 1.0  # uF/cm^2

    def set_I0(self, I0):
        self.I0 = I0

    def set_I_inj(self, I_inj):
        self.I_int = I_inj

    def set_t0(self, t0):
        self.t0 = t0

    def set_tf(self, tf):
        self.tf = tf

    def set_cm(self, cm):
        self.cm = cm

    @staticmethod
    @abstractmethod
    def dALLdt(t, X, self, clamp):
        return np.array([0.0])

    def integrate(self, ics, clamp=False):
        return solve_ivp(self.dALLdt, (self.t0, self.tf), ics, args=(self, clamp))


class AbstractHHLike(NeuronalModel):
    def __init__(self):
        super().__init__()
        # mS/cm^2
        self.g_Na = 120.0
        self.g_K = 36.0
        self.g_L = 0.3
        # mV
        self.E_Na = 50.0
        self.E_K = -77.0
        self.E_L = -54.387

    def alpha_m(self, V):
        return 0.1 * (V + 40.0) / (1.0 - np.exp(-(V + 40.0) / 10.0))

    def beta_m(self, V):
        return 4.0 * np.exp(-(V + 65.0) / 18.0)

    def alpha_h(self, V, hyperpolarization=0):
        V_half = -65.0 - hyperpolarization
        return 0.07 * np.exp(-(V - V_half) / 20.0)

    def beta_h(self, V, hyperpolarization=0):
        V_half = -35.0 - hyperpolarization
        return 1.0 / (1.0 + np.exp(-(V - V_half) / 10.0))

    def alpha_n(self, V):
        return 0.01 * (V + 55.0) / (1.0 - np.exp(-(V + 55.0) / 10.0))

    def beta_n(self, V):
        return 0.125 * np.exp(-(V + 65) / 80.0)

    def tau_h(self, V, hyperpolarization):
        return self.alpha_h(V, hyperpolarization) + self.beta_h(V, hyperpolarization)

    def tau_m(self, V):
        return self.alpha_m(V) + self.beta_m(V)

    def tau_n(self, V):
        return self.alpha_n(V) + self.beta_n(V)

    def m_inf(self, V):
        return self.alpha_m(V) / self.tau_m(V)

    def n_inf(self, V):
        return self.alpha_n(V) / self.tau_n(V)

    def h_inf(self, V, hyperpolarization=0):
        return self.alpha_h(V, hyperpolarization) / self.tau_h(V, hyperpolarization)

    def I_Na(self, V, m, h):
        return self.g_Na * m**3 * h * (V - self.E_Na)

    def I_K(self, V, n):
        return self.g_K * n**4 * (V - self.E_K)

    def I_L(self, V):
        return self.g_L * (V - self.E_L)


class HodgkinHuxley(AbstractHHLike):
    def __init__(self):
        super().__init__()

    @staticmethod
    def dALLdt(t, X, self, clamp=False):
        V, m, h, n = X

        if clamp:
            dVdt = 0
        else:
            dVdt = (
                self.I0
                + self.I_inj(t)
                - self.I_Na(V, m, h)
                - self.I_K(V, n)
                - self.I_L(V)
            ) / self.cm

        dmdt = self.alpha_m(V) * (1.0 - m) - self.beta_m(V) * m
        dhdt = self.alpha_h(V) * (1.0 - h) - self.beta_h(V) * h
        dndt = self.alpha_n(V) * (1.0 - n) - self.beta_n(V) * n
        return np.array([dVdt, dmdt, dhdt, dndt])


class Rinzel(AbstractHHLike):
    """
    inherit all the alpha/beta from HH
    use the approximation that h = h0 - n
    """

    def __init__(self):
        super().__init__()
        self.h0 = 0.8

    def set_h0(self, h0):
        self.h0 = h0

    def h(self, n):
        return self.h0 - n

    @staticmethod
    def dALLdt(t, X, self, clamp=False):
        V, n = X
        if clamp:
            dVdt = 0
        else:
            dVdt = (
                self.I_inj(t)
                + self.I0
                - self.I_Na(V, self.m_inf(V), self.h(n))
                - self.I_K(V, n)
                - self.I_L(V)
            ) / self.cm
        dndt = self.alpha_n(V) * (1.0 - n) - self.beta_n(V) * n

        return np.array([dVdt, dndt])


class Kepler(AbstractHHLike):
    def __init__(self):
        super().__init__()

    def dh_inf(self, V, dV=0.01):
        return (self.h_inf(V + dV) - self.h_inf(V - dV)) / (2 * dV)

    @staticmethod
    def dALLdt(t, X, self, clamp=False):
        V, Vh = X
        if clamp:
            dVdt = 0
        else:
            dVdt = (
                self.I_inj(t)
                + self.I0
                - self.I_Na(V, self.m_inf(V), self.h_inf(Vh))
                - self.I_K(V, self.n_inf(Vh))
                - self.I_L(V)
            ) / self.cm

        dVhdt = (self.h_inf(V) - self.h_inf(Vh)) / (
            self.dh_inf(Vh) / (self.alpha_h(V) + self.beta_h(V))
        )

        return np.array([dVdt, dVhdt])


class DestexhePare(NeuronalModel):
    def __init__(self):
        super().__init__()
        self.I0 = 0
        self.E_Na = 55
        self.E_k = -85
        self.g_K = 100
        self.gkm = 2

    def set_gkm(self, gkm):
        self.gkm = gkm

    V_t = -58
    V_s = -10

    def alpha_m(self, V):
        return -0.32 * (V - self.V_t - 13) / (np.exp(-(V - self.V_t - 13) / 4) - 1)

    def beta_m(self, V):
        return 0.28 * (V - self.V_t - 40) / (np.exp((V - self.V_t - 40) / 5) - 1)

    def alpha_h(self, V):
        return 0.128 * np.exp(-(V - self.V_t - self.V_s - 17) / 18)

    def beta_h(self, V):
        return 4 / (1 + np.exp(-(V - self.V_t - self.V_s - 40) / 5))

    def alpha_n(self, V):
        return -0.032 * (V - self.V_t - 15) / (np.exp(-(V - self.V_t - 15) / 5) - 1)

    def beta_n(self, V):
        return 0.5 * np.exp(-(V - self.V_t - 10) / 40)

    def alpha_km(self, V):
        return 0.0001 * (V + 30) / (1 - np.exp(-(V + 30) / 9))

    def beta_km(self, V):
        return -0.0001 * (V + 30) / (1 - np.exp((V + 30) / 9))

    def I_km(self, V, m):
        return self.gkm * m * (V - self.E_k)

    @staticmethod
    def dALLdt(t, X, self, clamp=False):
        V, m, h, n, mk = X

        if clamp:
            dVdt = 0
        else:
            dVdt = (
                self.I_inj(t)
                + self.I0
                - self.I_Na(V, m, h)
                - self.I_K(V, n)
                - self.I_L(V)
                - self.I_km(V, mk)
            ) / self.cm
        dmdt = self.alpha_m(V) * (1.0 - m) - self.beta_m(V) * m
        dhdt = self.alpha_h(V) * (1.0 - h) - self.beta_h(V) * h
        dndt = self.alpha_n(V) * (1.0 - n) - self.beta_n(V) * n
        dmkdt = self.alpha_km(V) * (1.0 - mk) - self.beta_km(V) * mk
        return np.array([dVdt, dmdt, dhdt, dndt, dmkdt])


# TODO check what to inherit from
class MorrisLecar(NeuronalModel):
    gk = 2.0
    gl = 0.5
    Vk = -0.7
    Vl = -0.5
    om = 1  # ?
    # TODO I_Ca has a V-1
    E_Ca = 1  # correct?
    # TODO Cm?

    def __init__(self):
        super().__init__()
        self.V1 = -0.01
        self.V2 = 0.15
        self.V3 = 0.1
        self.V4 = 0.145
        self.I0 = 0
        self.phi = 0.333
        self.gca = 1

    def set_phi(self, phi):
        self.phi = phi

    def m_inf(self, V):
        return 0.5 * (1 + np.tanh((V - self.V1) / self.V2))

    def n_inf(self, V):
        return 0.5 * (1 + np.tanh((V - self.V3) / self.V4))

    def lam_n(self, V):
        return self.phi * np.cosh((V - self.V3) / (2 * self.V4))

    def I_Ca(self, V):
        return self.gca * self.m_inf(V) * (V - self.E_Ca)

    def I_L(self, V):
        return self.gl * (V - self.Vl)

    def I_K(self, V, w):
        return self.gk * w * (V - self.Vk)

    def dwdt(self, V, w):
        return self.lam_n(V) * (self.n_inf(V) - w)

    @staticmethod
    def dALLdt(t, X, self, clamp=False):
        V, w = X

        if clamp:
            dVdt = 0
        else:
            dVdt = self.I0 - self.I_K(V, w) - self.I_L(V) - self.I_Ca(V) + self.I_inj(t)
        dwdt = self.dwdt(V, w)

        return np.array([dVdt, dwdt])
