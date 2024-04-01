import numpy as np
from scipy.integrate import solve_ivp
#from typing import override

# TODO make phi sophisticated
# TODO add Izhikevich
# TODO add the Fitzhugh Nagumo ones from homework 4?
# TODO add regular Fitzhugh Nagumo?
# TODO add generic integrate with reset
# TODO stability analysis
# TODO frequency curves? 

class HodgkinHuxley():

    #mS/cm^2
    g_Na = 120.0
    g_K  =  36.0
    g_L  =   0.3
    
    # mV
    E_Na =  50.0
    E_K  = -77.0
    E_L  = -54.387

    # TODO make phi do what phi is actually supposed to do
    phi = 1
    

    def __init__(self):
        self.t_on = 0
        self.I_p = 0
        self.pulse_width = 0    
        self.cm = 1.0 # uF/cm^2
        self.t0 = 0
        self.tf = 400
        self.I_0 = 0

    def set_I0(self, I_0):
        self.I_0 = I_0

    def set_t0(self, t0):
        self.t0 = t0

    def set_tf(self, tf):
        self.tf = tf

    def set_cm(self, cm):
        self.cm = cm
    
    def alpha_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.1*(V+40.0)/(1.0 - np.exp(-(V+40.0) / 10.0))

    def beta_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 4.0*np.exp(-(V+65.0) / 18.0)

    def alpha_h(self, V, hyperpolarization = 0):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_half = -65.0 - hyperpolarization
        return 0.07*np.exp(-(V-V_half) / 20.0)

    def beta_h(self, V, hyperpolarization = 0):
        """Channel gating kinetics. Functions of membrane voltage"""
        V_half = -35.0 - hyperpolarization
        return 1.0/(1.0 + np.exp(-(V-V_half) / 10.0))

    def alpha_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.01*(V+55.0)/(1.0 - np.exp(-(V+55.0) / 10.0))

    def beta_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.125*np.exp(-(V+65) / 80.0)

    def tau_h(self, V, hyperpolarization):
        return self.alpha_h(V, hyperpolarization) + self.beta_h(V, hyperpolarization)

    def tau_m(self, V):
        return self.alpha_m(V) + self.beta_m(V)

    def tau_n(self, V):
        return self.alpha_n(V) + self.beta_n(V)

    def m_inf(self, V):
        return self.alpha_m(V)/self.tau_m(V)

    def n_inf(self, V):
        return self.alpha_n(V)/self.tau_n(V)

    def h_inf(self, V, hyperpolarization=0):
        return self.alpha_h(V, hyperpolarization)/self.tau_h(V, hyperpolarization)

    def I_Na(self, V, m, h):
        return self.g_Na * m**3 * h * (V - self.E_Na)

    def I_K(self, V, n):
        return self.g_K  * n**4 * (V - self.E_K)
    
    def I_L(self, V):
        return self.g_L * (V - self.E_L)

    #injected current
    def set_current_pulse(self, I_p, t_on, pulse_width):
        self.t_on = t_on
        self.I_p = I_p 
        self.pulse_width = pulse_width

    def I_inj(self, t):
        """
        External Current

        |  :param t: time
        |  :return: step up to 10 uA/cm^2 at t>100
        |           step down to 0 uA/cm^2 at t>200
        |   this is a pulse of current and you should modify it
        """
        return self.I_p*(t>self.t_on) - self.I_p*(t>(self.t_on + self.pulse_width)) 

    @staticmethod
    def dALLdt(t, X, self, clamp = False):
        """

        |  :param X:
        |  :param t:
        |  :return: calculate membrane potential & activation variables
        """
        V, m, h, n = X

        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + self.I_0 - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V)) / self.cm
        dmdt = self.alpha_m(V)*(1.0-m) - self.beta_m(V)*m
        dhdt = self.alpha_h(V)*(1.0-h) - self.beta_h(V)*h
        dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n
        return np.array([dVdt, dmdt, dhdt, dndt])

    def integrate(self, ics = [-60, 0.09, 0.4, 0.4], clamp = False):
        return solve_ivp(self.dALLdt, (self.t0, self.tf), ics, args = (self, clamp))
        

class Rinzel(HodgkinHuxley):
    """
    inherit all the alpha/beta from HH
    use the approximation that h = h0 - n
    """
    def __init__(self):
        super().__init__()
        self.h0 = .8

    def set_h0(self, h0):
        self.h0 = h0

    def h(self, n):
        return self.h0 - n

    @staticmethod
    def dALLdt(t, X, self, I_0 = 0, clamp = False):
        V, n = X
        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + self.I_0 - self.I_Na(V, self.m_inf(V), self.h(n)) - self.I_K(V, n) - self.I_L(V)) / self.cm
            dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n

        return np.array([dVdt, dndt])


class Kepler(HodgkinHuxley):
    def __init__(self):
        super().__init__()

    def dh_inf(self, V, dV = .01):
        return (self.h_inf(V + dV) - self.h_inf(V - dV))/(2*dV)

    @staticmethod
    def dALLdt(t, X, self, I_0 = 0, clamp = False):
        V, Vh = X
        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + self.I_0 - self.I_Na(V, self.m_inf(V), self.h_inf(Vh)) - self.I_K(V, self.n_inf(Vh)) - self.I_L(V)) / self.cm
            dVhdt = (self.h_inf(V) - self.h_inf(Vh))/(self.dh_inf(Vh)/(self.alpha_h(V) + self.beta_h(V)))

        return np.array([dVdt, dVhdt])



class DestexhePare(HodgkinHuxley):
    def __init__(self):
        super().__init__()
        self.E_Na = 55
        self.E_k = -85
        self.g_K = 100
        self.gkm = 2

    def set_gkm(self, gkm):
        self.gkm = gkm

    V_t = -58
    V_s = -10

    def alpha_m(self, V):
        return -.32*(V-self.V_t-13)/(np.exp(-(V-self.V_t-13)/4)-1)

    def beta_m(self, V):
        return .28*(V-self.V_t-40)/(np.exp((V-self.V_t-40)/5)-1)

    def alpha_h(self, V):
        return .128*np.exp(-(V-self.V_t-self.V_s-17)/18)

    def beta_h(self, V):
        return 4/(1+np.exp(-(V-self.V_t-self.V_s-40)/5))

    def alpha_n(self, V):
        return -.032*(V-self.V_t-15)/(np.exp(-(V-self.V_t-15)/5)-1)

    def beta_n(self, V):
        return .5*np.exp(-(V-self.V_t-10)/40)

    def alpha_km(self, V):
        return .0001*(V+30)/(1-np.exp(-(V+30)/9))

    def beta_km(self, V):
        return -.0001*(V+30)/(1-np.exp((V+30)/9))

    def I_km(self, V, m):
        return self.gkm * m * (V-self.E_k)

    @staticmethod
    def dALLdt(t, X, self, clamp = False):
        """
        Integrate

        |  :param X:
        |  :param t:
        |  :return: calculate membrane potential & activation variables
        """
        V, m, h, n, mk = X

        if clamp:
            dVdt = 0
        else:
            dVdt = (self.I_inj(t) + self.I_0 - self.I_Na(V, m, h) - self.I_K(V, n) - self.I_L(V) - self.I_km(V, mk)) / self.cm
        dmdt = self.alpha_m(V)*(1.0-m) - self.beta_m(V)*m
        dhdt = self.alpha_h(V)*(1.0-h) - self.beta_h(V)*h
        dndt = self.alpha_n(V)*(1.0-n) - self.beta_n(V)*n
        dmkdt = self.alpha_km(V)*(1.0-mk) - self.beta_km(V)*mk
        return np.array([dVdt, dmdt, dhdt, dndt, dmkdt])


class MorrisLecar(HodgkinHuxley):
    gk = 2.0
    gl = 0.5
    Vk = -0.7
    Vl = -0.5
    om = 1 #?
    #TODO I_Ca has a V-1
    E_Ca = 1 # correct?
    #TODO Cm?

    def __init__(self):
        super().__init__()
        self.V1 = -0.01
        self.V2 = 0.15
        self.V3 = 0.1
        self.V4 = 0.145
        self.phi = .333
        self.gca = 1

    def set_phi(self, phi):
        self.phi = phi

    def activate_SNIC(self):
        self.V3 = .1
        self.V4 = 1.45
        self.gca = 1
        self.phi = .333
        self.I_0 = 0

    def activate_hopf(self):
        self.V3 = 0
        self.V4 = .3
        self.gca = 1.1
        self.phi = .2
        self.I_0 = 0

    def activate_homoclinic_below_bifurcation(self):
        self.V3 = 0.1
        self.V4 = 0.145
        self.gca = 1.0
        self.phi = 1.15
        self.I_0 = 0

    def activate_homoclinic_above_bifurcation(self):
        self.V3 = 0.1
        self.V4 = 0.145
        self.gca = 1.0
        self.phi = 1.15
        self.I_0 = 0.08

    def m_inf(self, V):
        return 0.5 * (1 + np.tanh((V - self.V1)/self.V2))

    def n_inf(self, V):
        return 0.5 * ((1 + np.tanh((V - self.V3)/self.V4)))

    def lam_n(self, V):
        return self.phi * np.cosh((V - self.V3)/(2*self.V4))

    def I_Ca(self, V):
        return self.gca * self.m_inf(V) * (V-self.E_Ca)

    def I_L(self, V):
        return self.gl * (V - self.Vl)

    def I_K(self, V, w):
        return self.gk * w * (V - self.Vk)

    def dVdt(self, V, w, t):
        return self.I_0 - self.I_K(V, w) - self.I_L(V) - self.I_Ca(V) + self.I_inj(t)

    def dwdt(self, V, w):
        return self.lam_n(V)* ( self.n_inf(V) - w )

    @staticmethod
    def dALLdt(t, X, self, clamp = False):
        V, w = X

        if clamp:
            dVdt = 0
        else:
            dVdt = self.dVdt(V, w, t)
        dwdt = self.dwdt(V, w)

        return np.array([dVdt, dwdt])


class CoupledML(MorrisLecar):
    def __init__(self, n):
        super().__init__()
        self.n = n
        self.het = 1.0
        self.Is = np.zeros(n)
        self.gc = 0

    def set_gc(self, gc):
        self.gc = gc

    def set_Is(self, Is):
        assert (len(Is) == self.n), "Length Is != number of coupled neurons n"
        self.Is = Is

    def dVdt(self, I0, V, w):
        return I0 - self.I_K(V, w) - self.I_L(V) - self.I_Ca(V)

    @staticmethod
    def dALLdt(t, X, self, clamp = False):
        Vs = X[:self.n]
        ws = X[self.n:]

        if clamp:
            dVsdt = np.zeros(n)
        else:
            dVsdt = np.array([self.dVdt(I, V, w) for I, V, w in zip(self.Is, Vs, ws)])

        # TODO generalize to more than two cells with a matrix or something
        # figure out how several cells couple
        print(f"coupling {self.gc*(Vs[1] - Vs[0])}")
        dVsdt[0] += self.gc*(Vs[1] - Vs[0])
        dVsdt[1] += self.gc*(Vs[0] - Vs[1])/self.het

        dwsdt = np.array([self.dwdt(V, w) for V, w in zip(Vs, ws)])

        return np.concatenate((dVsdt, dwsdt))


