"""
Main demo for the Hodgkin Huxley neuron model
"""
from neuronal_models import HodgkinHuxley
import matplotlib.pyplot as plt

hh = HodgkinHuxley()

X = hh.integrate()
ts = X.t
V = X.y[0]
m = X.y[1]
h = X.y[2]
n = X.y[3]

ina = hh.I_Na(V, m, h)
ik = hh.I_K(V, n)
il = hh.I_L(V)

plt.figure()

plt.subplot(4,1,1)
plt.title('Hodgkin-Huxley Neuron')
plt.plot(ts, V, 'k')
plt.ylabel('V (mV)')

plt.subplot(4,1,2)
plt.plot(ts, ina, 'c', label='$I_{Na}$')
plt.plot(ts, ik, 'y', label='$I_{K}$')
plt.plot(ts, il, 'm', label='$I_{L}$')
plt.ylabel('Current')
plt.legend()

plt.subplot(4,1,3)
plt.plot(ts, m, 'r', label='m')
plt.plot(ts, h, 'g', label='h')
plt.plot(ts, n, 'b', label='n')
plt.ylabel('Gating Value')
plt.legend()

plt.subplot(4,1,4)
i_inj_values = [hh.I_inj(t) for t in ts]
plt.plot(ts, i_inj_values, 'k')
plt.xlabel('t (ms)')
plt.ylabel('$I_{inj}$ ($\\mu{A}/cm^2$)')
plt.ylim(-1, 40)

plt.show()
