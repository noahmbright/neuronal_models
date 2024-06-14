# Pitt MATH 3370 Mathematical Neuroscience Models

An collection of models for use in Pitt Math 3370 Mathematical Neuroscience. :brain:

# Installation
How *do* you install this package?

## Dependencies

Numpy, Scipy, maybe a C compiler if I get around to the stochastic coupling

# Supported models
## Hodgkin-Huxley
4D model in $V, m, h, n$. Default initial conditions are $(V, m, h, n) = (-65,
0.05, 0.6, 0.317)$. 

The dynamics of the model are given by

$$
c\frac{dV}{dt} = I_0 + I_{\\text{inj}} - I_{\text{Na}} - I_{\text{K}} -
g_{\text{Leak}}(V - V_{\\text{Leak}}),
$$

$$
\frac{dm}{dt} = \alpha_{m}(V)(1 - m) - \beta_m(V)m,
$$

$$
\frac{dh}{dt} = \alpha_{h}(V)(1 - h) - \beta_h(V)h,
$$

$$
\frac{dn}{dt} = \alpha_{n}(V)(1 - n) - \beta_n(V)n,
$$

with the $\alpha$ and $\beta$ given by

$$
\alpha_m(V) = \phi \times .1 \times \frac{(V + 40)}{1 - \exp\left[-(V +
40)/10\right]}
\qquad
\beta_m(V) = \phi\times 4 \exp\left[-(V + 65)/18\right]
$$

$$
\alpha_h(V) = \phi \times 0.07\exp\left[-(V + 65)/20\right]
\qquad 
\beta_h(V) = \phi\frac 1 {1 + \exp\left[-(V + 35)/10\right]}
$$

$$
\alpha_n(V) = \phi\times 0.01 \frac{V + 55}{\exp\left[-(V + 55)/10\right]}
\qquad
\beta_n(V) = \phi\times0.125\exp\left[-(V + 65)/80\right].
$$

Finally, the currents are given by
$$
I_{\text{Na}} = g_{\text{Na}}h(V - V_{\text{Na}}) m^3, 
\quad 
I_{\text{K}} = g_{\\text{K}} (V - V_{\text{K}}) n^4.
$$



These are default values for parameterizing the $\alpha$ and $\beta$. 

TBD: should these be expressed in a more generic way? 

The reversal potentials and conductances default to

$$
V_{\text{Na}} = 50 \quad V_{\text{K}} = -77 \quad V_{\text{Leak}} = -54.387
$$

$$
g_{\text{Na}} = 120 \quad g_{\text{K}} = 36 \quad g_{\text{Leak}} = 0.3.
$$

$c$ and $\phi$ default to 1. 

## Rinzel Reduction

The Rinzel reduction of the HH model inherits the $\alpha$ and $\beta$
functions, as well as the reversal potentials, channel conductances, $c$ and
$\phi$. The dynamics for $V$ and $n$ are also inherited, but $h$ and $m$ are
modeled: 

$$
\frac{dm}{dt} = \frac{\alpha_{m}(V)}{\alpha_m(V) + \beta_m(V)},
$$

$$
h = h_0 - n.
$$

$h_0$ is a constant that defaults to 0.8.

## Kepler Reduction

Where did I get these equations from? Canvas?

## Integrate-and-fire
### Quadratic

## Destexhe-Pare

Destexhe-Pare is a 5D model in $V, m, h, n$ and $m_{\text{K}}$, with default
initial conditions 

$$
(V, m, h, n, m_{\text{K}}) = (-73.87,0,1,0.002,0.0075).
$$

The dynamics of $m, n,$ and $h$ are the same as in $HH$, with the dynamics of
$V$ and $m_{\text{K}}$ given by:

$$
c\frac{dV}{dt} = I_0 + I_\text{Inj}} - g_{\text{L}}(V - E_\text{L}}) -
I_{\text{Kdr}}(V,n) - I_{\text{Na}}(V, m, h) - I_{\text{Km}}
$$

$$
\frac{dm_{\text{K}}}{dt} = \alpha_{m_{\text{K}}}(V)(1 - m_{\text{K}}) - \beta_{m_{\text{K}}}(V)m_{\text{K}}.
$$

## Izhikevic

The Izhikevic model is 2D in $u$ and $V$, with dynamics

$$
\dot V = I_0 + I_{\text{inj}} + V^2 - u,
$$

$$
\dot u = a(bV - u).
$$

On $V$ crossing $V_{\text{threshold}}$ from below, the resetting

$$
V\to c, u \to u + d
$$

takes place. $V_{\text{threshold}}$, $a, b, c$, and $d$ are all adjustable
member variables of the model. 

TODO: Need reasonable defaults.


## Morris-Lecar

# How to couple models
How *do* you couple models?

# Other functionality

## Gillespie Simulations
Currently supports 
