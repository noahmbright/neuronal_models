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
\frac{dV}{dt} = I_0 + I_{\text{inj}} - g_{\text{Na}}h(V - V_{\text{Na}}) m^3 -
g_{\text{K}} (V - V_{\text{K}}) n^4 - \frac 1 {c} g_{\text{Leak}}(V -
V_{\text{Leak}}),
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
\beta_n(V) = \phi\times0.125\exp\left[-(v + 65)/80\right]
$$

These are default values for parameterizing the $\alpha$ and $\beta$. 

TBD: should these be expressed in a more generic way? 

The reversal potentials and conductances default to

$$
V_{text{Na}} = 50 \quad V_{text{K}} = -77 \quad V_{text{Leak}} = -54.387
$$

$$
g_{text{Na}} = 120 \quad g_{text{K}} = 36 \quad g_{text{Leak}} = 0.3.
$$

$c$ and $\phi$ default to 1. 

## Kepler Reduction

## Rinzel Reduction

## Integrate-and-fire
### Quadratic
## Izhikevic
## Morris-Lecar

# How to couple models
How *do* you couple models?

# Other functionality

## Gillespie Simulations
Currently supports 
