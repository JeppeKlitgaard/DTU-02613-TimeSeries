#import "@preview/physica:0.9.4": super-T-as-transpose
#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#import "utils.typ": mathformatter, mref

// Font fix because FruityFeedback is trash
#set text(font: "STIX Two Text")
#show math.equation: set text(font: "STIX Two Math")

#show: codly-init.with()
#codly(languages: codly-languages)


#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show: super-T-as-transpose
#set table.header()

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: true, upright: true)

#let Cov = math.op("Cov")

///////////////////////////////////////////////////////////////////////////////

#align(center, text(17pt)[
  *02417 Times Series Analysis - Assignment 2*
])

Authors:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-03-24

#pagebreak()

= Stability <sec:1_stability>

We are given an AR(2) process ${X_t}$ with residual term $ε_t$ arising from a stochastic process ${ε_t} ∈ 𝒩(0, σ_ε^2=1)$ :
$
  y_t + ϕ_1 y_(t-1) + ϕ_2 y_(t-2) = ε_t
$ <eq:1_ar2>

== Stationarity

To prove stationarity we must derive the _characteristic equation_ of the AR(2) process.
This is done by taking the Z-transform of the AR(2) process.

We first restrict ourselves to the _homogeneous_ part of the equation by removing $ε(t)$
and rewrite past operations using the backshift operator, $B y_t ≡ y_(t-1)$:

$
  y_t + ϕ_1 B y_t + ϕ_2 B^2 y_t = 0
$

We consider $B = z^(-1)$ where $B$ is understood to be defined in the time-domain of the process while $z$ belongs to the complex frequency domain. Additionally, we factor out $y_t$ and multiply through by $z^2$ to arrive at the _characteristic equation_:

$
  y_t (z^2 + ϕ_1 z + ϕ_2) = 0
$ <eq:1_characteristic_eq>

We then solve for the roots of @eq:1_characteristic_eq:

$
  z^2 + ϕ_1 z+ ϕ_2 = 0
$

This is a quadratic equation in $z$ and can be solved using the quadratic formula:

$
  z_± = (-ϕ_1 ± sqrt(ϕ_1^2 - 4ϕ_2)) / 2
$

Plugging in $ϕ_1 = -0.7, ϕ_2 = -0.2$ from the assignment we find:
$
  z_± ∈ {-0.218, 0.918}
$ <eq:1_roots>

Recalling the condition on stationarity to be that all roots of the process are within the unit circle, or equivalently:
$
|z| <= 1 wide ∀ med z
$

Which is satisfied by inspection of @eq:1_roots.

As such, we conclude that the process is stationary for $ϕ_1 = -0.7, ϕ_2 = -0.2$.

== Invertibility

Invertibility refers to the ability to represent current errors as a function of past observations,
which is inherently obtained for auto-regressive processes.
To realise this, we simply inspect @eq:1_ar2 and find that the error term has already been expressed as a function of past observations, thus satisfying the definition of invertibility.

== Autocorrelation

We recall the definition of _autocorrelation_ for a _stationary process_:
$
  ρ(k) = γ(k)/γ(0)
$ <eq:1_autocorrelation>

Where $γ(k)$ is the _autocovariance_ for a timeshift $k$:
$
  γ(k) = Cov[y_t, y_(t+k)]
$ <eq:1_autocovariance>

We consider the autocorrelation of an AR(2) process by solving for $y_t$ in @eq:1_ar2:
$
  y_t = -ϕ_1 y_(t-1) - ϕ_2 y_(t-2) + ε_t
$ <eq:1_ar2_solved>

And then inserting this in @eq:1_autocorrelation:
$
  ρ(k) = Cov[y_t, y_(t+k)] / γ(0)
    &= Cov[y_t, -ϕ_1 y_(t-1+k) - ϕ_2 y_(t-2+k) + ε_t] / γ(0)\
    &= Cov[y_t, -ϕ_1 y_(t-1+k)] / γ(0) + Cov[y_t, -ϕ_2 y_(t-2+k)] / γ(0) + cancel(Cov[y_t, ε_t] / γ(0))\
    &= -ϕ_1 Cov[y_t, y_(t-1+k)] / γ(0) - ϕ_2 Cov[y_t, y_(t-2+k)] / γ(0)\
    &= -ϕ_1 γ(k-1) / γ(0) - ϕ_2 γ(k-2) / γ(0)\
  ρ(k) &= -ϕ_1 ρ(k-1) - ϕ_2 ρ(k-2)\
$ <eq:1_autocorrelation_recursion>

Notably by stationary it follows $ρ(-k) = ρ(k)$ and from @eq:1_autocorrelation we find $ρ(0) = 1$,
which allows us to build a recursive relation for $ρ(k)$ from:

$
  ρ(0) &= 1\
  ρ(1) &= -ϕ_1 ρ(0) - ϕ_2 ρ(-1) = -ϕ_1 - ϕ_2ρ(1) = (-ϕ_1)/(1+ϕ_2)\
$

== Autocorrelation Plot

Computing the sequence of autocorrelations with different lags $k∈ {0, 1, ..., 30}$ is then trivial using @eq:1_autocorrelation_recursion, the outcome of which is shown in @fig:1_4_autocorrelation_plot.

#figure(
  image("output/1_4_autocorrelation_function.png"),
  caption: [Autocorrelation function of AR(2) process with $ϕ_1 = -0.7, ϕ_2 = -0.2$.
    Notice exponential decay.],
) <fig:1_4_autocorrelation_plot>

While the @fig:1_4_autocorrelation_plot shows the correct autocorrelation function
for the assignment description as can be trivially verified by inspection of @eq:1_ar2_solved,
it is not beyond the realm of possibility that the assignment description is incorrect and the
intended parameters were $ϕ_1 = 0.7, ϕ_2 = 0.2$, which we have also plotted in @fig:1_4_autocorrelation_plot_incorrect.

#figure(
  image("output/1_4_autocorrelation_function_incorrect.png"),
  caption: [Autocorrelation function of AR(2) process with $ϕ_1 = 0.7, ϕ_2 = 0.2$.
    Cycling behavior and exponential decay]
) <fig:1_4_autocorrelation_plot_incorrect>

These two figures yield an intuition for the role of the parameters,
where those that lead to a negated relation to the previous observations
give rise to the cyclical behaviour found in @fig:1_4_autocorrelation_plot_incorrect.
This will be investigated further in the following sections.

= Simulating seasonal processes <sec:2_simulating_seasonal_processes>

In order to simulate the seasonal processes we utilise the Python library `statsmodels`,
which offers an implementation of the `SARIMAX` (Seasonal AutoRegressive Integrated Moving Average with eXogenous regressors) model.

= Identifying ARMA(p,q) Models <sec:3_arma_identification>

