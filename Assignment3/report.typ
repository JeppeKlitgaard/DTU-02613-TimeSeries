#import "@preview/physica:0.9.4": super-T-as-transpose
#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#import "utils.typ": mathformatter, mref

#set page(
  paper: "a4",
  numbering: "1/1",
  header: context {
    if(counter(page).get().at(0)== 1) [

    ] else [
    #set align(center)
    #set text(10pt, fill: gray)
      02417 Times Series Analysis - Assignment 3
  ]
  }
)

// Font fix because FruityFeedback is trash
#set text(font: "STIX Two Text")
#show math.equation: set text(font: "STIX Two Math")

#show: codly-init.with()
#codly(languages: codly-languages)


#set heading(numbering: "1.1.a")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show: super-T-as-transpose
#set table.header()

#let vv = mathformatter(underbars: 0, bold: true, upright: false)
#let mm = mathformatter(underbars: 0, bold: true, upright: true)

#let Cov = math.op("Cov")

///////////////////////////////////////////////////////////////////////////////

#align(center, text(17pt)[
  *02417 Times Series Analysis - Assignment 3*
])

Authors:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-04-21

#pagebreak()

= Simulations <sec:1_stability>

We are given an AR(2) process ${X_t}$ with residual term $ε_t$ arising from a stochastic process ${ε_t} ∈ 𝒩(0, σ_ε^2=1)$ :
$
  X_t + ϕ_1 X_(t-1) + ϕ_2 X_(t-2) = ε_t
$ <eq:1_ar2>

Our implementation will be based on the `statsmodels` library in Python, which offers an implementation of the `SARIMAX` (Seasonal AutoRegressive Integrated Moving Average with eXogenous regressors) model.

Importantly, the `SARIMAX` model follows an alternative formulation of the model:
$
  y_t = φ_1 y_(t-1) + φ_2 y_(t-2) + ε_t
$ <eq:1_sarimax_statsmodels>

Where we have used the two different variants of the letter phi to denote the different formulations.
Comparing @eq:1_ar2 and @eq:1_sarimax_statsmodels we find that the two formulations are equivalent given a sign change:
$
  φ_1 = -ϕ_1, wide φ_2 = -ϕ_2
$ <eq:1_phi_sign_change>

== Realisations

=== 1.1 & 1.2
We simulate the given process 5 times using the ```SARIMAX``` module with $n=200$ observations and the coefficients set to: $ϕ_1 = -0.6$ and $ϕ_2 = 0.5$

#figure(
  // image("img/blabla.png"),
  [],
  caption: [Simulation and Autocorrelation function of AR(2) process with $ϕ_1 = -0.6, ϕ_2 = 0.5$.],
) <fig:1_1_sim-acf>

- plot Simulation and ACF together
- plot "empirical" ACF and rho(k) for k lag in {0, ..., 30}
- comment on the result

=== 1.3

- new simulations & ACFs with $ϕ_1 = -0.6$ and $ϕ_2 = -0.3$
- comment on each, focus on stationarity

=== 1.4

- with $ϕ_1 = 0.6$ and $ϕ_2 = -0.3$

=== 1.5

- with $ϕ_1 = -0.7$ and $ϕ_2 = -0.3$

=== 1.6

- with $ϕ_1 = -0.75$ and $ϕ_2 = -0.3$

=== 1.7

- comment on: is it enough to plot the ACF or is the time-series helpful?
- ...in @fig:1_1_sim-acf we can see...

#pagebreak()

= Predicting Monthly Solar Power <sec:2_predicting_monthly_solar_power>

We are given a seasonal AR model:
$
  (1 - ϕ_1B)(1 + Φ_1B^12)(log(Y_t)-µ) = ε_t
$ <eq:2_seasonal_ar>

where $Y_t$ is the monthly energy from the plant in MWh, ${ε_t}$ is a white-noise process with variance $σ_ε^2$ . The parameters $ϕ_1 = −0.38$, $Φ_1 = −0.94$ and $µ = 5.72$ are assumed to be known. Based on 36 observations, it is found that $σ_ε^2 = 0.22^2$.

== Forecasting and Residuals

=== Reformulation

Using a substitution $X_t = log(Y_t) - μ$ we can rewrite the model, @eq:2_seasonal_ar, as:
$
  (1 - ϕ_1B)(1 + Φ_1B^12)X_t = ε_t
$ <eq:2_rewritten_ar_1>

We are then asked to rewrite the model in order to calculate the residuals $hat(ε)_(t+1|t)$,
which we understand to be the estimated one-step-ahead forecast error, such that $hat(ε)_(t+1|t) ≔ X_(t+1) - hat(X)_(t+1)$.
Recalling $B^q X_t ≔ X_(t-q)$, we expand @eq:2_rewritten_ar_1 to find:
$
  (1 - ϕ_1B)(1 + Φ_1B^12)X_t &= ε_t\
  X_t - ϕ_1 B X_t + Φ_1 B^12 X_t - ϕ_1 Φ_1 B^13 X_t &= ε_t\
  X_t - ϕ_1 X_(t-1) + Φ_1 X_(t-12) - ϕ_1 Φ_1 X_(t-13) &= ε_t\
$

Multiplying through by $B^(-1)$ we find:
$
  X_(t+1) - ϕ_1 X_t + Φ_1 X_(t-11) - ϕ_1 Φ_1 X_(t-12) &= ε_(t+1)\
$

Solving for $X_(t+1)$ we find:
$
  X_(t+1) = ϕ_1 X_t - Φ_1 X_(t-11) + ϕ_1 Φ_1 X_(t-12) + ε_(t+1)\
$ <eq:2_forecast_1>

We then construct the _optimal predictor_ for $X_(t+1)$, @Madsen_2008[eq.~5.139]:
$
  hat(X)_(t+1)
  &= 𝔼[X_(t+1)|X_t, X_(t-1), …]\
  &= hat(ϕ)_1 X_t - hat(Φ)_1 X_(t-11) + hat(ϕ)_1 hat(Φ)_1 X_(t-12) wide wide 𝔼[ε] = 0\
$

Where we have used the fact that the $ε∼𝒩(0, σ^2_ε)$.
Thus the one-step-ahead forecast error becomes:
$
  hat(ε)_(t+1|t) &= X_(t+1) - hat(X)_(t+1)\
  &= X(t+1) - hat(ϕ)_1 X_t + hat(Φ)_1 X_(t-11) - hat(ϕ)_1 hat(Φ)_1 X_(t-12)\
$

=== Verifying Residuals in Dataset <sec:2_1_a>

As already noted in the section above, we assume that the residuals are independent and identically distributed (i.i.d.)
with a zero-mean normal distribution.

We can verify that this is indeed the case by constructing the
appropropriate seasonal AR model using the `SARIMAX` class from `statsmodels`.
Again we take careful notice of the signs of the parameters $ϕ_1$, $Φ_1$.
Rather strangely, these are given with two different conventions, so only $ϕ_1$
must be negated to match the `SARIMAX` implementation.

Utilising the built-in `plot_diagonstics` function of the fitted model, we produce @fig:2_residual_analysis.

#figure(
  image("output/2_residual_analysis.png"),
  caption: [Residual analysis of the solar power dataset using $ϕ_1=-0.38$, $Φ_1=-0.94$ and $σ_ε^2=0.22^2$
  in the notation established by @eq:2_rewritten_ar_1.
  Note that the residual analysis is performed on the transformed data $X_t = log(Y_t) - μ$.],
) <fig:2_residual_analysis>

By inspection of @fig:2_residual_analysis we find no evidence that the residuals violate
the assumption of i.i.d. with a zero-mean normal distribution.

== Forecasting

Using the same model that we fitted in @sec:2_1_a, we produce a 12 month forecast ($k=12$)
using the `forecast` method.

Importantly, the forecasted data must be transformed back into the original scale of the dataset:
$
  X_t & := log(Y_t) - μ wide wide &"Forward"\
  Y_t & = e^(X_t + μ) wide wide &"Backward"\
$

#let data_2_2= csv("output/2_2.csv")

#figure(
  table(
    columns: data_2_2.at(0).len(),
    table.header([*Time*], [*Power, $bold(Y_t)$ \[MWh*\]], $bold(X_t)$),
    ..data_2_2.slice(1).flatten(),
  ),
  caption: [Twelve month forecast using model given in @eq:2_seasonal_ar.]
) <table:2_2_forecast>

Which we can plot alongside the original data to produce @fig:2_2_forecast.
#figure(
  image("output/2_2_forecast.png"),
  caption: [Twelve month forecast using model given in @eq:2_seasonal_ar.]
) <fig:2_2_forecast>

We find that the forecasted data matches our intuition of the data well and
follows the seasonal pattern observed in the historical data of power generation.

== Prediction Intervals

#bibliography("report.bib")
