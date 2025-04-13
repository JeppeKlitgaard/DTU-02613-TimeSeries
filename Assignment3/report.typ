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
#let Var = math.op("Var")

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
  (1 + ϕ_1B)(1 + Φ_1B^12)(log(Y_t)-µ) = ε_t
$ <eq:2_seasonal_ar>

where $Y_t$ is the monthly energy from the plant in MWh, ${ε_t}$ is a white-noise process with variance $σ_ε^2$ . The parameters $ϕ_1 = −0.38$, $Φ_1 = −0.94$ and $µ = 5.72$ are assumed to be known. Based on 36 observations, it is found that $σ_ε^2 = 0.22^2$.

== Forecasting and Residuals

=== Reformulation

Using a substitution $X_t = log(Y_t) - μ$ we can rewrite the model, @eq:2_seasonal_ar, as:
$
  (1 + ϕ_1B)(1 + Φ_1B^12)X_t = ε_t
$ <eq:2_rewritten_ar_1>

We are then asked to rewrite the model in order to calculate the residuals $hat(ε)_(t+1|t)$,
which we understand to be the estimated one-step-ahead forecast error, such that $hat(ε)_(t+1|t) ≔ X_(t+1) - hat(X)_(t+1)$.
Recalling $B^q X_t ≔ X_(t-q)$, we expand @eq:2_rewritten_ar_1 to find:
$
  (1 + ϕ_1B)(1 + Φ_1B^12)X_t &= ε_t\
  X_t + ϕ_1 B X_t + Φ_1 B^12 X_t + ϕ_1 Φ_1 B^13 X_t &= ε_t\
  X_t + ϕ_1 X_(t-1) + Φ_1 X_(t-12) + ϕ_1 Φ_1 X_(t-13) &= ε_t\
$

Multiplying through by $B^(-1)$ we find:
$
  X_(t+1) + ϕ_1 X_t + Φ_1 X_(t-11) + ϕ_1 Φ_1 X_(t-12) &= ε_(t+1)\
$

Solving for $X_(t+1)$ we find:
$
  X_(t+1) = ε_(t+1) - ϕ_1 X_t - Φ_1 X_(t-11) - ϕ_1 Φ_1 X_(t-12)\
$ <eq:2_forecast_1>

We then construct the _optimal predictor_ for $X_(t+1)$, @Madsen_2008[eq.~5.139]:
$
  hat(X)_(t+1)
  &= 𝔼[X_(t+1)|X_t, X_(t-1), …]\
  &= -hat(ϕ)_1 X_t - hat(Φ)_1 X_(t-11) - hat(ϕ)_1 hat(Φ)_1 X_(t-12) wide wide 𝔼[ε] = 0\
$

Where we have used the fact that the $ε∼𝒩(0, σ^2_ε)$.
Thus the one-step-ahead forecast error becomes:
$
  hat(ε)_(t+1|t) &= X_(t+1) - hat(X)_(t+1)\
  &= X_(t+1) + hat(ϕ)_1 X_t + hat(Φ)_1 X_(t-11) + hat(ϕ)_1 hat(Φ)_1 X_(t-12)\
$ <eq:2_1_prediction_error>

=== Verifying Residuals in Dataset <sec:2_1_a>

As already noted in the section above, we assume that the residuals are independent and identically distributed (i.i.d.)
with a zero-mean normal distribution.

We can verify that this is indeed the case by constructing the
appropropriate seasonal AR model using the `SARIMAX` class from `statsmodels`.
Again we take careful notice of the signs of the parameters $ϕ_1$, $Φ_1$, both
of which must be negated to match the convention of the `SARIMAX` class.

Utilising the built-in `plot_diagonstics` function of the fitted model, we produce @fig:2_residual_analysis.

#figure(
  image("output/2_residual_analysis.png"),
  caption: [Residual analysis of the solar power dataset using $ϕ_1=-0.38$, $Φ_1=-0.94$ and $σ_ε^2=0.22^2$
  in the notation established by @eq:2_rewritten_ar_1.
  Note that the residual analysis is performed on the transformed data $X_t = log(Y_t) - μ$.],
) <fig:2_residual_analysis>

By inspection of @fig:2_residual_analysis we find no evidence that the residuals violate
the assumption of i.i.d. with a zero-mean normal distribution.
In particular, we find that the standardised residuals do not exhibit strong systematic patterns,
although the first summer (observations 5-10) does appear to be skewed slightly positive in the residuals.
Additionally, we find the residuals to be well-described by a normal distribution as seen in the Q-Q plot and histogram.

One could carry out a Lljung-Box test to further substantiate this claim, but we find that
this is not nececssary in this case.

== Forecasting

Using the same model that we fitted in @sec:2_1_a, we produce a 12 month forecast ($k=12$)
using the `forecast` method.

Importantly, the forecasted data must be transformed back into the original scale of the dataset:
$
  X_t & := log(Y_t) - μ wide wide &"Forward"\
  Y_t & = e^(X_t + μ) wide wide &"Backward"\
$ <eq:2_transformations>

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

We now wish to consider the _variance_ in our _predictions_, or rather we seek to
obtain the _variance of the prediction error_, @eq:2_1_prediction_error:
$
  Var(hat(ε)_(t+1|t)) = Var(X_(t+1) - hat(X)_(t+1))\
$

In order to carry out this calculation, we convert our process to a Moving-Average (MA) process:
$
  X_t
  &= ε_t/((1 + ϕ_1B)(1 + Φ_1B^12)) wide wide ⇔ &&#mref("eq:2_seasonal_ar")\
  &= (1 + ψ_1 B + ψ_2 B^2 + …)ε_t+ ... wide wide && "Definition of MA process"\
$ <eq:2_ma_process_1>

From @eq:2_ma_process_1 we recover the constraint:
$
  (1 + ψ_1 B + ψ_2 B^2 + …)(1 + ϕ_1 B)(1 + Φ_1 B^12) &= 1\
  (1 + ψ_1 B + ψ_2 B^2 + …)(1 + ϕ_1 B + Φ_1 B^12 + ϕ_1 Φ_1 B^13) &=1\
$ <eq:2_ma_process_constraints>

Which imposes that all coefficients of the power series of $B$ must vanish.
Comparing powers of $B^i$ we find the following constraints on the weights
$ψ_i$ of the transfer function:
$
  ψ_i = cases(
    (-ϕ_1)^i wide & 1 ≤ i ≤ 11\
    (-ϕ_1)^i - Φ_1(-ϕ_1)^(i-12) wide & 12 ≤ i ≤ 23\
    ⋮ wide & ⋮
  )
$ <eq:2_ma_weights>

Using our MA formulation we again consider the expression of $X_(t+1)$ using the treatment in @Madsen_2008[sec:~5.7.1]:
$
  X_(t+k) = ε_(t+k) + ψ_1 ε_(t+k-1) + … + ψ_k ε_(t) + ψ_(k+1) ε_(t-1) + …
$

We consider the expectation value of the residual $ε_(t+k)$.
For past observations, when $k ≤ 0$, have a realisation of the residual and thus this will be our expected value.
For future observations, when $k > 0$,
we have no realisation of the residual and thus the expected value is the mean of the distribution of $ε$,
which we recall to be zero:
$
  𝔼[ε_(t+k)|X_t, X_(t-1), …] = cases(
    ε_(t+k) wide & k ≤ 0\
    0 wide & k > 0
  )
$ <eq:2_residual_expectation>

From which we can obtain our $k$-step predictions, $X_(t+k|t)$:
$
  X_(t+k|t)
  &= 𝔼[X_(t+k)|X_t, X_(t-1), ...]\
  &= cancel(𝔼[ε_(t+k)|X_t, X_(t-1), ...] + ψ_1 𝔼[ε_(t+k-1)|X_t, X_(t-1), ...] + …) wide wide && ⇐ #mref("eq:2_residual_expectation")\
  &quad + ψ_k 𝔼[ε_(t)|X_t, X_(t-1), ...] + ψ_(k+1) 𝔼[ε_(t-1)|X_t, X_(t-1), ...] + …\
  &= ψ_k ε_(t) + ψ_(k+1) ε_(t-1) + …\

$

Using our previous definition of the prediction error, @eq:2_1_prediction_error, we find:
$
  hat(ε)_(t+k|t)
    &= X_(t+k) - hat(X)_(t+k)\
    &= ε_(t+k) + ψ_1 ε_(t+k-1) + … + cancel(ψ_k ε_(t) + ψ_(k+1) ε_(t-1) + …)\
    &quad cancel(-(ψ_k ε_(t) + ψ_(k+1) ε_(t-1) + …))\
    &= ε_(t+k) + ψ_1 ε_(t+k-1) + … + ψ_(k-1) ε_(t+1)
$ <eq:2_prediction_error_k>

Where the variance of this prediction error at time $t+k$ is given by:
$
  σ^2_(k)
  &= Var(hat(ε)_(t+k|t))\
  &= Var(ε_(t+k) + ψ_1 ε_(t+k-1) + … + ψ_(k-1) ε_(t+1))\
  &= Var(ε_(t+k)) + ψ_1^2 Var(ε_(t+k-1)) + … + ψ_(k-1)^2 Var(ε_(t+1))\
$

We consider that $hat(ε)_(t+k|t) ∼ 𝒩(0, σ^2_ε)$ and thus find $Var(ε_(t+k)) = σ^2_ε$ for all $k>0$.
Plugging this into @eq:2_prediction_error_k yields:
$
  σ^2_(k) = σ^2_ε (1 + ψ_1^2 + … + ψ_(k-1)^2)
$

This enables us to calculate the $(1-α)$ confidence interval of our $k$-step predictions:
$
  hat(X)_(t+k|t) ± u_(α\/2) σ_(k)
  = hat(X)_(t+k|t) ± u_(α\/2) σ_ε sqrt(1 + ψ_1^2 + … + ψ_(k-1)^2)\
$ <eq:2_prediction_interval>

Where $u_(α\/2)$ is the $α\/2$ quantile of the standard normal distribution.
For $α = 5%$ we find $u_(α\/2) ≈ 1.96$.

With this, we are able to calculate the prediction intervals for our forecasted data
using by leveraging the expression for the weights $ψ_i$ in @eq:2_ma_weights,
which yields the 95% prediction intervals shown in @fig:2_3_forecast_ci.

#figure(
  image("output/2_3_forecast_ci.png"),
  caption: [Twelve month forecast using model given in @eq:2_seasonal_ar with 95% prediction intervals.]
) <fig:2_3_forecast_ci>

We note that the for $1 <= k <= 11$ the seasonal component does not contribute to the prediction error
as also hinted in the assignment description,
but for $k=12$ we find a linear dependence on the coefficient of the seasonal component. $Φ_1$.
This is a significant difference when contrasted against not including the seasonal component.


== Forecast Commentary

We consider the forecast found in @fig:2_3_forecast_ci to be a relatively accurate description
of the expected power generation of the solar plant during the year of 2011.
In the original dataset we observe a slight decline in the peak power generation as a function of time,
which is also observed in the forecast.
This can be attributed to expected degradation of the photo-voltaic solar panels over time.

Notably, the prediction intervals are not symmetric over the expected value,
which is a result of the non-linear transformation of the data as described in @eq:2_transformations.

Additionally, we observe that the prediction intervals are quite wide,
which again can be understood in part by referring to the exponential transformation of the data.
Additionally we observe large variations in the shape of the observations of the 3 year history of the solar plant,
which pushes up the estimated variance of the residuals, $σ_ε^2 = 0.22^2$.
This strongly influences the width of the prediction intervals, which is particularly notable
at power generation observations due to the transformation.
This is simply a consequence of linearising the data in order to fit it to our linear model.

From inspection, we find that the prediction intervals are likely too wide during the summer months
where power generation is at its peak, due to overly large effect of the variance during the winter months on the residual error estimates.

The implementation of the prediction intervals can be found in the attached Jupyter Notebook `JK_2.ipynb`.

#bibliography("report.bib")
