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
#set text(font: "STIX Two Text", size: 10pt)
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

== <sec:1_1>
We simulate 5 realisations of the process up to $n=200$ observations using a burn-in period of $N_B = 10000$ with parameters $ϕ_1 = -0.6$ and $ϕ_2 = 0.5$ and plot the results in @fig:1_1.

#figure(
  image("output/1_1.png"),
  caption: [Simulation and Autocorrelation function of AR(2) process with $ϕ_1 = -0.6, ϕ_2 = 0.5$.],
) <fig:1_1>

We find that the results look as expected.

== <sec:1_2>

We recall the definition of _autocorrelation_ for a _stationary process_:
$
  ρ(k) = γ(k)/γ(0)
$ <eq:1_autocorrelation>

Where $γ(k)$ is the _autocovariance_ for a timeshift $k$:
$
  γ(k) = Cov[X_t, X_(t+k)]
$ <eq:1_autocovariance>

We consider the autocorrelation of an AR(2) process by solving for $X_t$ in @eq:1_ar2:
$
  X_t = -ϕ_1 X_(t-1) - ϕ_2 X_(t-2) + ε_t
$ <eq:1_ar2_solved>

And then inserting this in @eq:1_autocorrelation:
$
  ρ(k) = Cov[X_t, X_(t+k)] / γ(0)
    &= Cov[X_t, -ϕ_1 X_(t-1+k) - ϕ_2 X_(t-2+k) + ε_t] / γ(0)\
    &= Cov[X_t, -ϕ_1 X_(t-1+k)] / γ(0) + Cov[X_t, -ϕ_2 X_(t-2+k)] / γ(0) + cancel(Cov[X_t, ε_t] / γ(0))\
    &= -ϕ_1 Cov[X_t, X_(t-1+k)] / γ(0) - ϕ_2 Cov[X_t, X_(t-2+k)] / γ(0)\
    &= -ϕ_1 γ(k-1) / γ(0) - ϕ_2 γ(k-2) / γ(0)\
  ρ(k) &= -ϕ_1 ρ(k-1) - ϕ_2 ρ(k-2)\
$ <eq:1_autocorrelation_recursion>

Notably, by stationary it follows $ρ(-k) = ρ(k)$ and from @eq:1_autocorrelation we find $ρ(0) = 1$,
which allows us to build a recursive relation for $ρ(k)$ from:

$
  ρ(0) &= 1\
  ρ(1) &= -ϕ_1 ρ(0) - ϕ_2 ρ(-1) = -ϕ_1 - ϕ_2ρ(1) = (-ϕ_1)/(1+ϕ_2)\
$ <eq:1_autocorrelation_recursion_2>

We can compute the empirical autocorrelation function (ACF) using the `plot_acf` function from `statsmodels`,
which we compute for each of the realisations shown in @fig:1_1.

This is plotted in @fig:1_2_acf along with the theoretical autocorrelation function, $ρ(k)$, as given by @eq:1_autocorrelation_recursion_2.

#figure(
  image("output/1_2_acf.png"),
  caption: [
    Theoretical and empirical autocorrelation functions of an AR(2) process with $ϕ_1 = -0.6, ϕ_2 = 0.5$.
    Shaded blue region is the 95% confidence interval for the empirical ACF.
  ],
) <fig:1_2_acf>

We observe good agreement between the theoretical and empirical autocorrelation functions.
The deviations observed after $k=4$ are results of the finite sample size of the simulation,
which is also reflected in the shaded blue confidence intervals.
The damped oscillations of the autocorrelation functions are expected, as also indicated
in @table:order_identification on #ref(<table:order_identification>, form: "page", supplement: "Page").


#pagebreak()
== <sec:1_3>

We now repeat the simulation with $ϕ_1 = -0.6$ and $ϕ_2 = -0.3$ and plot the results in @fig:1_3 and @fig:1_3_acf.

#figure(
  image("output/1_3.png"),
  caption: [Simulation and Autocorrelation function of AR(2) process with $ϕ_1 = -0.6, ϕ_2 = -0.3$.],
) <fig:1_3>

#figure(
  image("output/1_3_acf.png"),
  caption: [
    Theoretical and empirical autocorrelation functions of an AR(2) process with $ϕ_1 = -0.6, ϕ_2 = -0.3$.
    Shaded blue region is the 95% confidence interval for the empirical ACF.
  ],
) <fig:1_3_acf>

We observe that the process remains stationary, but the characteristic correlation times
are now much longer than for the previous case in @sec:1_1 and @sec:1_2.
In @fig:1_3 we can clearly see that variations in the time series are much more pronounced than in @fig:1_1.

#pagebreak()
== <sec:1_4>

We now repeat the simulation with $ϕ_1 = 0.6$ and $ϕ_2 = -0.3$ and plot the results in @fig:1_4 and @fig:1_4_acf.

#figure(
  image("output/1_4.png"),
  caption: [Simulation and Autocorrelation function of AR(2) process with $ϕ_1 = 0.6, ϕ_2 = -0.3$.],
) <fig:1_4>

#figure(
  image("output/1_4_acf.png"),
  caption: [
    Theoretical and empirical autocorrelation functions of an AR(2) process with $ϕ_1 = 0.6, ϕ_2 = -0.3$.
    Shaded blue region is the 95% confidence interval for the empirical ACF.
  ],
) <fig:1_4_acf>

We again observe that the process is stationary, but find that the sign change in $ϕ_2$ with respect to @sec:1_3 has introduced a fast oscillation in the time series.

We find that for realisations 1 and 3 the autocorrelation does decays slower than expected by
theory. This is again attributed to the finite sample size of the simulation,
the artifacts of which become more pronounced for processes with longer characteristic
correlation times such as the one observed in @fig:1_4.

We propose increasing the number of simulated observations, $n$, to a larger number
for such processes if the empirical autocorrelation function is to be used for identification
or interpretation.

#pagebreak()
== <sec:1_5>

We again repeat the simulation with $ϕ_1 = -0.7$ and $ϕ_2 = -0.3$ and plot the results in @fig:1_5 and @fig:1_5_acf.

#figure(
  image("output/1_5.png"),
  caption: [Simulation and Autocorrelation function of AR(2) process with $ϕ_1 = -0.7, ϕ_2 = -0.3$.],
) <fig:1_5>

#figure(
  image("output/1_5_acf.png"),
  caption: [
    Theoretical and empirical autocorrelation functions of an AR(2) process with $ϕ_1 = -0.7, ϕ_2 = -0.3$.
    Shaded blue region is the 95% confidence interval for the empirical ACF.
  ],
) <fig:1_5_acf>

In figure @fig:1_5 we observe that the process is exactly at the critical point
of stationarity, where the theoretical autocorrelation function does not decay to zero.

This can be understood by considering the stationarity conditions for an AR(2) process, which
require that the _roots_, $z_(1,2)$ of the _characteristic polynomial_ lie outside the unit circle, such that $|z_(1, 2)| > 1$.

For the AR(2) process in @eq:1_ar2 the characteristic polynomial is given by:
$
  P(z) = 1 + ϕ_1 z + ϕ_2 z^2
$

Where the roots are then given by:
$
  z_(1,2) &= (ϕ_1 ± sqrt(ϕ_1^2 + 4ϕ_2))/(-2 ϕ_2)\
$ <eq:ar2_roots>

We find:
$
  |z_(1,2)| = {1, 3.overline(33)}
$

Where clearly, $|z_1| = 1 ≯ 1$, violates the strict inequality of the stationarity condition,
the consequence of which is a process where shocks do not decay away over time.

We note that the process does not explode because we are exactly at the point of criticality.

#pagebreak()
== <sec:1_6>

Once more we repeat the simulation, this time with $ϕ_1 = -0.75$ and $ϕ_2 = -0.3$ and plot the results in @fig:1_6 and @fig:1_6_acf.

#figure(
  image("output/1_6.png"),
  caption: [Simulation and Autocorrelation function of AR(2) process with $ϕ_1 = -0.75, ϕ_2 = -0.3$.],
) <fig:1_6>

#figure(
  image("output/1_6_acf.png"),
  caption: [
    Theoretical and empirical autocorrelation functions of an AR(2) process with $ϕ_1 = -0.75, ϕ_2 = -0.3$.
  ],
) <fig:1_6_acf>

Using @eq:ar2_roots we find:
$
    |z_(1,2)| ≈ {0.963, 3.463}
$

Which reveals that the process is decidedly not stationary and should blow up
as time progresses, which is also what we observe in @fig:1_6 and @fig:1_6_acf.

@fig:1_6_acf fails to compute the empirical ACF, while the theoretical ACF quickly increases with time.

== <sec:1_7>

While a lot of information may be gleaned from the autocorrelation function alone, it is useful to also plot the time series itself to get a better understanding of the process, particularly
for processes with long characteristic correlation times.
For instance, the autocorrelation does not give information about the scale amplitude of the process.

For higher order AR processes, the autocorrelation function may be difficult to interpret,
in which case it is again useful to plot the time series itself.

Lastly, the autocorrelation function provides no information about the process realisations
at specific points in time, which is often useful for understanding or interpreting the underlying process.

The code used for the simulations can be found in the attached notebook `JK_1.ipynb`.

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

One could carry out a Ljung-Box test to further substantiate this claim, but we find that
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
    align: (center, right, right),
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
where power generation is at its peak,
due to the overly large effect of the variance during the winter months on the residual error estimates.

The implementation of the prediction intervals can be found in the attached Jupyter Notebook `JK_2.ipynb`.

#pagebreak()
= An ARX Model for the Heating of a Box <sec:3>

We are given a data set of hourly measurements from a box instrumented with a heater,
an internal temperature sensor and an external temperature sensor.

The box has window on its south-facing wall, onto which the vertical solar radiation is also recorded in the data set.

The internal temperature of the box was kept approximately constant using a thermostatic control of the internal heater.

In this section, we will be working AutoRegressive models with eXogenous variables (ARX),
for which we will adopt a notation in which $upright("AR")(p)"-"upright("X") (e_1, e_2, …)$ refers to
an ARX where $p$ determines the order of the AutoRegressive (AR) part of the model while $e_j$ for $j ∈ ℕ_+$ refers to the order of the $j$th exogenous variable $X_j$.

This yields the following generic ARX model:

$
  Y_t
    = c
    - underbrace(sum_(i = 1)^p ϕ_i Y_(t-i), "AR Part")
    + underbrace(sum_(j = 1)^J sum_(k=0)^(e_j) ω_(j, k) X_(j, t-k), "Exogenous Part")
    + ε_t
$

Where $c$ is a constant offset, $ϕ_i$ are the AR coefficients, $ω_(j, k)$ are the exogenous coefficients and $ε_t$ is a white-noise such that $ε_t ∼ 𝒩(0, σ_ε^2)$. $J$ refers to the number of exogenous variables such that $j ∈ {1, 2, …, J}$ become the indices of the exogenous variables.

To further unify notation, we assume a model, where ${ T_t }$ represents the series $T_(upright("delta"))$, ${ G_t }$ is the series $G_v$ and ${P_t }$ is the series $P_h$ for $t in { 1 \, . . . \, 167 }$ in the training dataset.

== Exploratory Dependency Analysis <sec:3_1>

We load the data set and seek to visually explore the dependency between the heater power $P_h$, temperature delta between inside and outside, $Δ T$, and the vertical solar radiation $G_v$.

For the purposes of visualisation, we standardize the data by subtracting the mean and dividing by the standard deviation,
which we then plot in @fig:3_1_analysis.

#figure(
  image("output/3_1_analysis.png"),
  caption: [Standardized data set of heater power $P_h$, temperature delta $Δ T$ and vertical solar radiation $G_v$.]
) <fig:3_1_analysis>

@fig:3_1_analysis clearly shows an inverse relationship between the heater power and solar radiation,
as we would intuitively expect on physical grounds – for a constant temperature difference, we would expect the sum of the heater power and solar radiation to be constant, thus requiring the heater to compensate for any changes in solar heating.

We additionally observe a positive correlation between the heater power and the temperature difference,
which again matches our physical understanding that a larger temperature difference would lead to a larger heat loss from the box, assuming the inside temperature is larger than the external temperature. Thermostaticity would then require a larger heater power to compensate for the increase in heat loss.

== Test/Train Split <sec:3_2>

We now split the dataset into two parts, one for training and one for testing.
The cut-off point is set to 2013-02-06 00:00, which is the last observation of the training set.

== Further Dependency Investigation <sec:3_3>

Following on from the investigation in @sec:3_1,
we plot the three variables against each other in the pair plots found in @fig:3_3_pairplot.

#figure(
  image("output/3_3_pairplot.png"),
  caption: [Pair plot of the training data set.]
) <fig:3_3_pairplot>

In @fig:3_3_pairplot we again observe that the heater power $P_h$ and the temperature difference $Δ T$ are positively correlated, with the heater power being larger when the temperature difference is larger.
Interestingly, we find that the temperature difference as a function of the heater power roughly follows
one of two lines, which we can further understand when inspecting the pair plot between the temperature difference and the solar radiation $G_v$, which appears to be roughly bimodal.
This does not lend itself to immediate interpretation, but does suggest that there is an underlying structure in
the data that may be elucidated with further analysis.

Additionally, we confirm the observation from @sec:3_1 that the heater power and solar radiation $G_v$ are negatively correlated.

Additionally we plot the autocorrelations and cross-correlations in a similar matrix to the one found in @fig:3_3_pairplot,
which reveals characteristic timescales for the three variables and their correlations.

#figure(
  image("output/3_3_acf_ccf.png"),
  caption: [Autocorrelation and cross-correlation of the training data set.]
) <fig:3_3_acf_ccf>

In @fig:3_3_acf_ccf we find a seasonality of approximately 24 hours in the heater power $P_h$
and solar radiation $G_v$, which we would expect from the diurnal solar cycle.

The temperature difference does not follow this pattern, partly because it has a much longer timescale in the data set.
This could be attributed to more complex weather patterns, such as different wind patterns or nightly cloud cover.
The characteristic timescale of the autocorrelation function of the temperature difference
reveals that the external temperature varies rather slowly compared to the other two variables.

The inverse relationship between solar radiation and heater power is also confirmed in the cross-correlation plot,
where we find a negative correlation between the two variables.
We find only a small correlation between the temperature difference and the solar radiation.

Lastly, we inspect the partial autocorrelation functions of the three variables in @fig:3_3_pacf.

#figure(
  image("output/3_3_pacf.png"),
  caption: [Partial autocorrelation of the training data set.]
) <fig:3_3_pacf>

Comparing @fig:3_3_pacf with the ARMA order identification table, @table:order_identification, we
find that the heater power $P_h$ and temperature difference $Δ T$ may be well described
by an AR(2) process, while the solar radiation $G_v$ appears to be well described by an AR(1) process.

The heater power $P_h$ and solar radiation $G_v$ additionally show evidence of seasonality, suggesting
a seasonality of 24 hours may be appropropriate for modelling their behaviour.

#figure(
  table(
    columns: (1.5fr, 4fr, 4fr),
    table.header("", [ACF $ρ(k)$], [PACF $φ_(k k)$]),
    $"AR"(p)$, [Damped exponential and/or sine functions], $φ_(k k) = 0 "for" k > p$,
    $"MA"(q)$, $ρ(k) = 0 "for" k > q$, [Dominated by damped exponential and/or sine functions],
    $"ARMA"(p,q)$, [Damped exponential and/or sine functions after lag $q - p$], [Dominated by damped exponential and/or sine functions after lag $p - q$s],
  ),
  caption: [Reproduction of Table 6.1 from @Madsen_2008[p.~155] showing the expected behaviour of the autocorrelation function for different ARMA processes],
) <table:order_identification>

Notably, we cannot reasonably estimate the parameters of, for example, an ARX model from the figures above,
will instead require an outright fitting of the model(s) to the data.

== Impulse Response <sec:3_4_impulse_response>

While the assignment description does not explicitly state the orders of the ARX model,
which naturally need to be determined prior to carrying out an impulse response analysis,
we consider an $"AR"(1)"-"X(1, 1)$ model to be a reasonable choice, which is supported by the analysis of the PACF in @fig:3_3_pacf above.

We consider the heater power, $P_h$, to be the _endogenous variable_, while the solar radiation, $G_v$, and temperature difference, $Δ T$, are the _exogenous variables_, which yields the following model:
$
  P_t = c - ϕ_1 P_(t-1) + ω_1 T_t + ω_2 G_t + ε_t\
$

Where parameters $c, ϕ_1, ω_1, ω_2 ∈ ℝ$, noting that we have also used the convention describing the different time series described in the introduction in @sec:3.

While the phrasing in the assignment description is somewhat ambiguous, we interpret
a "lag up to 10" to mean the that the impulse response function should be calculated
for lags up to and including 10 time steps, which importantly is different to the _lag_ given by
the order of the AR(p) process.

TODO: FINISH REFACTOR HERE

There are different options for modelling an impulse response for a model with 2 exogenous variables: $P_t$ with $G_t$ or $G_t$ as exogenous variable, as well as a model for $P_t$ with both as exogenous variables.

TODO: I DON'T REALLY UNDERSTAND THIS SENTENCE ↑

$
upright("AR(p)-X(e1, e2)") = c + sum_(i = 0)^p phi.alt_i B^i P_t + sum_(i = 0)^(e_1) omega_i B^i T_t + sum_(i = 0)^(e_2) beta_i B^i G_t + epsilon_t
$<eq:3_4_context_arx>

for both as exogenous, given parameters/coefficients $phi.alt_i \, beta_i \, omega_i in bb(R)$.

Now there are different options how to interpret the phrasing "up to lag 10". This could refer to the lagged samples of the modeled series $P_t$, hence the AR components, or it could refer to the depth of recursion for the impulse response function $upright("IRF")(k)$.
It is ambiguous from which of the exogenous variables, $G_t$ or $T_t$ the unit impulse should come from. Either from both together or each gives an impulse separately.

Therefore, we need to make reasonable choices for the model. We choose an AR(1)-X(1,1) model (based on the findings in @sec:3_3) with a given lag of of $k=10$ for the impulse response function IRF(k).

While we interpret the task, such that we model two impulse responses, from $G_t$ and $T_t$ separately, the unit impulse given from the variable $T_t$ seems much more meaningful. As @fig:3_3_pairplot shows, it has a positive correlation with the target series $P_t$, while @fig:3_1_analysis clearly indicates, the variable $G_t$ has an inverse relationship with the target $P_t$.
Even before estimating anything, we can assume that a unit impulse from $T_t$ will be impactful, while for $G_t$ the response will decay extremely fast, as the original series have opposing effects.

Having such a relatively simple model, one could potentially read-off maybe the first $phi.alt_1$ parameter of the AR part via an ACF and PACF plot to conclude that probably $|phi.alt_1|<1$. Maybe even find an argument for alternating signs throughout all coefficients, based on the PACF plot.

In summary, the final AR(1)-X(1,1) model becomes:

$
P_t = - phi.alt_1 P_(t-1) + omega_1 T_t + beta_1 G_t + epsilon_t
$<eq:3_4_ar1_x1_1>

We will use the the ```AutoReg```class from ```statsmodels.tsa.ar_model``` in Python, which allows to provide an exogenous variable and has built in parameter estimation (via OLS and conditional MLE).
To double check, we also did our own classic OLS fit on the design matrix: $X = [P_(t-1), T_t, G_t]$ noted in the form of column vectors of the corresponding series, which yields a parameter vector of $Theta = [- phi.alt_1, omega_1, beta_1]^T in bb(R)^(p+e_1+e_2)$.

$
  arrow.r.double &  & P_(t-p, n) & = X dot.op Theta + epsilon_t = Y \
  arrow.l.r.double &  & hat(Theta) & = (X^T X)^(- 1) X^T Y \
$<eq:3_4_ols_parameter_est>

The resulting parameters are given as:

#figure(
  table(
    columns: 3,
    table.header($phi.alt_1$, $omega_1$, $beta_1$),
    [0.4127], [2.3352], [-0.0836],
  ),
  caption: [AR(1)-X(1,1) coefficients with $T_d$ and $G_v$ exogenous of order 1],
) <table:3_4_irf_arx_coeff>

For the impulse response function, we give the unit impulse first from $T_t$ and use the recursion

$
upright("IR") (0) & = 0 phi.alt_1+ 1 omega_1 + 0 beta_1\
upright("IR") (1) & = upright("IR") (0) phi.alt_1 + upright("IR") (0) omega_1 + upright("IR") (0) beta_1 \
upright("IR") (2) & = upright("IR") (1) phi.alt_1 + upright("IR") (0) phi.alt_2 + 0 omega_1 + 0 beta_1 \
upright("IR") (3) & = upright("IR") (2) phi.alt_1 + upright("IR") (1) phi.alt_2 + upright("IR") (0) phi.alt_3 + 0 omega_1 + 0 beta_1\
 & #h(0em) #h(0em) dots.v\
upright("IR") (k = 15) & = upright("IR") (14) phi.alt_1 + upright("IR") (13) phi.alt_2 + ... + upright("IR") (5) phi.alt_10 + 0 omega_1 + 0 beta_1\
 & #h(0em) #h(0em) dots.v\
$<eq:3_4_ir_recursion>

Or in shorter sum notation this would be with $upright("IRF")(0) = 1 omega_1 + 0 beta_1$ set as initial value:

$
  upright("IRF")(k) = sum_(i=1)^k sum_(j=0)^(min(i,p)) [ phi.alt_j upright("IRF") (i-j-1) bb(I)_(i-j-1≥0) ]
$<eq:3_4_ir_summation>

for $k in bb(N)$ the lag and $p$ the AR model order.

#figure(
  image("output/3_4_ar1_x1_1_impulse_response_Td.png"),
  caption: [Impulse Response of the AR(1)-X(1,1) with the unit impulse from $T_d$]
) <fig:3_4_ir_Tt>

Now we give a second plot of a unit impulse from $G_v$:

#figure(
  image("output/3_4_ar1_x1_1_impulse_response_Gv.png"),
  caption: [Impulse Response of the AR(1)-X(1,1) with the unit impulse from $G_v$]
) <fig:3_3_ir_Gv>

The impulse response coefficients $v_k$ are:

#figure(
  table(
    columns: 11,
    table.header("impulse", $v_1$, $v_2$, $v_3$, $v_4$, $v_5$, $v_6$, $v_7$, $v_8$, $v_9$, $v_10$,),
    $T_t arrow.r P_t$, [2.3352], [-0.9637], [0.3977], [-0.1641], [0.0677], [-0.0279], [0.0115], [-0.0048], [0.002], [-0.0008],
    $G_t arrow.r P_t$, [-0.0836], [0.0345], [-0.0142], [0.0059], [-0.0024], [0.001], [-0.0004], [0.0002], [-0.0001], [0.0]
  ),
  caption: [impulse response coefficients for AR(1)-X(1,1)],
) <table:3_4_irf_coeff>

As already explained based on @fig:3_1_analysis, the two impulse responses counter-act. The IR (impulse response) for this model decays rather quickly, especially with the impact of an exogenous variable that is reciprocal to the target variable. Most of these effects were already assumed based on available information.

We can also confirm, that this process is stable, since $sum_(k=1)^(infinity) |v_k| < infinity$.

Because of this stability the model is rather robust to unit shocks from both exogenous variables, more robust from shocks of $G_v$ than from $T_t$ in fact. As a consequence, with this model, we can accept a certain degree of fluctuation in the exogenous variables, without having to worry about a big impact on the predictive performance of our model. For example, the effect of a measurements error or sensor defect in $T_d$ or $G_v$ that acts as a shock to our system, will not last for very long.
Thus, having a larger order for the AR part of our model pretects against shocks from the other variables. On the other hand, the model is more vulnerable to sudden changes in auto-regressive behaviour.

== Linear Regression Model <sec:3_5>

We are given a linear model with the following form:
$
  P_(h,t) = ω_1 T_("delta",t) + ω_2 G_(v,t) + ε_t\
$

Where $ε_t ~ 𝒩(0, σ^2)$ and assumed to be i..i.d. (independent and identically distributed).

We fit such a model to the data and produce the estimation shown in @fig:3_5_forecast.

#figure(
  image("output/3_5_forecast.png"),
  caption: [Linear regression model with $T_"delta"$ and $G_v$ as exogenous variables.]
) <fig:3_5_forecast>

Additionally we evaluate the residuals during the training period and one-step prediction erorrs
for the test data set in @fig:3_5_errors.

#figure(
  image("output/3_5_errors.png"),
  caption: [Residual analysis of the linear regression model.]
) <fig:3_5_errors>

We observe that the errors, here understood to be the residuals for the training data set and the one-step prediction errors for the test data set, look relatively good.
We have chosen to plot the test data series and predictions here as well in order
to be able to gauge the generalisation performance of the model.
We conduct an analysis of the autocorrelation functions in @fig:3_5_acf.

#figure(
  image("output/3_5_acf.png"),
  caption: [Autocorrelation of the residuals of the linear regression model.]
) <fig:3_5_acf>

In @fig:3_5_acf we observe there remain significant correlations in the residuals,
which suggests that the model is not fully capturing the underlying structure of the data.

To further investigate this, we employ a cross-correlation analysis of the residuals and the exogenous variables in @fig:3_5_ccf.

#figure(
  image("output/3_5_ccf.png"),
  caption: [Cross-correlation of the residuals of the linear regression model with the exogenous variables.]
) <fig:3_5_ccf>

Here we observe that some significant correlations remain between the errors and the exogenous variables, which suggests that the model may not be fully capturing the underlying structure of the data and a more complex model may be required.

In order to more formally substantiate this, we employ the Ljung-Box test on the residuals of the model, which reveals that significant correlations remain in the residuals of the model for all choices of lags $k$.

The fact that the autocorrelation in @fig:3_5_acf shows significant correlations with past values
strongly suggests that a model with an autoregressive component may be more appropriate.
Generally such models of linear time-invarant systems (LTI) are well-described by
the use of _transfer functions_, which are a powerful tool for modelling the relationship between the input and output of a system. In the case of an ARX model, its transfer function
contains the appropriate AR components expressed as backshifts of the input series.

In order to build up a model which better captures the underlying structure of the data,
we may use the idea of _transfer function building_, in which we build up a model from the residuals of a simpler model by employing _pre-whitening_, in which models are used as _filters_ and spikes in the CCF plots are used to identify a candidate lag structure.

Pre-whitening would involve fitting an ARMA model on each of the input variables $T_t$ and $G_t$.

$
  phi.alt (B) G_t &= theta (B) epsilon_t
  phi.alt_2 (B) T_t &= theta_2 (B) epsilon_t
$

Which when re-arranged will be transformed into the white-noise component:

$
  phi.alt (B) (theta (B))^(-1) G_t &=  epsilon_t
  phi.alt_2 (B) (theta_2 (B))^(-1) T_t &=  epsilon_(2,t)
$

Which we can then apply as a 'filter' to the input series $P_t$:

$
  phi.alt (B) (theta (B))^(-1) P_t = gamma_t
  phi.alt_2 (B) (theta_2 (B))^(-1) P_t = gamma_(2,t)
$

And then calculate the CCF of $epsilon_(2,t)$ and $gamma_(2,t)$ to get the coefficients of a transfer function from $T_t$ to $P_t$ and analogously for $G_t$.

This is beyond the scope of this assignment, where we will instead follow the model
building outlined in the assignment description.


== AR(1)-X(1,1) Model <sec:3_6>

This time we employ a first order ARX model:
$
  P_(h,t) = -ϕ_1 P_(h,t-1) + ω_1 T_("delta",t) + ω_2 G_(v,t) + ε_t\
$

Where $ϕ_1$ is the AR coefficient, $ω_1$ and $ω_2$ are the exogenous coefficients and $ε_t$ is a white-noise process as before.

We first perform the forecasting and estimation, as shown in @fig:3_6_forecast.

#figure(
  image("output/3_6_forecast.png"),
  caption: [AR(1)-X(1,1) model with $T_"delta"$ and $G_v$ as exogenous variables.]
) <fig:3_6_forecast>

With subsequent error analysis shown in @fig:3_6_errors.

#figure(
  image("output/3_6_errors.png"),
  caption: [Residual analysis of the AR(1)-X(1,1) model.]
) <fig:3_6_errors>

With ACF and CCF analysis shown in @fig:3_6_acf and @fig:3_6_ccf respectively.

#figure(
  image("output/3_6_acf.png"),
  caption: [ACF of the residuals of the AR(1)-X(1,1) model.]
) <fig:3_6_acf>

#figure(
  image("output/3_6_ccf.png"),
  caption: [CCF of the residuals of the AR(1)-X(1,1) model with the exogenous variables.]
) <fig:3_6_ccf>

Comparing @fig:3_5_acf and @fig:3_6_acf we observe that the autocorrelation of the residuals has been significantly reduced, the same of which is true for the cross-correlations shown in @fig:3_6_ccf. This suggests that this only slightly more complex model describes the underlying structure of the process much better than the linear regression model. A Ljung-Box test suggests that the residuals do not have significant autocorrelation left.

== AR(2)-X(2,2) Model <sec:3_7>

Lastly, we are given a second order ARX model:
$
  P_(h,t) =
    -ϕ_1 P_(h,t-1) - ϕ_2 P_(h,t-2)
    + ω_(1,0) T_("delta",t) + ω_(1,1) T_("delta",t-1)
    + ω_(2,0) G_(v,t) + ω_(2,1) G_(v,t-1)
    + ε_t
$

This yields the following estimations:

#figure(
  image("output/3_7_forecast.png"),
  caption: [AR(2)-X(2,2) model with $T_"delta"$ and $G_v$ as exogenous variables.]
) <fig:3_7_forecast>

Yielding the following error analysis:

#figure(
  image("output/3_7_errors.png"),
  caption: [Residual analysis of the AR(2)-X(2,2) model.]
) <fig:3_7_errors>

With ACF and CCF analysis shown in @fig:3_7_acf and @fig:3_7_ccf respectively:

#figure(
  image("output/3_7_acf.png"),
  caption: [ACF of the residuals of the AR(2)-X(2,2) model.]
) <fig:3_7_acf>

#figure(
  image("output/3_7_ccf.png"),
  caption: [CCF of the residuals of the AR(2)-X(2,2) model with the exogenous variables.]
) <fig:3_7_ccf>

We cannot directly observe an improvement when only looking at the ACF and CCF plots, and
it may be difficult to understand whether the increased model complexity is justified.
Again a Ljung-Box test suggests that the residuals do not have significant autocorrelation left.

=== AIC and BIC <sec:3_7_ic>

Now we want to delve a bit deeper into other metrics for model selection, specifically the BIC (Bayesian Information Criterion) and AIC (Akaike Information Criterion). Generically they are calculated as:

$
  upright("AIC") & = - 2 log L (upright(bold(Y)) \| upright(bold(hat(psi)))) + 2 p\
  upright("BIC") & = - 2 log L (upright(bold(Y)) \| upright(bold(hat(psi)))) + p log n\
$<eq:3_7_aic_bic>

with the log-likelihood as:

$ log L (upright(bold(Y)) \| upright(bold(hat(psi)))) = - n / 2 #scale(x: 180%, y: 180%)[\[] log (2 pi) + log #scale(x: 180%, y: 180%)[\(] (upright(bold(Y)) - upright(bold(hat(Y))))^2 / n #scale(x: 180%, y: 180%)[\)] + 1 #scale(x: 180%, y: 180%)[\]] $

for $upright(bold(hat(psi)))$ as the MLE estimated parameters,
$upright(bold(Y))$ as the target/modeled time-series, $n$ as the number
of observations and $p$ as the number of parameters. In our case for the
ARX model, since we only have an AR part, the MLE becomes an OLS
estimation analogous to @eq:3_4_ols_parameter_est, hence $upright(bold(hat(psi))) = hat(Theta)$.

#figure(
  image("output/3_7_aic_bic_model_comparison.png"),
  caption: [AIC, BIC comparison of different models orders from @sec:3_5, @sec:3_6, @sec:3_7]
) <fig:3_7_aic_bic_model_comparison>

This is a classic Elbow curve. In various modelling fields, information theoretical measures or straight residual measures (such as MSE, etc.) are plotted against different model variations. The general idea is, to deduce at which point, there is a good trade-off between parameters of the model and the model performance.
In this case, we have the x-axis as increasing model order, hence increasing complexity and proneness to overfitting. Hence, we are looking for a trade-off between performance and model complexity. It is an interesting choice of metric to do this kind of plot. AIC and BIC inherently punish model complexity, as long as $log(n) > 2$, the BIC does even more so; this also explains the slight difference between those two curves. Thus they theoretically already account for one of the decision parameter that such Elbow-curve is intended to help with.
If we look at the plot, there is a clear trend towards the more complex AR(2)-X(2,2) model. Traditionally, the decision here would be to select the AR(1)-X(1,1) model, as it already performs relatively well and the performance gain to the more complex model is very small. Hence, favour simplicity before performance, in the hope of better generalisation.

However, as the AIC and BIC both already punish the increased model complexity of AR(2)-X(2,2) model, but still produce a lower score, the interpretation would yield in: "more complex, but worth it". Therefore, the model selection would fall onto the AR(2)-X(2,2) model.


== RMSE metric and Prediction on test-dataset<sec:3_8>

We now use the on the training dataset estimated parameters $hat(Theta)$ and use those for one-step predictions on the test dataset with an accordingly constructed design matrix.

#figure(
  image("output/3_8_rmse_model_comparison.png"),
  caption: [RMSE comparison of different models orders from @sec:3_5, @sec:3_6, @sec:3_7 on concatenation of train and test dataset]
) <fig:3_8_rmse_model_comparison_test>

Reflecting upon the question: Does @fig:3_8_rmse_model_comparison_test yield the same model selection as via the AIC, BIC metrics / criteria?

Yes. We already observed on the AIC and BIC on the training dataset that higher model order improves the predictive performance. Even on the distribution of residuals, we could observe increasing statistical evidence for independence of those (via Box-Ljung-Test).
Now we have validated via RMSE on the test dataset that our model generalises well. The AR(2)-X(2,2) model also performs the best in this setting, showing that we have not overfitted with increasing complexity.


== $k$-Step Predictions<sec:3_9>

We selected the AR(2)-X(2,2) model, mainly because of the AIC, BIC selection criteria. From experience, they produce a very reliable selection process.
We provide a selection of $k$ step-widths, which are usually choices that would make sense in real-life settings: 12h, 24h, 48h.

#figure(
  image("output/3_9_kstep_predictions_selection_ols.png"),
  caption: [k-step predictions for intervals]
) <fig:3_9_kstep_predictions_selection>

Beyond that we look at a plot of the RMSE against the step-width $k$ with the intend to draw some conclusions about the error behaviour over longer periods of prediction.

#figure(
  image("output/3_9_rmse_kstep_pred_wo_b.png"),
  caption: [k-step predictions for intervals]
) <fig:3_9_rmse_kstep_pred>

As expected, the RMSE reacts heavily to the amount of steps into the future. Naturally, in the beginning there is a huge error in the predictions, as for k-step predictions, we follow a recursive approach, making predictions upon predicted values. For low $k$ this means, there are not many original values of the target series included, resulting in predictions quickly drifting off.

$P_h$ clearly shows seasonal drops in power (@fig:3_1_analysis). In some intervals however, there are exceptions, where the drop is not as sudden or not as steep. We assume this results in the small buckles in @fig:3_9_rmse_kstep_pred. It is a step width, that goes in phase with the seasonal drops in $P_h$.

In context the most reasonable k-step intervals, as in @fig:3_9_kstep_predictions_selection, yield reasonable RMSEs until $k≤40$, however, not quite as good as the one-step fits. Yet, it does provide the advantage of being able to forecast k-steps into the future.

One could try different starting points for the k-step predictions, instead of just going through the entire dataset. With each further prediction the uncertainty rises, since the newest predictions depend on the older predictions.

Overall, the model performs well; the residuals approach a white noise like distribution. There is still some visible auto-regressive behaviour present, but statistically the evidence is almost sufficient for independence of residuals. Trying a higher order could resolve that.

Apart from the prediction accuracy (which could potentially be improved), in an operational setting the model would work just fine with multi-step predictions. It is not too computationally expensive and would be easily integrated into a prediction pipeline.
Naturally, the further into the future we predict, the higher the uncertainty; that notion is present without calculating the prediction intervals. Beyond that, we saw certain points of the series (towards the drops of $P_h$) where the predictions tend to deviate (@fig:3_9_rmse_kstep_pred), the seasonality of the RMSE over $k$ step-width.
This does not improve if we simply choose a $k$ step-width, that produces a good RMSE, since we would work on a continually predicting on estimated values. At some point the system will take a state (the power drops), where the predictions have shown to be weaker. Hence, by the nature of the series, this point of drop in $P_h$ will naturally come, regardless of the k-step, so this cannot be avoided.

Modelling-wise an RLS model with higher "forgetting" coefficient may be a viable option, to account for these know effects.

Given that knowledge on the historic data, it would be wiser to simply accept an increased uncertainty in specific periods of the series and contextualize consequences: Reduce heating at the drops, but also have power available in case of slight mis-predictions.

It could also be an option to reduce the measurement interval from hourly to 5-min intervals, to be able to react quicker and counter prediction uncertainty.

== Conclusions <sec:3_10_conclusions>

*Some Reflections on the model selection process (@sec:3_5 to @sec:3_7):*

The construction of the selection process was a step by step increase in model complexity.

The AIC, BIC are calculated on the training dataset, because we needed the log-likelihood. We basically test if AIC and BIC yield generalisable conclusions, hence, if our decision on model complexity paid of for the test dataset.

Our addition of the Box-Ljung-Test, fostered the results and contributed to the decision of a more complex model.
It would be interesting to see, when the tipping point arrives, where additional model complexity is not worth it anymore. In the Elbow plot @fig:3_7_aic_bic_model_comparison we could already see a decreasing gradient.

#bibliography("report.bib")
