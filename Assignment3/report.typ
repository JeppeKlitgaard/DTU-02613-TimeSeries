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
where power generation is at its peak,
due to the overly large effect of the variance during the winter months on the residual error estimates.

The implementation of the prediction intervals can be found in the attached Jupyter Notebook `JK_2.ipynb`.

= An ARX Model for the Heating of a Box <sec:3>

We are given a data set of hourly measurements from a box instrumented with a heater,
an internal temperature sensor and an external temperature sensor.

The box has window on its south-facing wall, onto which the vertical solar radiation is also recorded in the data set.

The internal temperature of the box was kept approximately constant using a thermostatic control of the internal heater.


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

TODO ELABORATE ON WHAT THIS DOES NOT SHOW?

== ARX Models

In this section, we will be working a different flavours of ARX models (extended auto-regressive models), thus it is key, to define a common notation especially for the eXogenous variables.

Therefore, our clear and generic notation for ARX models with
multiple exogenous variables:

$ upright("AR") (p) - upright("X") (e_1 \, e_2 \, . . .) $

where $p$ is as in most literature, the order of the AR part of the
model, $e_1$ is the order of the first exogenous variable, $e_2$ the
order of the second exogenous variable and so on. Formally, we have:

$ 
Y_t = c - sum_(i = 1)^p phi.alt_i Y_(t - i) + sum_(i = 0)^(e_1) omega_(1 \, i + 1) X_(1 \, t - i) + sum_(i = 0)^(e_2) omega_(2 \, i + 1) X_(2 \, t - i) + . . . + epsilon_t 
$<eq:3_4_generic_arx>

To further unify notation, we assume a model, where ${ T_t }$ represents the series $T_(upright("delta"))$, ${ G_t }$ is the series $G_v$ and ${P_t }$ is the series $P_h$ for $t in { 1 \, . . . \, 167 }$ in the training dataset.


== Impulse Response

There are different options for modelling an impulse response for a model with 2 exogenous variables: $P_t$ with $G_t$ or $G_t$ as exogenous variable, as well as a model for $P_t$ with both as exogenous variables.

$
upright("AR(p)-X(e1, e2)") = c + sum_(i = 0)^p phi.alt_i B^i P_t + sum_(i = 0)^(e_1) omega_i B^i T_t + sum_(i = 0)^(e_2) beta_i B^i G_t + epsilon_t 
$<eq:3_4_context_arx>

for both as exogenous, given parameters/coefficients $phi.alt_i \, beta_i \, omega_i in bb(R)$.

Now the instructions are not entirely clear, because an impulse response needs a system for which it is calculated. It is neither indicated what the ARX model order $p$ should be, neither the order of the exogenous series (i.e. which values do we choose for $e_1, e_2$).
Due to ambiguation throughout the course and this assignment, it is also unclear what is referred to by "up to lag 10". This could refer to the lagged samples of the modeled series $P_t$, hence the AR components, or it could refer to the depth of recursion for the impulse response function.
Further it is unknown from which of the exogenous variables, $G_t$ or $T_t$ the unit impulse should come from.

Therefore, we need to make reasonable choices for the model. We choose an AR(10)-X(1,1) model with a given lag of of $k=20$ for the impulse response function IR(k). On a training dataset of 167 samples it seems reasonable to plot the impulse response on a sizable lag, to actually see some impact and behaviour.
An AR order of $p=10$ also seems reasonable to actually have some coefficients that capture the model nuances, without immediately decaying the response.
The unit impulse shall be given from the variabele $T_t$ as @fig:3_3_pairplot shows, it has a positive correlation with the target series $P_t$, while @fig:3_1_analysis clearly indicates, the variable $G_t$ has an inverse relationship with the target $P_t$. This does carry quite valuable information, however, the impulse response will decay extremely fast, when having opposing effects.

Thus the final AR(10)-X(1,1) model becomes:

$ 
P_t = c - sum_(i = 1)^10 phi.alt_i P_(t-i) + omega_1 T_t + beta_1 G_t + epsilon_t 
$<eq:3_4_ar10_x1_1>

For the impulse response function, we give the unit impulse from $T_t$ and use the recursion

$ 
upright("IR") (0) & = 0 phi.alt_1+ 1 omega_1 + 0 beta_1\
upright("IR") (1) & = upright("IR") (0) phi.alt_1 + upright("IR") (0) omega_1 + upright("IR") (0) beta_1 \
upright("IR") (2) & = upright("IR") (1) phi.alt_1 + upright("IR") (0) phi.alt_2 + 0 omega_1 + 0 beta_1 \
upright("IR") (3) & = upright("IR") (2) phi.alt_1 + upright("IR") (1) phi.alt_2 + upright("IR") (0) phi.alt_3 + 0 omega_1 + 0 beta_1\
 & #h(0em) #h(0em) dots.v\
upright("IR") (k = 15) & = upright("IR") (14) phi.alt_1 + upright("IR") (13) phi.alt_2 + ... + upright("IR") (5) phi.alt_10 + 0 omega_1 + 0 beta_1\
$<eq:3_4_ir_recursion>

In sum notation this would be with $upright("IRF")(0) = 1 omega_1 + 0 beta_1$ set as initial value:

$
  upright("IRF")(k) = sum_(i=1)^k sum_(j=0)^(min(i,p)) [ phi.alt_j upright("IRF") (i-j-1) bb(I)_(i-j-1≥0) ]
$<eq:3_4_ir_summation>

for $k in bb(N)$ the lag and $p$ the AR model order.

Having such a relatively complex model, there is no way to "estimate" it with sufficient accuracy. One could potentially read-off maybe the first two $phi.alt_1, phi.alt_2$ parameters of the AR part via an ACF and PACF plot to conclude that the $phi.alt_i$ coefficients are all $|phi.alt_i|<1$. Maybe even find an argument for an alternating sign of he coefficients, based on the PACF plot. However, none of this will yield a usable result. Thus, we need to actually fit the model to the data. 

We will use the the ```AutoReg```class from ```statsmodels.tsa.ar_model``` in Python, which allows to provide an exogenous variable and has built in parameter estimation (via OLS and conditional MLE).
To double check, we also did our own classic OLS fit on the design matrix:

$
  X = [1, P_(t-1), P_(t-2), P_(t-3), ..., P_(t-10), T_t, G_t]
$

noted in the form of column vectors of the corresponding series, which yields a parameter vector of $Theta = [c, - phi.alt_1 ..., - phi.alt_(10), omega_1, beta_1]^T in bb(R)^(p+2)$.

$  
  arrow.r.double &  & P_(t - k) & = X dot.op Theta + epsilon_t = Y \
  arrow.l.r.double &  & hat(Theta) & = (X^T X)^(- 1) X^T Y \
$<eq:3_4_ols_parameter_est_ar10_x1_1>

The resulting parameters are given as:

#figure(
  table(
    columns: 13,
    table.header($c$, $phi.alt_1$, $phi.alt_2$, $phi.alt_3$, $phi.alt_4$, $phi.alt_5$, $phi.alt_6$, $phi.alt_7$, $phi.alt_8$, $phi.alt_9$, $phi.alt_(10)$, $omega_1$, $beta_1$),
    [3.3073], [0.3294], [0.0069], [-0.0176], [0.086], [0.0124], [0.0308], [0.0171], [-0.0191], [0.0081], [-0.0022], [2.0537], [-0.0925], 
  ),
  caption: [AR(10)-X(1,1) coefficients with $T_d$ and $G_v$ exogenous of order 1],
) <table:3_4_irf_arx_coeff>

Now we can calculate the impulse response via the recursive or the short-sum notation above as in @eq:3_4_ir_summation. 

#figure(
  image("output/3_4_impulse_response_Tt.png"),
  caption: [Impulse Response of the AR(10)-X(1,1) with the unit impulse from $T_d$]
) <fig:3_4_ir_Tt>

If the phrasing "Present it for both variables" is interpreted as presenting the impulse response with a unit impulse from both variables separately, we can give a second plot of just an impulse from $G_v$:

#figure(
  image("output/3_4_impulse_response_Gv.png"),
  caption: [Impulse Response of the AR(10)-X(1,1) with the unit impulse from $G_v$]
) <fig:3_3_ir_Gv>

As already explained based on @fig:3_1_analysis, the two impulse responses counter-act. The IR (impulse response) for this model decays rather quickly, especially with the impact of an exogenous variable that is reciprocal to the target variable. Most of these effects were already assumed based on available information.

We can interpret, that the model with order $10$ is rather robust to unit shocks from both exogenous variables. As a consequence, with this model, we can accept a certain degree of fluctuation in the exogenous variables, without having to worry about a big impact on the predictive performance of our model. For example, the effect of a measurements error or sensor defect in $T_d$ or $G_v$ that acts as a shock to our system, will not last for very long.
Thus, having a larger order for the AR part of our model pretects against shocks from the other variables. On the other hand, the model is more vulnerable to sudden changes in auto-regressive behaviour.

== Model Selection (3.5 - 3.8)

In this section, we will stick to the notation convention for ARX models introduced above and create a selection of models for comparison.

=== Linear Regression

For the simple linear regression (OLS) model we have:

$
  P_t = omega_1 T_t + beta_1 G_t + epsilon_t
$<eq:3_5_ols_model>

with the design matrix: $X = [T_t \, G_t]$ for the column-vectors
and with $Theta = [omega_1 \, beta_1]^T$ as parameter vector. The estimation (fitting the model) for $hat(Theta)$ will be analogous to @eq:3_4_ols_parameter_est_ar10_x1_1.

The resulting parameters are $hat(Theta) = [3.8948 \, -0.1099]^T$ with a corresponding RMSE of $approx 69.87$. Thus we have our linear regression predictions as $X dot.op hat(Theta) = hat(P_t)$.
Further we can shall do one-step predictions (as in @sec:3_8_OSPred_RMSE), which is a way of iteratively re-fitting the model. It is predicting the values of $hat(P_t)$ one step ahead, with the parameters $hat(Theta_(t-1))$ fitted onto the previous time-step, hence, also the design matrix $X_(t-1) = [T_(0, t-1), G_(0,t-1)]$ constructed until the previous time-step. Consequently, for $k=1$:

$
  hat(P)_(t+k|t) = X_(t+k) dot.op hat(Theta_(t))
$<eq:3_5_os_pred>

with

$
  X_(t+k) &= [T_(0, t+k), G_(0,t+k)] quad upright("series from observations 0 to ") t+k
  hat(Theta_t) &= (X_t^T X_t)^(-1) X_t^T P_t
  hat(epsilon)_(t) &= P_(t+k) - hat(P)_(t+k|t)
$

PLOTTING

- observations of PLOTTING
- need for a transfer function?

=== One-Step Predictions & RMSE (3.8)<sec:3_8_OSPred_RMSE>

When producing one-step predictions, a critical factor is to mind the burn-in period. As the model parameters are re-fitted at each step of the prediction, taking into account the new predicted values of the series, in the first few steps, there is little data to fit on. Thus, the predictions are not very strong. This is called the burn-in period.

As we must give one-step predictions on the test-dataset (which is already quite few observations), it makes sense to include some of the last observations of the train-dataset, to counter the burn-in.
This way, there is hopefully not much of an impact on the predictions of the test-dataset. 
Additionally, this would also be a more fair comparison to the previous exercises, as those are OLS fitted on the entire train-dataset, thus do not suffer from a burn-in, as do one-step predictions.

To test out, how much of a spill from the train-dataset we need, we can check fitting the one-step prediction on *only* the test-data. 

INSERT FIGURE

The decision rule on how much burn-in we grant the model is: Increase the burn-in $b$ (meaning, exclude the first $b$ observations and predictions) until there are no more negative values in the one-step predictions of $P_h$. This is a context based rule, as negative wattage for the heater, the variable $P_h$ represents, does not make sense.
The FIGURE shows the results for $b=11$. Therefore, we include the last 11 samples/observations of the train-dataset into the design matrix $X$ for the one-step predictions.

Now including 11 samples from the train-data as spill over, we observe the following:

On the new concatenated series, we really only need $b ≥ 6$ to not have absurd values in $P_h$ (negative or >1000), BUT with $b=8$, we now effectively have 64 one-step predictions as asked for in the task.

We added 11 samples from the training data to the 64 samples of the test-data. So we fitted step-wise on a total of 64+11=75 observations. Now we allow a burn-in of 8 and remain with 64 one-step samples... so $64+11-8=64$?!

That begs the questions: Where have the 3 samples gone lost?

- the highest model order of 2 for the AR(2)-X(2,2) swallows 2 samples
  in padding/cutting of the series
- 1 samples is swallowed by the process of one-step predictions, since
  we always look one sample back for $hat(y_(t \| t - 1))$
- in fact, out of 64 test-data samples, we would only get 63 one-step
  predictions, as we need to start ON the 2nd sample (t=2) to fit the
  model on the 1st sample at:
  $hat(y_(t \| t - 1)) = hat(y_(2 \| 2 - 1)) = hat(y_(2 \| 1))$
- otherwise we would hit $t - 1 = 0$ which is not feasible

- therefore the RMSE formula given in the task is actually incorrect, as it averages over 3 more values than there are one-step predictions available (64 test-data samples; 3 casualties, 2 for model order and 1 for OS-predictions)
- that is not even considering a burn-in!
- hence the effective number of observations fair for comparison to the previous models is much lower than 64
- please, fucking think about what you write Peder, at least 64 fucking times, before you confuse the shit out of people for hours!

Reflecting upon the question: Does this yield the same model selection as via the AIC, BIC metrics / criteria?

Yes and no. Only going via the RMSE, we could conclude from the first approach, where we only used pure test-dataset samples, that we should select the AR(2)-X(2,2) model.
However, after adjusting the setting (to have a spill of samples from the train-dataset) for a technically more fair comparison, the RMSE would yield the conclusion to select the AR(1)-X(1,1).
However, to put this into perspective, the RMSE is overall much lower with the second modelling approach. This makes sense makes sense, since there is much more information available to the model. This lets us conclude that there is a lot of relevant information in the further past. 
Yet, there seems to be a notion that a model with shorter lag, an AR(1) component, captures this better than a longer lag. This could be interpreted as that past information is overall valuable, but the series behaviour is ultimately very short-term oriented.

TODO

#bibliography("report.bib")
