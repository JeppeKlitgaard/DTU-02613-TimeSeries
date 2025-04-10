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


#set heading(numbering: "1.1")
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

== Realisations

=== 1.1 & 1.2
We simulate the given process 5 times using the ```SARIMAX``` module with $n=200$ observations and the coefficients set to: $ϕ_1 = -0.6$ and $ϕ_2 = 0.5$

#figure(
  image("img/blabla.png"),
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
  (1 - ϕ_1B)(1 + Φ_1B^{12})(log(Y_t)-µ) = ε_t
$ <eq:2_seasonal_ar>

where $Y_t$ is the monthly energy from the plant in MWh, ${ε_t}$ is a white-noise process with variance $σ_ε^2$ . The parameters $ϕ_1 = −0.38$, $Φ_1 = −0.94$ and $µ = 5.72$ are assumed to be known. Based on 36 observations, it is found that $σ_ε^2 = 0.22^2$.


= OLD PART!!!!!!!!!!!!!!!!!!!!

In order to simulate the seasonal processes we utilise the Python library `statsmodels`,
which offers an implementation of the
`SARIMAX` (Seasonal AutoRegressive Integrated Moving Average with eXogenous regressors) model.

Looking at the definition of the SARIMAX implementation and comparing against @eq:2_seasonal_arima, we confirm that the parameters are defined in a similar fashion,
meaning there should not be any sign transformations needed for the parameters.

We employ initial conditions $vv(y_t) = vv(0)$ and avoid burn-in effects by
simulating the models for $N_"burn-in" = 10000$ before simulating the process.

For all simulations, we use $n=1000$ observations and let $ε_t ∼ 𝒩(0, 1)$.

When computing the (partial) autocorrelation functions we use `plot_(p)acf`
from `statsmodels` with $N_"lags"=30$ and a significance level of $p=0.05$ for the confidence intervals.

In order to understand the behaviour of the simulations,
we refer to Table 6.1 in @Madsen_2008[p.~155], which we reproduce here:

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

==

#figure(
  image("output/plot_2_1.png"),
  caption: [Simulated AR(1) process with $ϕ_1 = 0.6$.],
) <fig:2_1>

In @fig:2_1 we note the exponential decay of the autocorrelation function and single significant value in the partial autocorrelation function.
We also find that the exponential decay of the autocorrelation function matches the
expected shape from @eq:1_autocorrelation_recursion.

==

#figure(
  image("output/plot_2_2.png"),
  caption: [Simulated Seasonal AR(1) process with seasonality $s=12$ and $Φ_1 = -0.9$.]
) <fig:2_2>

In @fig:2_2 we note that the autocorrelation function is periodic with period $s=12$ while still decaying exponentially. Because of the negative sign of $Φ_1$ we find that the value of $ρ(k)$ cycles around $0$.

Additionally, we find that the significant point in the partial autocorrelation function has shifted to $k=12$. Additional points sticking out of the confidence interval are assumed to be statistical noise as would be expected with a confidence level $p=0.05$.

==

#figure(
  image("output/plot_2_3.png"),
  caption: [Simulated Seasonal ARMA process with seasonality $s=12$ and parameters $ϕ_1 = 0.9$ and $Θ_1 = -0.7$.]
) <fig:2_3>

Referring to @fig:2_3, we note that the AR part of the process contributes to a single significant value in the partial autocorrelation function at $k=1$ that is then convolved with a decaying 'Dirac brush' with periodicity $s=12$.

The autocorrelation function can be observed to decay away as expected for the AR(1) part of the model,
but interestingly also exhibits some periodicity with period $s=12$ arising from the MA(1) seasonal component.

==

#figure(
  image("output/plot_2_4.png"),
  caption: [Simulated Seasonal AR process with seasonality $s=12$ and parameters $ϕ_1 = -0.6$ and $Φ_1 = -0.8$.]
) <fig:2_4>

The picture emerging from @fig:2_4 gets more muddied for this more complex model, though using @table:order_identification we would expect the autocorrelation function to have a exponentially decaying sine envelope that is repeated every $s=12$ observations, also decaying for each repetition.

For the partial autocorrelation we would expect to see a single significant value at $k=1$ and then additional significant values around $k=12$. Due to the interaction between the regular and seasonal AR(1) processes, we would expect there to be several significant values, though the exact number is difficult to deduce using intuition alone.

==

#figure(
  image("output/plot_2_5.png"),
  caption: [Simulated Seasonal MA model with seasonality $s=12$ and parameters $θ_1=0.4, Θ_1=-0.8$.]
) <fig:2_5>

From @table:order_identification we find that we would expect an exponential decay of the partial autocorrelation function for an MA(q) process with a limited number of significant values in the autocorrelation function, which is also reflected by @fig:2_5.
With moderate difficulty we observe a single significant value in the initial part of the autocorrelation function as we would expect from the regular MA(1) process. This is then repeat once around $k=12$ by the seasonal MA(1) process. We note that both the $k=0$ and $k=1$ values from the regular MA(1) process are convolved with the single expected peak from the seasonal MA(1) process to produce 3 significant values at $k∈{11, 12, 13}$.

The parital autocorrelation function is more difficult to interpret, but here we would expect
a decaying exponential envelope, potentially over a sine function.
While the the _inner_ exponential function understood to arise from the regular MA(1) process with $θ_1=0.4$ decays very quickly as observed around $k=1$, we observe the expected periodicity of $s=12$ effectuated by the seasonal MA(1) component.

==

#figure(
  image("output/plot_2_6.png"),
  caption: [Simulated Seasonal ARMA model with seasonality $s=12$ and parameters $θ_1 = -0.4, Φ_1=0.7$.]
) <fig:2_6>

Lastly, in @fig:2_6 we observe an exponential decay in the partial autocorrelation function which we attribute to the regular MA(1) progress. We also observe it repeated at $k=12$, though not at $k=24$, as would be expected by the AR(1) seasonal component.

In the autocorrelation function we find a decaying exponential envelope over a 3-element sequence of alternating signs repeated at $k=12$ and $k=24$, which again is
consistent with the rules outlined in @table:order_identification.

== Summary of Identifications

We have perhaps cheated a bit by using the rules from @table:order_identification to identify the processes, as the conclusions outlined in @table:order_identification should have been
deduced and presented here.

However, we can add additional commentary - for instance, we realise that more complex processes would be very difficult to identify simply by inspection of the autocorrelation and partial autocorrelation functions. Instead, we propose that a model is iteratively built in order to identify the parameters of such processes. Here a single AR(1) or MA(1) model may be extended appropriately by first fitting the model to a realisation of the process and then inspecting the residuals for any remaining autocorrelation. A new simple model may be fit on the residuals and combined with the original model to construct a better fitting model, which may then again be improved iteratively by inspection of the residuals.

The seasonality parameter will generally be relatively easy to deduce from autocorrelations, especially if the process only features a single periodicity with little variance in the periods.

Overall it is probably best to use parameter estimation techniques. In a
real world scenario, one would likely have the time-series as actual
data and not just as a plotted graph. There are two main ways to
estimate especially $phi.alt$ and $theta$ for ARMA models:

+ OLS regression for AR models
  $upright(bold(hat(phi.alt))) = (X^T X)^(- 1) X^T upright(bold(y_(t - k)))$
  with a feature matrix X constructed from the $t - k$ to the $t - N$ th
  samples, with order of the AR model $k$
+ MLE for MA models to estimate $upright(bold(hat(theta)))$ based on a
  likelihood function that assumes Gaussian $upright(bold(epsilon))_t$

Since none of the given model display any trends, but rather are all
stationary, a de-trending via differencing on real world data may be
appropriate. Otherwise combining a standard trend model like OLS, RLS or
WLS as a prior model, could be useful.


= Identifying ARMA(p,q) Models <sec:3_arma_identification>

Now, using the rules established and tested in @sec:2_simulating_seasonal_processes,
we are able to infer the order of the 3 processes presented in the assignment description.

Additionally, we use the same function used to generate the realisations in @sec:2_simulating_seasonal_processes to inspect a realisation of the proposed ARMA structure for each of the 3 given processes.

We note that we are informed that the processes are ARMA(p,q) models
and as such not expected to feature any seasonal or differencing components.

== Process 1

In the given plots of the first process, we observe no significant values in the autocorrelation nor partial autocorrelation function. As such, we conclude that the process is an ARMA(0,0) process, which is equivalent to a white noise process. From the timeseries plot, we find that the variance is likely unitary. Recreating a realisation of such a process, we obtain @fig:3_1.

#figure(
  image("output/plot_3_1.png"),
  caption: [Simulated ARMA(0,0) process with $ε_t ∼ 𝒩(0, 1)$]
) <fig:3_1>

== Process 2

For the second process we find two significant values in the partial autocorrelation function
and a double-exponential decay in the autocorrelation function, which we immediately identify as an ARMA(2,0) process in accordance with @table:order_identification.

A realisation of such a process is shown in @fig:3_2.

#figure(
  image("output/plot_3_2.png"),
  caption: [Simulated ARMA(2,0) process with $ϕ_1 = 0.45, ϕ_2 = 0.3$]
) <fig:3_2>

== Process 3

The third process is rather more difficult to identify.
Squinting slightly, we find the autocorrelation function to appear as if it is decaying
exponentially with two distinct characteristic decay constants, which would suggest an AR order of $p=2$. Finding the envelope of the absolute partial autocorrelation function to resemble a decaying exponential, we find that the MA order may be $q=1$.

The cyclical nature of the partial autocorrelation function suggests that one of the autoregressive parameters must be negative.

Fiddling a bit with the magnitude of the parameters, we obtain @fig:3_3, which matches the the plot of the third process in the assignment description well.

#figure(
  image("output/plot_3_3.png"),
  caption: [Simulated ARMA(2,1) process with $ϕ_1 = 1.1, ϕ_2 = -0.2, θ_1 = 0.9$]
) <fig:3_3>


=== Commentary:

In addition, it would be a good idea to start some residual analysis on the given models. As mentionend in Ex 2.7: in a real world scenario, we would have a time-series as actual data.

The main difficulty is the estimate the order of models, so the parameters $(p,d,q)$. From there, building a model can be started done by fitting and identifying parameters numerically.
Ideally, we can start with models of lower order, thus lower complexity and then analyse the residuals between model and data. As metrics, we suggest AIC and BIC, as those also penalize complexity of a model, which prevents overfitting. This will be especially useful when combining with regression methods (OLS, WLS, RLS).

Furthermore, the given chart indicates some notion of auto-regressive seasonality. Yet, the period of seasonality would be too long to actually fit a model of reasonable complexity (as the ACF and PACF do not allow estimations for such long periods; we are talking about 100 periods towards the end of the series).
Therefore, we stick to the simpler ARMA(2,0,1) model. As it does fit the model well.

#bibliography("report.bib")
