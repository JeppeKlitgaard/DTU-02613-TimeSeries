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
      02417 Times Series Analysis - Assignment 2
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
  *02417 Times Series Analysis - Assignment 2*
])

Authors:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-03-24

#pagebreak()

= Stability <sec:1_stability>

We are given an AR(2) process ${Y_t}$ with residual term $ε_t$ arising from a stochastic process ${ε_t} ∈ 𝒩(0, σ_ε^2=1)$ :
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

We are given a seasonal ARIMA model:
$
  ϕ(B) Φ(B^s) ∇^d ∇^D_s y_t = θ(B) Θ(B^s) ε_t
$ <eq:2_seasonal_arima>

To grasp a better understanding, it is helpful to explain the role of some variables in the model:

$
s & arrow.r upright("seasonal shift") &  & upright("periodicity of seasonal ARIMA")\
p & arrow.r upright("lag of time-series") &  & upright("for AR, polynomial order")\
 &  &  & phi.alt (B) Y_t = Y_t (1 + phi.alt_1 B + phi.alt_2 B^2 . . . phi.alt_p B^p) = epsilon_t\
q & arrow.r upright("lag of random noise") &  & upright("for MA, polynomial order")\
 &  &  & theta (B) epsilon_t = epsilon_t (1 + theta_1 B + theta_2 B^2 . . . theta_p B^p) = Y_t\
phi.alt & arrow.r upright("AR") (p) &  & upright("coeff. for the auto-regressive part")\
Phi & arrow.r upright("AR") (P) &  & upright("coeff. for the seasonal auto-regressive part")\
theta & arrow.r upright("MA") (q) &  & upright("coeff. for the moving-average part")\
Theta & arrow.r upright("MA") (Q) &  & upright("coeff. for the seasonal moving-average part")\
 &  &  & upright("hence shift ") B^s\
nabla^d & arrow.r nabla^d Y_t = Y_t - Y_(t + d) &  & upright("difference shift of normal model")\
nabla_s^D & arrow.r nabla_s^D Y_t = Y_(t + s) - Y_(t + s + d) &  & upright("difference shift of seasonal model")\
$


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
  caption: [Simulated ARMA(2,1) process with $ϕ_1 = 0.45, ϕ_2 = 0.3, θ_1 = -0.5$]
) <fig:3_3>

#bibliography("report.bib")
