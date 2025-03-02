#import "@preview/physica:0.9.4": super-T-as-transpose
#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *

#show: codly-init.with()
#codly(languages: codly-languages)


#set heading(numbering: "A.1.1")
#set math.equation(numbering: "(1)", supplement: [Eq.])
#show: super-T-as-transpose
#set table.header()

#let below(x, b) = {
  math.scripts(math.attach(math.limits(x), b: b))
}

#let underbar(x) = {
  below(x, math.macron)
}

#let mm(x) = {
  let out = math.upright(x)
  out = math.bold(out)
  out = underbar(out)
  out = underbar(out)

  out
}

#let vv(x) = {
  let out = math.upright(x)
  out = math.bold(out)
  out = underbar(out)

  out
}

#let mref(key) = {
  let lbl = label(key)
  math.text([(#ref(lbl, supplement: none))])
}

// Math defines
#let θ = $vv(θ)$
#let θ_pred = $hat(#θ)$
#let X = $mm(X)$
#let X_pred = $hat(#X)$
#let y = $vv(y)$
#let y_pred = $hat(#y)$
#let ε = $vv(ϵ)$
#let W = $mm(W)$

= Assignment 1

== Plot Data

We are given a tabular data set "`BIL54`" describing the total number of registered motor driven vehicles in Denmark. The data given as a time series with monthly observations starting in January 2018.

We divide the data set into a _training set_ containing data from January 2018 to December 2023, and a _test set_ containing data from January 2024 to December 2024.

=== Construct time variable

The data set is indexed by an `ISO8601` date-time column, which we transform into a floating point time variable, $x$ as follows:

$
  x = "Year" + ("Month")/12 wide wide "Month" &∈ {0, 1, …, 11}\
                                       "Year" &∈ {2018, 2019, …, 2024}\
$

When referring to a vector of time values, we will follow the vector notation $vv(x)$. Other vectors are denoted similarly, while matrices are denoted with a double underline, $mm(X)$.

It is important to distinguish between the time-like vector-valued variable $vv(x)$ and the design matrix $mm(X)$, which will be introduced later.

=== Plotting observations

Using only the training data, we plot the total number of registered vehicles in Denmark as a function of time in @fig:plot_observations.

#figure(
  image("output/plot_observations.png"),
  caption: [Training data describing the total number of registered vehicles in Denmark as a function of time.]
) <fig:plot_observations>

Apart from a significant jump between 2021 and 2022, the data follows roughly a linear trendline, especially between 2018 and the beginning of 2021.
After the mentioned jump, the linear trend continues for approximately 6 months into 2021.
Then the dataset shows signs of stagnating growth for the total number of cars.
From end 2021 to 2024, the data roughly follows linearly trend again, though notably with a smaller slope than previously.
This can be interpreted as a saturation in car registration across Denmark.

We pose that reasons for this could be:
1. Most of the population already has a car and no further need to register a new car
2. Due to some legislation or other factor it has become less popular of viable to register a new car, which causes weaker demand for new car registrations

== Linear Trend Model <sec:2_linear_trend_model>

We are given the General Linear Model (GLM):
$
  Y_t = θ_i + θ_2 ⋅ x_t + ϵ_t
$ <GLM>

=== Matrix Form

We immediately rewrite @GLM as a matrix, observing the tensor notation outlined previously:

$
  #y = #X #θ + #ε
$ <GLM_matrix>

Where $#X$ denotes the _design matrix_, $#θ$ the _parameter vector_, and $#ε$ represents a stochastic noise term, $#ε ∼ 𝒩(0, σ^2)$.

In our case, the _design matrix_ is constructed with an intercept $θ_1$ and a single trend parameter $θ_2$:
$
  #X = mat(delim: "[",
    1, x_1;
    1, x_2;
    ⋮, ⋮;
    1, x_n
  )
$

Using the first 3 data points of the given data, we can construct the following
model:

$
  #y &= #X #θ + #ε \
  mat(delim: "[",
    y_1;
    y_2;
    y_3
  ) &= mat(delim: "[",
    1, x_1;
    1, x_2;
    1, x_3;
  ) mat(delim: "[",
    θ_1;
    θ_2;
  ) + mat(delim: "[",
    ϵ_1;
    ϵ_2;
    ϵ_3;
  ) \
  mat(delim: "[",
    2930483;
    2934044;
    2941422;
  ) &= mat(delim: "[",
    1, 2018.000;
    1, 2018.083;
    1, 2018.167;
  ) mat(delim: "[",
    θ_1;
    θ_2;
  ) + mat(delim: "[",
    ϵ_1;
    ϵ_2;
    ϵ_3;
  )
$

=== Parameter Estimation <sec:2_2_paramater_estimation>

We can estimate the parameters $#θ$ leveraging the _normal equations_.

This is done by minimizing the _residual sum of squares_ (RSS) between the observations and the model predictions, that is:

$
  #θ_pred = limits("argmin")_θ norm(#y - #X #θ)^2
$ <theta_pred_optimisation>

We state without proof that the solution to the above optimization problem is given by the _normal equations_:

$
  mref("theta_pred_optimisation") quad ⇒ quad
  #X^T #y = #X^T #X #θ_pred quad ⇔ quad
  #θ_pred = (#X^T #X)^(-1) #X^T #y
$ <theta_estimation>

Where $#X$ is assumed to be invertible.

We additionally consider the standard errors of the parameter estimates, $#θ_pred$.

Under the assumption that the observations are described by @GLM, that is, the data is drawn from a _Simple Linear Model_ overlaid with a stochastic (i.i.d) noise term,
we find that the residuals of the mdoel are exactly the noise term, $#ε$.

$
  #ε = #y - #y_pred
$

Where

$
#y_pred = #X #θ_pred
$ <y_pred>


For the benefit of the reader,
we reproduce the relationship between the residuals and the covariance matrix,
$mm(Σ)$ @Madsen_2008[p.~36-37]:

From @theta_estimation, we have:
$
  𝔼[#θ_pred]
  &= 𝔼[(#X^T #X)^(-1) #X^T #y] \
  &= 𝔼[(#X^T #X)^(-1) #X^T (#X #θ + #ε)]    wide wide && mref("GLM_matrix") \
  &= 𝔼[(#X^T #X)^(-1) (#X^T #X) #θ]    wide wide && #ε ∼ vv(𝒩)(0, σ^2) \
  &= #θ
$

We consider the following:
$
  #θ_pred - 𝔼[#θ_pred]
  &= (#X^T #X)^(-1) #X^T #y - #θ \
  &= (#X^T #X)^(-1) #X^T (#X #θ + #ε) - #θ \
  &= (#X^T #X)^(-1) #X^T #ε
$

We can now evaluate the covariance of the predicted parameters, #θ_pred:
$
  "Cov"[#θ_pred]
  &= 𝔼[(#θ_pred - 𝔼[#θ_pred]) (#θ_pred - 𝔼[#θ_pred])^T]\
  &= 𝔼[((#X^T #X)^(-1) #X^T #ε) ((#X^T #X)^(-1) #X^T #ε)^T]\
  &= 𝔼[(#X^T #X)^(-1) #X^T #ε #ε^T #X ((#X^T #X)^(-1))^T]\
  &= (#X^T #X)^(-1) #X^T 𝔼[#ε #ε^T] #X ((#X^T #X)^T)^(-1)\
  &= (#X^T #X)^(-1) #X^T "Var"[#ε #ε^T] #X ((#X^T #X)^T)^(-1)\
  &= (#X^T #X)^(-1) #X^T σ^2 #X (#X^T #X)^(-1)    wide wide wide wide && #ε ∼ vv(𝒩)(0, σ^2) \
  &= σ^2 (#X^T #X)^(-1) \
$ <eq:cov_theta_pred>

Where we understand the variances of the predicted variables to be
the diagonal elements of the covariance matrix:

$
  "Var"[#θ_pred]
  &= "diag"("Cov"[#θ_pred])
$

Which gives the following standard errors for the parameter estimates:

$
  σ_(hat(θ)_i) = sqrt("Var"[hat(θ)_i]) = sqrt(σ^2 (#X^T #X)^(-1)_(i)) wide wide i ∈ {1, 2}
$

Where $σ$ is calculated using the _residual sum of squares_ (RSS) with $N$ and $p$ denoting the number of rows and parameters in the design matrix #X respectively:

$
  σ = norm(#y - #X #θ_pred)^2/sqrt(N - p)
$ <eq:sigma>

Computing these for the entire training dataset gives:

$
  hat(θ)_1 &= (-110 ± 4) × 10^6\
  hat(θ)_2 &= (56.1 ± 1.8) × 10^3\
$

Using these predicted parameters,
we are able to produce an estimate of the vehicle registrations for the training dataset
using the General Linear Model as seen in @fig:plot_observations_predicted.

#figure(
  image("output/plot_observations_predicted.png"),
  caption: [Estimation of vehicle registrations using a Linear Trend Model.]
) <fig:plot_observations_predicted>

=== Prediction <sec:2_3_prediction>

We now wish to predict future vehicle registrations using our simple model.
From @y_pred, we see that doing so simply requires the construction of an appropriate
design matrix, #X_pred.
This may be constructed simply as:

$
  #X_pred _i = mat(delim: "[",
    1, 2024 + (i-1)/12
  ) wide wide i ∈ {1, 2, …, 12}
$

Where $i$ indexes the rows of the design matrix.

Carrying out the forecasting as described by @y_pred
along with appropriate estimation of the confidence interval of our prediction,
we obtain @table:forecast_2_3.
It should be noted that the estimation of the confidence interval is valid only
in the case where the observations are described by the General Linear Model.

#let forecast_2_3 = csv("output/ols_forecast.csv")
#figure(
  table(
    columns: 4,
    table.header(..forecast_2_3.remove(0).map(strong)),
    ..forecast_2_3.flatten()
  ),
  caption: [12 month forecast of vehicle registrations in Denmark.]
) <table:forecast_2_3>

=== Plot of Forecast

We can now present the forecasted data in @table:forecast_2_3
along with the training, test, and prediction data sets

#figure(
  image("output/plot_observations_forecast.png"),
  caption: [12 month forecast of vehicle registrations in Denmark.]
)

=== Commentary on Forecast <sec:2_5_commentary>

We find that the prediction is a relatively poor match against the test data,
from which we understand that a Simple Linear Model is likely not an appropriate
model for the data.

In particular, we observe significant local deviations from the model,
which can be understood by considering the modelling domain.
It is reasonable to assume that the number of registered vehicles will
be influenced by market conditions,
such as government subsidies, registration fees, and taxation schemes.
Additionally, we would expect supply chain disruptions to have strongly influence
the contemporary pricing and availability of vehicles.

A better model would likely incorporate these factors.
A model that incorporates locality without apriori domain knowledge could be
a _Weighted Least Squares_ model with local weights.

=== Residual Analysis

In order to substantiate @sec:2_5_commentary,
we perform a residual analysis on the prediction and forecasting data.

Recalling the assumptions of our model,
we expect the residuals to be normally distributed around zero:

$
vv(r) = #y - #y_pred ∼ 𝒩(0, σ^2)
$

It is readily apparent in @fig:plot_2_6_forecast_residuals that the residuals
are not normally distributed, nor do they average to zero:

$
  𝔼[vv(r)] = -8739
$

It should be noted, that the sum of the residuals is relatively close to zero,
which is by construction as can be seen in @sec:2_2_paramater_estimation.

We conclude that the model is not appropriate to describe the data.

#figure(
  image("output/plot_2_6_forecast_residuals.png"),
  caption: [
    Residual analysis on prediction and forecasting of total vehicle registrations in Denmark.
    Dashed black line delineates the transition from prediction to forecasting.]
) <fig:plot_2_6_forecast_residuals>

== WLS - Local Linear Trend Model

In order to mitigate some of the issues observed in the Linear Trend Model of @sec:2_linear_trend_model, we propose a _Weighted Least Squares_ model with local weights.

That is, more recent observations are weighted higher than older observations.

This is facilitated by repurposing the variance-covariance matrix of the residuals. For Ordinary Least Squares (OLS), we assumed the residuals to be modelled as $#ε ∼ vv(𝒩)(0, σ^2)$. That is, we have imposed an assumption of _homoscedasticity_ on the residuals #ε. In Weighted Least Squares (WLS) modelling, we relax this assumption and instead consider the variances of the residuals to be _heterogeneous_ – that is, the residuals are _heteroscedastic_:

$
  #ε ∼ vv(𝒩)(0, vv(σ)^2)
$

For Global Weighted Least Squares, one would ideally know the variances of the residuals a priori and simply weigh the residuals by the inverse of the variances.

We instead choose to _exploit_ the Weighted Least Squares model to introduce _locality_ in the estimation by letting the covariances of the residuals be modelled by $σ^2 #W$. We recall that the solution to the OLS problem (@theta_estimation) was obtained by optimization of the _residual sum of squares_ (RSS) between the observations and the model predictions.
If we choose to weigh the residuals by their recency, we able to obtain an
estimator that values recent information more highly than older information.

A concrete example of a #W that would induce locality would be one given by

$
#W = mat(delim: "[", lambda_(N), 0, 0, dots.h.c, 0; 0, lambda_(N - 1), 0, dots.h.c, 0; 0, 0, lambda_(N - 2), dots.h.c, 0; dots.v, dots.v, dots.v, dots.down, dots.v; 0, 0, 0, dots.h.c, lambda_(0))
$

If we set $λ_N, λ_(N-1) = 0$, we can appreciate that the model is now more local in that
the two oldest observations are not considered in the estimation. Typically one would choose #W to be a diagnoal matrix of exponentially decaying weights, as we will see in @sec:3_1_variance_covariance_and_weight_matrices.

Referencing @Madsen_2008[p.~38-39], we state without proof that the _normal equations_ for the WLS model are given by:

$
  (#X^T #W #X) #θ_pred = #X^T #W #y
$

Note that we have changed the notation slightly, denoting the weight matrix as #W rather than perverting the sigma as is done in @Madsen_2008.
That is, we have made the substitution $#W = mm(Σ)^(-1)$.

Assuming invertibility of $#X^T #W #X$, we obtain the parameter estimator:
$
  #θ_pred = (#X^T #W #X)^(-1) #X^T #W #y
$ <eq:3_theta_estimation>

We again follow the proof given in @sec:2_2_paramater_estimation, though noting that the variances of the residuals are now given by $𝔼[#ε #ε^T] = σ^2 #W$. As such, @eq:cov_theta_pred becomes:

$
  "Cov"[#θ_pred] = σ^2 (#X^T #W #X)^(-1)
$

A typical choice of the weight matrix would be to have exponentially decaying weights:

$
  #W = sum_(i=1)^N sum_(j=1)^N δ_(i j) λ^(N-i)
  wide wide δ_(i j) =
  cases(
    1 quad & i = j\
    0 quad & i ≠ j
  )
$

Where $δ_(i j)$ is _Kroneckers delta_, $λ$ is the _forgetting factor_, and $N$ denotes the dimension of #y.

To recap, we obtaing _locality_ in our model by gradually _forgetting_ older observations.

=== Variance-Covariance and Weight Matrices <sec:3_1_variance_covariance_and_weight_matrices>

Encouraging the reader to grit their teeth we summarise the 'variance-covariance' matrices for the OLS and WLS models as follows:

$
  mm(Σ) &= σ^2 mm(𝕀) wide && "OLS"\
  #W &= σ^2"diag"(vv(λ)) wide vv(λ) = mat(delim: "[", λ^(N-1), λ^(N-2), … , λ^0) wide && "WLS"\
$

Stating these more explicitly, they become:
$
  mm(Σ) &= mat(delim: "[",
  σ^2;
  ,⋱,;
  ,,σ^2;
  ) wide && "OLS"\

  #W &= mat(delim: "[",
  σ^2 λ^(N-1);
  ,σ^2 λ^(N-2);
  ,,⋱,;
  ,,,σ^2 λ^1;
  ,,,,σ^2;
  ) wide && "WLS"\
$

For the rest of the analysis, we will consider the _forgetting factor_ $λ ≡ 0.9$.

=== Weighting Regimen

Considering 72 time steps, we visualise the decay of the weights in @fig:plot_3_2_weights

#figure(
  image("output/plot_3_2_weights.png"),
  caption: [Decay of weights in the Weighted Least Squares model for $λ≡0.9$.]
) <fig:plot_3_2_weights>

To better visualise the weights, we plot the training data where the opacity of the individual observations is given by their weight in the WLS model parameterised by $λ=0.9$:

#figure(
  image("static_output/plot_3_2_decaying_weights.png", width: 50%),
  caption: [Forgetting in WLS Model visualised using training data set for $λ≡0.9$.]
) <fig:plot_3_2_decaying_weights>

=== Sums of weights

We find that the sum of the weights is given by a geometric sequence:

$
  sum_(i=1)^N λ^(N-i) = (1 - λ^N)/(1 - λ)
$

For the limit $N → ∞$ and $λ=0.9$, we find:

$
  lim_(N → ∞) sum_(i=1)^N λ^(N-i) = 10
$

We consider $N=72$ and find the truncation gives:

$
  sum_(i=1)^72 λ^(72-i) ≈ 9.995 wide "WLS"
$

The corresponding sum for the OLS model is clearly not bounded and just becomes:

$
  sum_(i=1)^N 1 = N wide "OLS"
$

For $N=72$, this simply becomes $72$.

=== Parameter Estimation

We can now estimate the parameters of the WLS model using @eq:3_theta_estimation:

$
  hat(θ)_1 &= (-1159200 ± 600) × 10^2\
  hat(θ)_2 &= (5173.5 ± 2.6) × 10^2\
$

Where the standard error is computed using Equation (3.43) in @Madsen_2008[p.~39].

=== Forecasting

We now perform a 12 month forecast using the WLS model in similar to fashion
to the one carried out in @sec:2_3_prediction, as seen in @fig:plot_3_5_forecast.

#figure(
  image("output/plot_3_5_forecast_wls_ols.png"),
  caption: [Comparison of the OLS and WLS ($λ=0.9$) forecasts
  of total number of registered vehicles in Denmark.]
) <fig:plot_3_5_forecast>

It should be noted that the standard error estimate used in the WLS model has
been corrected with respect to @eq:sigma for the truncating effect of the weights by calculating the
degrees of freedom using the _total memory_, $T$ of the model instead:

$
  σ = norm(#y - #X #θ_pred)^2/sqrt(T - p)
$ <eq:sigma_wls>

Where the _total memory_ is given by the sum of the weights:

$
  T = sum_(i=1)^N λ^(N-i) = (1 - λ^N)/(1 - λ)
$

For $λ<1$ and $N∈ℤ_+$ we observe that $T < N$, which agrees with intuition – a model with less memory _should_ be less _confident_ in its predictions,
all other things being equal.

Even with this correction reducing the model confidence,
we find an improved prediction and associated prediction interval
as seen in @fig:plot_3_5_forecast. This can be understood as the sum of squared residuals being reduced by the weighting scheme.

While the prediction does exhibit an improved fit to the test data
when compared to the prediction carried out with an OLS model,
it remains important to note that the model
may not generally be expected to forecast well, particularly in the event of
drastic changes in legislation or governance of motor driven vehicles.

Given a choice between the two models, short-term forecasts would be better served by the WLS model. It is not possible to confidently state which model would be preferable for long-term forecasting.

=== Different forgetting factors

Our choice of $λ=0.9$ was somewhat arbitrary and wholly unjustified.
In order to gauge which choice of forgetting factor may best fit out test data,
we run the WLS analysis repeatedly. It should be noted that the choice of forgetting factor will inevitably be a qualitative matter and there is no guarantee that the $λ$ that minimises the sum of squared residuals is necessarily a good or sensible choice for prediction for data sets that are not the current test set.

The outcome of the analysis can be seen in @fig:plot_3_6_forgetting_factors, where we have omitted the confidence interval estimates for clarity.

#figure(
  image("output/plot_3_6_forgetting_factors.png"),
  caption: [12 month forecasts of total number of registered vehicles in Denmark using Local Weighted Least Squares with a variety of forgetting factors.]
) <fig:plot_3_6_forgetting_factors>

We find that $λ=0.9$ is a reasonable choice for the test data set, though it is not necessarily the best choice for other data sets.

== Recursive Estimation and Optimization of $λ$ <sec:4_recursive_estimation>

The OLS problem whose parameters may be obtained as given by @theta_estimation can be reformulated as an _iterative_ problem:

For a timestep $t$ we find from @theta_estimation:

$
  #θ_pred _t &= (#X _t^T #X _t)^(-1) #X _t^T #y _t
  &= mm(R)^(-1)_t vv(h)_t
$ <eq:4_1_theta_pred>

Having used the definitions:
$
  mm(R) _t &≡ #X _t^T #X _t &&= sum_(i=1)^t vv(x)_i vv(x)_i^T\
  vv(h) _t &≡ #X _t^T #y _t &&= sum_(i=1)^t vv(x)_i y_i
$ <eq:4_1_R_h_sums>

=== Update Equations <sec:4_1_update_equations>

Inspection of @eq:4_1_R_h_sums readily reveals the objects at timestep $t$ as a function of their values at timestep $t-1$:
$
  mm(R) _t &= mm(R) _(t-1) + vv(x)_t vv(x)_t^T &&= sum_(i=1)^(t-1) vv(x)_i vv(x)_i^T + vv(x)_t vv(x)_t^T\
  vv(h) _t &= vv(h) _(t-1) + vv(x)_t y_t &&= sum_(i=1)^(t-1) vv(x)_i y_i + vv(x)_t y_t
$ <eq:4_1_R_h_recursive>

We can then rewrite @eq:4_1_theta_pred as:

$
  #θ_pred _t = #θ_pred _(t-1) + mm(R)^(-1)_t vv(x)_t (y_t - vv(x)_t^T #θ_pred _(t-1))
$ <eq:4_1_theta_pred_recursive>

#pagebreak()

Somewhat less neatly typeset, we also include the update equations as a pen-and-paper scan to complete the assignment objective.

#figure(
  image("static_output/jeppe_4_1_1.jpg", width: 80%),
  caption: [Update equations on pen and paper by Jeppe Klitgaard (s250250).]
) <fig:plot_4_1_recursive_estimation>

#figure(
  image("static_output/jeppe_4_1_2.jpg", width: 65%),
  caption: [Calculations of $R_1$ and $R_2$ on pen and paper by Jeppe Klitgaard (s250250).]
) <fig:plot_4_1_recursive_estimation>

#pagebreak()

#figure(
  image("static_output/yunis_4_1.svg"),
  caption: [Update equations on pen and paper by Yunis Wirkus (s250700).]
)
=== Computational Implementation <sec:4_2_computational_implementation>

We implement the Recursive Least Squares as given in @eq:4_1_R_h_recursive and @eq:4_1_theta_pred_recursive in Python as follows:

```python
def iterate(R, theta, x_t, y_t):
    R_t = R + x_t @ x_t.T
    theta_t = theta + (inv(R_t) @ x_t) * (y_t - x_t.T @ theta)

    return R_t, theta_t
```

Using initial conditions $mm(R)_0 = 0.1 ⋅ mm(𝕀)$ and $vv(θ)_0 = 0$, we iterate over the data set to $t'=t-3$ (not inclusive) and obtain the parameter estimates #θ_pred. These are presented in @fig:plot_4_2_rls_parameters.

#figure(
  image("output/plot_4_2_rls_parameters.png"),
  caption: [Parameter estimates obtained using Recursive Least Squares.]
) <fig:plot_4_2_rls_parameters>

We find that the initial parameter estimates are very poor,
but gradually improve as more data is included in the estimation.
Notably, we can see that as the slope, $θ_2$, increases as part of the iterative correction, the intercept, $θ_1$, shifts down.

While the matrix operations in @eq:4_1_theta_pred_recursive are somewhat opaque, careful inspection of the parameter update equation reveals a residual term $y_t - vv(x)_t^T #θ_pred _(t-1)$ where $vv(x)_t^T #θ_pred$ gives a prediction $hat(y)_t$ using the model found in the previous timestep.

The term $mm(R)_t^(-1) vv(x)_t$ then scales the residual appropriately, where we understand $mm(R)$ to take the role of $#X^T #X$ in the OLS model.

Now we will draw a comparison to the regular, non-recursive OLS
model. Imagine setting

$
ε_t = y_t - vv(x)_t^T #θ_pred _(t-1)
$

Which gives:

$
#θ_pred _t = #θ_pred _(t-1) + mm(R)^(-1)_t vv(x)_t ε_t
$

It becomes clear upon comparison with the regular OLS model that the term $vv(h)_t ≡ vv(x)_t ε_t$ may be interpreted as height of the parameter projection. Given
$upright(bold(R))_t = upright(bold(X))_t^T upright(bold(X))_t$ is a
projection of the design matrix $mm(X)_t$ into a symmetric, invertible
matrix, $vv(h)_t$ can be thought of as a scaled residual term.

Then overall the term becomes

$  & upright(bold(R))_t^(- 1) vv(x)_t (y_t - vv(x)_t^T #θ_pred _(t - 1))\
 & = upright(bold(R))_t^(- 1) vv(x)_t ε\
 & = upright(bold(R))_t^(- 1) vv(h)_t\
 & = upright(bold(hat(theta)))_t^(\*) &  & upright("from OLS")\
\
arrow.r.double upright(bold(hat(theta)))_t & = upright(bold(hat(theta)))_(t - 1) + upright(bold(hat(theta)))_t^(\*)
$

So we can think of the new $upright(bold(hat(theta)))_t$ as an updated
parameter estimate, based on the previous estimate _and_ a correction term arising from the residual of the current observation compared against a prediction of this observation using the previous model.

=== Calculation and Comparison of Parameter Estimates <sec:4_3_parameter_estimation>

We now calculate the estimates $#θ_pred _N$.
The meaning of $#θ_pred _N$ is ambiguous in the task.
We interpret it to mean $hat(θ) _1$ and $hat(θ) _2$ at time $t'=t$.

This gives us:

#figure(
  table(
    columns: 3,
    align: center,
    table.header([*Model*], [*$hat(θ)_1$*], [*$hat(θ)_2$*]),
    [*OLS*], [$-110.365428 × 10^6$], [$56.14 × 10^3$],
    [*RLS*], [$-58.32 × 10^3$], [$1.568 × 10^3$],
  ),
  caption: [Comparison of parameter estimates obtained using Ordinary Least Squares (OLS) and Recursive Least Squares (RLS).]
) <table:4_3_comparison_parameters>

From @table:4_3_comparison_parameters,
we see that the RLS model has produced a significantly different parameter estimate than the OLS model.
In order to gauge the performance of the two models, we compare the predictions in @fig:plot_4_3_ols_rls.

#figure(
  image("output/plot_4_3_ols_rls.png"),
  caption: [Comparison of the OLS and RLS forecasts of total number of registered vehicles in Denmark.]
) <fig:plot_4_3_ols_rls>

It is readily apparent that the RLS model is a poor fit to the data,
which we attribute to a very poor initial guess of the parameters $vv(θ)$.

From the intuition we gathered in @sec:4_2_computational_implementation,
we understand that $mm(R)^(-1)_t$ determines the scale of updates to parameters $#θ_pred _t$.
We can allow for larger updates by reducing the initial information matrix $mm(R)_t$,
though notably it cannot be $mm(0)$ as this would not be invertible.
Additionally, we can improve the performance of our RLS model by having a better estimate for $vv(hat(θ))_0$
.

Based on this, we use the OLS model to come up with a reasonable $vv(hat(θ))_0$ and set a smaller initial information matrix $mm(R)_0$:
$
  mm(R)_0 &= mm(𝕀) ⋅ 10^(-8)\
  vv(hat(θ))_0 &= vec(delim: "[", -100 × 10^6, 50 × 10^3)
$

This increases the initial 'eagerness' of the model to update the parameters when fed new information, while also giving a decent guess for the initial parameters.

In order to get an intuition for the RLS model, we now reproduce @fig:plot_4_2_rls_parameters with the updated initial conditions:

#figure(
  image("output/plot_4_3_good_rls_parameters.png"),
  caption: [Parameter estimates obtained using Recursive Least Squares with updated initial conditions.]
) <fig:plot_4_3_good_rls_parameters>

Reproducing @fig:plot_4_3_ols_rls with the updated initial conditions, we observe a significantly improved fit to the data:

#figure(
  image("output/plot_4_3_good_ols_rls.png"),
  caption: [Comparison of the OLS and RLS forecasts of total number of registered vehicles in Denmark with updated initial conditions.]
) <fig:plot_4_3_good_ols_rls>

From @fig:plot_4_3_good_ols_rls, we can realise that the recursive model could be improved with forgetting
as the parameters of a linear model are not constant with respect to time for the given data set.

=== Recursive Estimation of Forgetting Factor <sec:4_4_recursive_estimation_forgetting>

We first seek to get some intuition for the recursive model with _forgetting_.
Having realised previously that the $mm(R)_t$ matrix represents some _knowledge_ of the observations accumulated over past iterations, we can understand the _forgetting factor_ as a way to _decay_ this knowledge.

As such, we introduce the _forgetting factor_, $λ_t$:

$
  mm(R)_t = λ_t mm(R)_(t-1) + vv(x)_t vv(x)_t^T
$ <eq:4_4_R_forgetting>

It should be noted that one usually choses $λ_t$ to be a constant and it should not be confused with the power series that appears along the diagonal in Weighted Least Squares with exponential forgetting.
The exponential nature of the decay of knowledge is encoded into the iterative algorithm, as can be realised by inspection when expanding $mm(R)_(t-1)$ in @eq:4_4_R_forgetting.

Using the RLS model as given by the update equations in @sec:4_1_update_equations with the forgetting factor modification in @eq:4_4_R_forgetting, we can now estimate the parameters $θ_1$, $θ_2$ as a function of time $t$ for forgetting factors $λ ∈ {0.6, 0.99}$.

#figure(
  image("output/plot_4_4_rlsf_2_params_1.png"),
  caption: [Estimates of parameter $θ_1$ obtained using Recursive Least Squares with forgetting factors $λ ∈ {0.6, 0.99}$.]
) <fig:plot_4_4_rlsf_2_params_1>

#figure(
  image("output/plot_4_4_rlsf_2_params_2.png"),
  caption: [Estimates of parameter $θ_2$ obtained using Recursive Least Squares with forgetting factors $λ ∈ {0.6, 0.99}$.]
) <fig:plot_4_4_rlsf_2_params_2>

Inspecting @fig:plot_4_4_rlsf_2_params_1 and @fig:plot_4_4_rlsf_2_params_2,
we observe higher variation in the parameter estimates for $λ=0.6$ than for $λ=0.99$.
We can understand this by appreciating the increase in _locality_ that arises from a lower forgetting factor.
As expected, the more local model varies its parameter estimates more rapidly than a less local model.
Both models have been initialised with $mm(R)_0 = 𝕀 ⋅ 10^(-8)$ and $vv(θ)_0 = vv(0)$.

#let math_eval(x) = eval(x, mode: "markup")
#let wls_rls_params = csv("output/wls_rls_params.csv")
#figure(
  table(
    columns: 4,
    table.header(..wls_rls_params.remove(0).map(math_eval).map(strong)),
    ..wls_rls_params.flatten().map(math_eval)
  ),
  caption: [
    Parameter estimates using WLS and RLS models with different forgetting factors.
    Note that the parameters for the RLS model are for $t=N$.
  ]
) <table:4_4_wls_rls_params>

From @table:4_4_wls_rls_params, we observe that the RLS model with $λ=0.7$ produces final parameter estimates that are very close to the WLS model.

=== One-Step Predictions with RLS <sec:4_5_one_step_predictions>

We introduce the concept of a prediction horizon, $k$, which denotes the number of steps ahead of given data set we wish to predict. To ease notation, we employ $t+k|t$ to mean a prediction for time $t+k$ given data available at time $t$. As such:

$
  hat(y) _(t+k|t) = vv(x)_(t+k|t)^T #θ_pred _t
$

We may then define the residual as:

$
  hat(ε)_(t+k|t) = y_(t+k) - hat(y) _(t+k|t)
$

We employ the RLS model with forgetting factor $λ ∈ $ to produce residuals against the
training data set for one-step predictions ($k ≡ 1$) and present the results in @fig:plot_4_5_rls_residuals.

#figure(
  image("output/plot_4_5_rls_residuals.png"),
  caption: [
    Residuals obtained using Recursive Least Squares with forgetting factors $λ ∈ {0.6, 0.99}$ for one-step predictions.
  ]
) <fig:plot_4_5_rls_residuals>

We observe that the model with shorter memory ($λ=0.7$) is generally better at adaption to the
changes in trend found in data set.
Again, this matches with our expectations that a more local model will be able to better match the local trends in the data.
As a result, the magnitude of the residuals is generally lower for the $λ=0.7$ model.

=== Optimal Forgetting Factor for a Prediction Horizon <sec:4_6_optimal_forgetting>

We now seek to understand how the optimal forgetting factor $λ^*$ arises for a given prediction horizon characterised by a $k$-step prediction.

We consider the instructive cases of the extrema, $λ→0$ and $λ→1$. For $λ→0$, we have a model that is _very_ local and forgets past information almost immediately.
We would expect this to perform relatively poorly, even for a somewhat non-linear model.
For $λ→1$, we have a model that is _very_ global and retains past information almost indefinitely,
will perform well for a truly linear model, but for coloured noise, it will perform poorly.

Appreciating these two logical extrema,
we can begin to see how our forgetting factor should be matched to the horizon we are trying to predict. If we just want to predict the next timestep, it is intuitive that a more local model would in generally perform better since longer-term non-linear trends will not have a significant impact on the prediction. On the other hand, if we are trying to predict further into the future,
an overly local model would not have sufficient understanding of longer-term trends that may dominate the residual for a prediction horizon of this width.

To investigate this, we run the RLS model with forgetting factors $λ ∈ {0.50, 0.51, …, 0.99}$ for prediction horizons $k ∈ {1, 2, … , 12}$. and collect the Root Mean Square Error as a function of $k$ and $λ$:

$
  "RMSE"(k, λ) = sqrt(1/(N - k) sum_(t=1)^(N-k) hat(ε)_(t+k|t)^2(λ))
$

Where $hat(ε)_(t+k|t)(λ)$ denotes the residuals for an RLS model with forgetting factor $λ$ and prediction horizon $k$.

#figure(
  image("output/plot_4_6_rmsek.png"),
  caption: [
    Root Mean Square Error as a function of prediction horizon $k$ and forgetting factor $λ$.
    The minimal RMSE for each $k$ is marked with a dashed vertical line of the corresponding colour.
  ]
) <fig:plot_4_6_rmsek>

From @fig:plot_4_6_rmsek, we observe that the optimal forgetting factor for the very highest value of $k=12$ is $λ=0.99$, while shorter prediction horizons are better served by a more local models. A choice of around $λ=0.65$ is generally a good choice across different horizons.

Taking note of the shape of the curves in @fig:plot_4_6_rmsek, we can validate our speculation above. For short prediction horizons, the optimal forgetting factor is lower, while for longer prediction horizons, the optimal forgetting factor is higher.

=== Predictions Using RLS and Other Models
<sec:4_7_predictions_model_comparison>

#figure(
  image("output/plot_4_7_ols_ewls_rlsf_comparison.png"),
  caption: [
    Prediction of total car sales, compared by model (OLS, WLS, RLS) for horizon $k=12$ and forgetting factor $λ=0.57$.
  ]
) <fig:plot_4_7_ols_ewls_rlsf_comparison>

Here the weight factor $λ=0.57$ is chosen loosely based on the previous @fig:plot_4_6_rmsek, as an estimated good value (in the middle of optimal values for $k \in {1,...,11}$). Therefore, it is not surprising that for $k=12$ the RLS-F (RLS with forgetting) does perform very poorly. For the chosen time horizon, the RLS-F model would perform best for $λ \to 1$, at which point it basically becomes an OLS model, better at predicting more global trends.

In this szenario the E-WLS (exponential/local WLS) model performs best, simply by visual analysis. Yet, all models perform quite poorly, not finding a good compromise between local trends and global trend.
This could be overcome by different choices of $λ$. However, the general notion of "forgetting" models (such as E-WLS and RLS-F) being better at predicting short term horizons (e.g. $k \in {1,2,3,4}$) remains. Nonetheless, this is highly dependent on noise and variance in the given time-series.

In the short following exploration, we investigate superficially, whether the choice of $λ$ fundamentally changes that notion or just marginally improves the models predictive power.
This time, we consider a range for $λ$, where it does not have its absolute minimum, but given @fig:plot_4_6_rmsek for $k=12$, the RMSE is close to minimal values.

#figure(
  image("output/plot_4_7_ols_ewls_rlsf_lambda_comparison.png"),
  caption: [
    Zoomed in prediction graph of total car sales, compared by model (OLS, WLS, RLS) for horizon $k=12$ and forgetting factor $λ \in [0.5, 0.75]$.
  ]
) <fig:plot_4_7_ols_ewls_rlsf_lambda_comparison>

As we can see in @fig:plot_4_7_ols_ewls_rlsf_lambda_comparison, the RLS-F with high "forgetting" (low $λ$) does capture extremely short-term trends; in fact, goes opposite the global trend. For higher $λ$ it approximates the WLS model, for $λ > 0.87$ (according to @fig:plot_4_6_rmsek) it will approximate the OLS solution.
Overall, the predictive power does improve for certain parameter $λ$, but only for select time horizons $k \in {1,2}$.



=== Reflexions on time adaptive models
<sec:4_8_reflexions>


1. **Overfitting VS Underfitting**

The problem with these simple linear models is definitely underfitting. All constructed models performed somewhat poorly in explaining the data and predicting it.
There would be mainly 3 knobs to consider:
- increase the amount of parameters for the model (thus making it a multivariate linear regression; resulting in hyperplanes as predictors) - the dataset has more to offer
- increase the model complexity to for example polynomial models - however, these come with different problems (especially for time-series; asymptotic behaviour towards boundries)
- hyperparameter tuning of $λ$ or even the structure of the entire weight matrix $Σ$


2. **Test Set Challenges**

There are challenges in time-series due to covariance and independence assumptions. So far, each observation is treated as an independent random variable (with the exception of an exponential weight matrix, which somewhat introduces the notion of more recent samples being more important than older ones). 

    However, time series can present trend-changes (as is the example here) or even seasonality behaviours. This is different from a dataset, that simply maps one variable or event to another, where order of the events is not relevant.
    
    Thus, the choice of split of a test da can dramatically impact the model's predictive capability. The time frame chosen can include of not include trend changes, that may or may not be present in the other split of data.

    One way to handle this would be by cross-validation. A common method, that would also be applicable to time-series is an adjusted K-Fold cross-validation. Basically, creating different models by sequentially leaving out some parts of the training data and then weighting by the confidence of each model. Bootstrapping is another famous method but due to its nature not applicable.


3. **Recursive Estimation as a Saviour?**

Actually the step-wise prediction (for only the test dataset though) was supperior to global model parameters.
Yet, it was not applicable or helpful for **actual** prediction, as the time horizon to predict into the future may vary and the last parameter estimation cuts of at the very end of the training dataset. Hence it is does not have the advantages of continually updating parameters any longer.

    However, in general the RLS models suffer from severe 'burn-in', initial guess uncertainty and more delicate hyperparameter tuning. Additionally, they are more computationally expensive.
    These properties make then not a good fit so far, but it might depend on the origin of the time series.


4. **Additional adaptive estimation techniques**

Some ideas for adaptive (maybe recursive) estimation on time series:
- Bayesian models, taking new information to update likelihood and posterior; would also allow for a good initial estimate depending on the origin of the data for the prior distribution
- momentum models; taking sudden large jumps in time-series as a penalty to a weight, not to update much -> would be great for noisy time series
- averaging different weighted models to capture notions of local trend changes and reversion to a global trend
- dynamic weight models; adjust the size of a WLS dynamically/recursively, based on the variance, noise or momentum of data
- recurrent neural-network based models -> most often a good predictor but hard to interpret or modify without continuous re-training (very computationally expensive)
- regularized models: Ridge-regression (L2) or ISTA (L1), the latter does not work well on most time-series because the L1 norm penalizes mostly "edges"


5. **Additional Comments**

- RLS models are quite the pain in the a** to get right, especially choosing initial values for $θ$, $R_t$ and $λ$
- a huge pitfall is to think that: modelling a time series, that is (visually) very obviously not showing a linear behavior, with a simple linear model...is a good idea
- as for the initial description and visual analysis, it became apparent that there are significant trend changes present; so much in fact, that one could possibly also just take the last few data points to predict and save oneself all the hustle - a wiser choise of train/test split would significantly improve the models


#bibliography("report.bib")
