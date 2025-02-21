#import "@preview/physica:0.9.4": super-T-as-transpose

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

= Assignment 1

== Plot Data

=== Construct time variable

=== Plotting observations

#image("output/plot_observations.png")

== Linear Trend Model

We are given the General Linear Model (GLM) in sloppy notation:
$
  Y_t = θ_i + θ_2 ⋅ x_t + ϵ_t
$ <GLM>

=== Matrix Form

We rewrite @GLM as a matrix:

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

=== Parameter Estimation

We can estimate the parameters $#θ$ by deriving the _normal equations_:

$
  mref("GLM_matrix") quad ⇒ quad
  #X^T #y = #X^T #X #θ_pred quad ⇒ quad
  #θ_pred = (#X^T #X)^(-1) #X^T #y
$ <theta_estimation>

Where $#X$ has assumed to be invertible.

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
$

Where we understand the variances of the predicted variables to be
the diagonal elements of the covariance matrix:

$
  "Var"[#θ_pred]
  &= "diag"("Cov"[#θ_pred])
$

Which gives the following standard errors for the parameter estimates:

$
  σ_(hat(θ_i)) = sqrt("Var"[hat(θ_i)]) = sqrt(σ^2 (#X^T #X)^(-1)_(i)) wide wide i ∈ {1, 2}
$

Computing these for the entire training dataset gives:

$
  hat(θ_1) &= (-110360 ± 60) × 10^3\
  hat(θ_2) &= (3593.6 ± 1.8) × 10^3\
$

Using these predicted parameters,
we are able to produce an estimate of the vehicle registrations for the training dataset
using the General Linear Model as seen in @fig:plot_observations_predicted.

#figure(
  image("output/plot_observations_predicted.png"),
  caption: [Estimation of vehicle registrations using the General Linear Model.]
) <fig:plot_observations_predicted>

=== Prediction

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
It should be noted that the estimation of the confidence interval valid only
in the case where the observations are described by the General Linear Model.

#let forecast_2_3 = csv("output/forecast.csv")
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
Additionally we would expect supply chain disruptions to have strongly influence
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

As such, we conclude that the model is not appropriate to describe the data.

#figure(
  image("output/plot_2_6_forecast_residuals.png"),
  caption: [
    Residual analysis on prediction and forecasting of total vehicle registrations in Denmark.
    Dashed black line delineates the transition from prediction to forecasting.]
) <fig:plot_2_6_forecast_residuals>

#bibliography("report.bib")
