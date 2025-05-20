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
#let wider = h(3em)

///////////////////////////////////////////////////////////////////////////////

#align(center, text(17pt)[
  *02417 Times Series Analysis - Assignment 3*
])

Authors:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-04-21

#pagebreak()

= Parameter Estimation in State-Space Model <sec:1>

We are given a state-space model represented by the scalar system:

$
  X_t = A X_(t-1) + B + e_(1,t)
$

Where:
- $A$ is the scalar state transition coefficient
- $B$ is a scalar bias term
- $e_(1,t) ∼ 𝒩(0, σ_1^2)$ is a Gaussian noise term

The associated observation model is given as:
$
  Y_t = X_t + e_(2,t)
$

Where $e_(2,t) ∼ 𝒩(0, σ_2^2)$ is another Gaussian noise term associated with observations.

== Realisations of State Vector <sec:1.1>
We are asked to perform 5 independent realisations of the state-space model with parameters $A = 0.9, B = 1, σ_1^2 = 1$ using the initial state $X_0 = 5$. The length of each realisation is $n = 100$.

This is done using Python and NumPy, with the relevant coding being found in `1.ipynb`.

#figure(
  image(
    "output/1_1_realisations.png",
  ),
  caption: [Realisations of the state-space model with $A = 0.9, B = 1, σ_1^2 = 1$ and $X_0 = 5$],
) <fig:1.1>

The 5 independent realisations of the state can be seen in @fig:1.1. Note that
we show the state vector $X_t$ as opposed to the observation vector $Y_t$. While the assignment is somewhat ambiguous regarding whether the desired realisations are those of the state vector or the observation vector, we reasonably assume it to be the state vector given that the parameter $σ_2$ which is required to generate the observation vector is not given..

== Realisation of Observation Vector <sec:1.2>

We now generate another realisation, this time also calculating the observation vector $Y_t$, which suffers from additional _observation noise_ $e_(2,t)$. The parameters are the same as before, with the addition of $σ_2^2 = 1$:

$
X_t &= a X_(t - 1) + b + e_(1,t) &&wider e_(1,t) ~ 𝒩(0, sigma_1^2) \
Y_t &= X_t + e_(2,t) &&wider e_(2,t) ~ 𝒩(0, sigma_2^2 = 1),
$ <eq:1.2>

This yields @fig:1.2, in which we can see both the latent state vector $X_t$ and the observation vector $Y_t$, which clearly is affected by the observation noise.

#figure(
  image(
    "output/1_2_realisation.png",
  ),
  caption: [Realisation of the latent state vector $X_t$ and the associated observation vector $Y_t$ with $A = 0.9, B = 1, σ_1^2 = 1, σ_2^2 = 1$ and $X_0 = 5$],
) <fig:1.2>

We note in @fig:1.2 that the observation vector appears to hover around the latent state vector with a symmetric residual, as would be expected from @eq:1.2.

This corresponds to a real-world scenario where the act of _observing_ something will itself introduce a noise that is independent of the underlying process.

== Kalman Filter <sec:1.3>

Next we implement a simple Kalman Filter using the function template provided in `kalmanfilter.R`. We convert the function to Python and fill in the blank assignment operations as follows:

$
  "State" &&wider hat(X)_(1|0) &= X_"prior" &&wider (10.79)\
  "State Variance" &&wider Σ_(1|0)^(x x) &= P_"prior" &&wider (10.80)\

  "State" &&wider hat(X)_(t+1|t) &= A hat(X)_(t|t) + B &&wider (10.63)\
  "State Variance" && Var[hat(X)_(t+1|t)] &= Var[tilde(X)_(t+1|t)] = A^2 Var[hat(X)_(t|t)] + sigma_1^2 &&wider (10.54), (10.67)\
  \
  "Innovation" &&wider hat(Y)_(t+1|t) &= C hat(X)_(t+1|t) &&wider (10.64)\
  "Innovation Variance" &&wider Var[tilde(Y)_(t+1|t)] &= C^2 Var[tilde(X)_(t+1|t)] + R &&wider (10.68)\
  \
  "Kalman Gain" &&wider K_t &= C Var[tilde(X)_(t+1|t)] / Var[tilde(Y)_(t+1|t)] &&wider (10.75)\
  \
  "Filtered State" &&wider hat(X)_(t|t) &= hat(X)_(t|t-1) + K_t (Y_t - C hat(Y)_(t|t-1)) &&wider (10.73)\
  "Filtered State Variance" &&wider Var[hat(X)_(t|t)] &= (1 - K_t) Var[hat(X)_(t|t-1)] &&wider (10.74)\
$

Where rather than regurgitating the lengthy derivations, we refer to the relevant equations in the course textbook @Madsen_2008[Chapt.~10].

Implementing this, we are able to use the Kalman filter on a realisation similar to that outlined in @sec:1.2 and shown in @fig:1.2.

Helpfully, our Kalman Filter implementation already uses the latent state variance, which makes it particularly simple to compute a 95% confidence interval overlaid on the predicted state in @fig:1.3

#figure(
  image(
    "output/1_3_kalman_filter.png",
  ),
  caption: [
    Scalar Kalman Filter prediction with 95% confidence intervals.
    We have again used parameters $A = 0.9, B = 1, σ_1^2 = 1, σ_2^2 = 1, X_0 = 5$ and use the process given in @eq:1.2.
    ],
) <fig:1.3>

We find that the filter is able to predict the latent state using the observations well, noting that the noise variances are passed in explicitly as opposed to determined by the filter itself on the supplied data.

While further analysis may be performed, we note by inspection of @fig:1.3 that the confidence intervals appear to stabilise quickly after the first few observations, which is to be expected as the Kalman filter converges to the true state, which here is described by a stationary process.

The prediction intervals necessarily reach the observed size due to the inherent noise of underlying process and the observation noise. As such, the confidence intervals may not be reduced further by using a more complex model,
as they are an irreducible property of the underlying process.

They may, however, be reduced somewhat by decreasing the observation noise, which may be possible in a real-world scenario by improving the experimental design, for example by using more precise expermental equipment.

== Maximum Likelihood Framework <sec:1.4>

We now move into the Maximum Likelihood Framework, in which we formulate
a _likelihood function_ over the prediction of the Kalman Filter.

Building on the approach and functions derived in @sec:1.3 and using the template function provided in `myLogLikFun.R`, we seek a Python implementation of the log-likelihood function.

By using the prior and posterior distributions of the state vector, we are able to derive the likelihood distribution as @Madsen_2008[Sec.~10.3.3]:
$
  L = (tilde(Y)_(t+1|t) |X_(t+1), 𝒴_t))
  &∼ 𝒩(C(X_(t+1) - A hat(X)_(t|t) - B), σ_2)\
  &∼ 𝒩(hat(X)_(t+1|t+1), σ_2)\
  &= 1/(σ_2 sqrt(2π)) exp(- (Y_(t+1) - hat(X)_(t+1|t+1))^2 / σ^2)
$

From which we are able to construct the log-likelihood:
$
  log L = -log(σ_2) -log(2π)/2 - (Y_(t+1) - hat(X)_(t+1|t+1))^2/σ_2^2
$

Using this, we are able to estimate parameters using the Kalman Filter framework
through the maximum likelihood formulation based on observations.

The associated code may be found in the attached Jupyter Notebook named `1.ipynb`.

Using the implementations from previous sections, we construct $N=100$ realisations
of the observation vector with $n=100$ samples in each and subsequently employ SciPy's optimisation suite to
find parameters that maximize the log-likelihood function through minimisation of the negative log-likelihood.

Summary statistics are then calculated on the obtained parameter estimates to give
insight into the estimation process.

We perform a round of simulations and estimations with parameters $a=1, b=0.9, σ_1=1$,
which yields @fig:1.4.1.

#figure(
  image(
    "output/1_4_experiment_a_0.9_b_1_sigma_1_1.png",
    width: 45%
  ),
  caption: [Parameter estimates for $A, B, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $A = 0.9, B = 1, σ_1^2 = 1$],
) <fig:1.4.1>

We perform another two rounds of these simulations to produce @fig:1.4.2 and @fig:1.4.3, which show the parameter estimates for $A, B, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $A = 0.9, B = 0.9, σ_1^2 = 1$ and $A = 1, B = 0.9, σ_1^2 = 5$, respectively.

#grid(
  columns:2,
  [#figure(
    image(
      "output/1_4_experiment_a_0.9_b_0.9_sigma_1_1.png",
    ),
    caption: [Parameter estimates for $A, B, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $A = 0.9, B = 0.9, σ_1^2 = 1$],
  ) <fig:1.4.2>],
  [#figure(
    image(
      "output/1_4_experiment_a_1_b_0.9_sigma_1_5.png",
    ),
    caption: [Parameter estimates for $A, B, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $A = 1, B = 0.9, σ_1^2 = 5$],
  ) <fig:1.4.3>]
)

Comparing across the three figures, we find that the Kalman Filter does a fairly good job of estimating the parameters, although a large variance in the estimate of the bias term $B$ is observed. We can understand this intuitively, as the observable effect of the bias term, which simply shifts the state vector and thus the observation vector, is difficult to distinguish from a shift arising from the _random walk_ arising from the integration of the noise term $e_(1,t)$.

This is particularly true when the noise term is large, as in @fig:1.4.3, where we observe a longer-tailed distribution of the bias term $B$.

We find that for all three rounds and across all three parameters, the Kalman Filter is able to estimate the parameters with a reasonable degree of accuracy, although the bias term $B$ is difficult to estimate accurately.


#bibliography("report.bib")
