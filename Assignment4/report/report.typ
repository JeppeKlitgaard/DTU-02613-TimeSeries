#import "@preview/physica:0.9.4": super-T-as-transpose
#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#import "utils.typ": mathformatter, mref
#import "@preview/wrap-it:0.1.1":wrap-content

#set page(
  paper: "a4",
  numbering: "1/1",
  header: context {
    if(counter(page).get().at(0)== 1) [

    ] else [
    #set align(center)
    #set text(10pt, fill: gray)
      02417 Times Series Analysis - Assignment 4
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
  *02417 Times Series Analysis - Assignment 4*
])

Authors:
- Jeppe Klitgaard `<s250250@dtu.dk>`
- Yunis Wirkus `<s250700@dtu.dk>`

Date: 2025-05-23

// #pagebreak()

= Parameter Estimation in State-Space Model <sec:1>

We are given a state-space model represented by the scalar system:

$
  X_t = a X_(t-1) + b + e_(1,t)
$ <eq:1.1>

Where:
- $a$ is the scalar state transition coefficient
- $b$ is a scalar bias term
- $e_(1,t) ∼ 𝒩(0, σ_1^2)$ is a Gaussian noise term

The associated observation model is given as:
$
  Y_t = X_t + e_(2,t)
$

Where $e_(2,t) ∼ 𝒩(0, σ_2^2)$ is another Gaussian noise term associated with observations.

== Realisations of State Vector <sec:1.1>
We are asked to perform 5 independent realisations of the state-space model with parameters $a = 0.9, b = 1, σ_1^2 = 1$ using the initial state $X_0 = 5$. The length of each realisation is $n = 100$.

This is done using Python and NumPy, with the relevant coding being found in `1.ipynb`.

#figure(
  image(
    "output/1_1_realisations.png",
  ),
  caption: [Realisations of the state-space model with $a = 0.9, b = 1, σ_1^2 = 1$ and $X_0 = 5$],
) <fig:1.1>

The 5 independent realisations of the state as given by @eq:1.1 can be seen in @fig:1.1. Note that
we show the state vector $X_t$ as opposed to the observation vector $Y_t$. While the assignment is somewhat ambiguous regarding whether the desired realisations are those of the state vector or the observation vector, we reasonably assume it to be the state vector given that the parameter $σ_2$ which is required to generate the observation vector is not given.

== Realisation of Observation Vector <sec:1.2>

We now generate another realisation, this time also calculating the observation vector $Y_t$, which suffers from additional _observation noise_ $e_(2,t)$. The parameters are the same as before, with the addition of $σ_2^2 = 1$:

$
X_t &= a X_(t - 1) + b + e_(1,t) &&wider e_(1,t) ~ 𝒩(0, σ_1^2) \
Y_t &= X_t + e_(2,t) &&wider e_(2,t) ~ 𝒩(0, σ_2^2 = 1),
$ <eq:1.2>

This yields @fig:1.2, in which we can see both the latent state vector $X_t$ and the observation vector $Y_t$, which clearly is affected by the observation noise.

#figure(
  image(
    "output/1_2_realisation.png",
  ),
  caption: [Realisation of the latent state vector $X_t$ and the associated observation vector $Y_t$ with $a = 0.9, b = 1, σ_1^2 = 1, σ_2^2 = 1$ and $X_0 = 5$],
) <fig:1.2>

We note in @fig:1.2 that the observation vector appears to hover around the latent state vector with a symmetric residual, as would be expected from @eq:1.2.

This corresponds to a real-world scenario where the act of _observing_ something will itself introduce a noise that is independent of the underlying process.

== Kalman Filter <sec:1.3>

Next we implement a simple Kalman Filter using the function template provided in `kalmanfilter.R`. We convert the function to Python and fill in the blank assignment operations as follows:

$
  "State" &&wider hat(X)_(1|0) &= X_"prior" &&wider (10.79)\
  "State Variance" &&wider Σ_(1|0)^(x x) &= P_"prior" &&wider (10.80)\

  "State" &&wider hat(X)_(t+1|t) &= a hat(X)_(t|t) + b &&wider (10.63)\
  "State Variance" && Var[hat(X)_(t+1|t)] &= Var[tilde(X)_(t+1|t)] = a^2 Var[hat(X)_(t|t)] + σ_1^2 &&wider (10.54), (10.67)\
  \
  "Innovation" &&wider hat(Y)_(t+1|t) &= C hat(X)_(t+1|t) &&wider (10.64)\
  "Innovation Variance" &&wider Var[tilde(Y)_(t+1|t)] &= C^2 Var[tilde(X)_(t+1|t)] + R &&wider (10.68)\
  \
  "Kalman Gain" &&wider K_t &= C Var[tilde(X)_(t+1|t)] / Var[tilde(Y)_(t+1|t)] &&wider (10.75)\
  \
  "Filtered State" &&wider hat(X)_(t|t) &= hat(X)_(t|t-1) + K_t (Y_t - C hat(Y)_(t|t-1)) &&wider (10.73)\
  "Filtered State Variance" &&wider Var[hat(X)_(t|t)] &= (1 - K_t) Var[hat(X)_(t|t-1)] &&wider (10.74)\
$ <eq:1.3_Kalman_filter>

Where rather than regurgitating the lengthy derivations, we refer to the relevant equations in the course textbook @Madsen_2008[Chapt.~10].

Implementing this, we are able to use the Kalman filter on a realisation similar to that outlined in @sec:1.2 and shown in @fig:1.2.

Helpfully, our Kalman Filter implementation already uses the latent state variance, which makes it particularly simple to compute a 95% confidence interval overlaid on the predicted state in @fig:1.3

#figure(
  image(
    "output/1_3_kalman_filter.png",
  ),
  caption: [
    Scalar Kalman Filter prediction with 95% confidence intervals.
    We have again used parameters $a = 0.9, b = 1, σ_1^2 = 1, σ_2^2 = 1, X_0 = 5$ and use the process given in @eq:1.2.
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
  &∼ 𝒩(C(X_(t+1) - a hat(X)_(t|t) - b), σ_2)\
  &∼ 𝒩(hat(X)_(t+1|t+1), σ_2)\
  &= 1/(σ_2 sqrt(2π)) exp(- (Y_(t+1) - hat(X)_(t+1|t+1))^2 / σ^2)
$ <eq:1.4_implicit_Kalman_likelihood>

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
    "output/1_4_experiment_N100_a1_b0.9_σ_11.png",
    width: 50%
  ),
  caption: [Parameter estimates for $a, b, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $a = 1, b = 0.9, σ_1^2 = 1$],
) <fig:1.4.1>

We perform another two rounds of these simulations to produce @fig:1.4.2 and @fig:1.4.3, which show the parameter estimates for $a, b, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $a = 0.9, b = 0.9, σ_1^2 = 1$ and $a = 1, b = 0.9, σ_1^2 = 5$, respectively.

#grid(
  columns:2,
  [#figure(
    image(
      "output/1_4_experiment_N100_a0.9_b0.9_σ_11.png",
    ),
    caption: [Parameter estimates for $a, b, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $a = 0.9, b = 0.9, σ_1^2 = 1$],
  ) <fig:1.4.2>],
  [#figure(
    image(
      "output/1_4_experiment_N100_a1_b0.9_σ_15.png",
    ),
    caption: [Parameter estimates for $a, b, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $a = 1, b = 0.9, σ_1^2 = 5$],
  ) <fig:1.4.3>]
)

Note that $a=0.9$ was used for @fig:1.4.2 as opposed to the $a=5$ given in the assignment, which was reportedly a typo, as announced on the EdStem Forum by Justinas.

Comparing across the three figures, we find that the Kalman Filter does a fairly good job of estimating the parameters, although a large variance in the estimate of the bias term $b$ is observed. We can understand this intuitively, as the observable effect of the bias term, which simply shifts the state vector and thus the observation vector, is difficult to distinguish from a shift arising from the _random walk_ arising from the integration of the noise term $e_(1,t)$.

This is particularly true when the noise term is large, as in @fig:1.4.3, where we observe a wider distribution of the bias term $b$.

We find that for all three rounds and across all three parameters, we are able to estimate the parameters with a reasonable degree of accuracy, though we note that
the distribution of the estimates is quite large and the true value is only well-approximated when averaging across multiple realisations. We find that $a$ is the easiest parameter to estimate, while $b$ is the most difficult, as discussed above.

Alternatively, fewer or a single realisation may be used if the number of observations $N$ is sufficiently large.

== Non-Gaussian Noise <sec:1.5>

Our Kalman Filter was derived using the assumption of Gaussian noise, which
is often reasonable in practice due to the manifestation of the Central Limit Theorem. However, there are many cases where the underlying process may not be Gaussian, which we seek to explore by employing the Kalman Filter without modification on non-Gaussian processes.

To this end, we make a number of realisations using the Student's t-distribution in place of a normal distribution. Notably, we let the observation noise $e_(2,t)$ remain normally distributed, as this will often be the case in practice per the
Central Limit Theorem.

As such, our new process is given as:
$
  X_t &= a X_(t - 1) + b + e_(1,t) &&wider e_(1,t) ~ λ_t (ν, σ_1^2) \
  Y_t &= X_t + e_(2,t) &&wider e_(2,t) ~ 𝒩(0, σ_2^2 = 1),
$ <eq:1.5>

Where $λ_t$ is a Student's t-distribution with $ν$ degrees of freedom and scale parameter $σ_1^2$.

To appreciate the effect of the non-Gaussian noise, we investigate the probability density function of the Student's t-distribution with varying degrees of freedom, $ν$, and a normal distribution, which is shown in @fig:1.5_distributions.

#figure(
  image(
    "output/1_5_1_distributions.png",
  ),
  caption: [Probability density function of the Student's t-distribution with varying degrees of freedom $ν∈{1, 2, 5, 100}$ and a standard normal distribution.],
) <fig:1.5_distributions>

We note that as the degrees of freedom $ν$ increases, the Student's t-distribution approaches the normal distribution as expected.
For $ν=1$, the distribution is very heavy-tailed, which is expected to have a significant effect on the Kalman Filter. In effect, the underlying process is more likely to have noise realisations that are far from the mean, which can
lead to large jumps in the latent state vector $X_t$.

Given that the Kalman Filter models the latent state vector as having Gaussian noise, it will be unable to account for the large number of outliers that are likely to occur in low-$ν$ Student's t-distribution.
As a result, Kalman Filter will estimate an inflated estimate of the variance of the latent state vector.

With this in mind, the Kalman Filter will still usable in many circumstances,
thought with reduced performance. If the underlying process noise is known,
the filter could be modified to account for the non-Gaussian noise, though this is outside the scope of this assignment.

Additionally, we note that the confidence intervals calculated in @sec:1.3 are no longer valid, as they are based on the assumption of Gaussian noise.

Next, we can investigate the effect of the non-Gaussian noise on the process
described by @eq:1.5 by simulating realisations of the process with varying degrees of freedom $ν$ and comparing against the Gaussian process, which is the limit $ν=∞$.

#figure(
  image(
    "output/1_5_observations.png",
  ),
  caption: [Observations $Y_t$ of the processes given by @eq:1.5 with t-distributed process noise. Parameters are $a = 0.9, b = 1, σ_1^2 = 1, σ_2^2 = 1$ and $X_0 = 5$],
) <fig:1.5_observations>

The heavy-tailed nature of the t-distribution is clearly visible in @fig:1.5_observations, where we see large, sudden jumps in the observations $Y_t$.
This is particularly pronounced for $ν=1$, where the process is very heavy-tailed and the observations are erratic.

=== Parameter Estimation on Non-Gaussian Processes <sec:1.5.1>

We now perform the same parameter estimation as in @sec:1.4, but using the non-Gaussian process given by @eq:1.5 with $ν=∈{1, 2, 5, 100}$ and remaining parameters $a = 0.9, b = 1, σ_1^2 = 1$. We additionally perform the estimation using the Gaussian process $𝒩(0, σ_1^2)$ for comparison.
We simulate a total of $N=100$ realisations of each process, with each realisation consisting of $n=100$ observations.

The unmodified Kalman Filter is used to estimate the parameters through the
maximum likelihood framework as in @sec:1.4, with the results shown in @fig:1.5_experiment.

#figure(
  image(
    "output/1_5_experiment.png",
  ),
  caption: [Parameter estimates for $a, b, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $a = 0.9, b = 0.9, σ_1^2 = 1$],
) <fig:1.5_experiment>

As seen in @fig:1.5_experiment, we find reduced performance of the Kalman Filter for the very heavy-tailed
cases with $ν ∈ {1, 2}$, though the reduction in parameter estimation performance when compared against the Gaussian process is negligible for cases $ν≥5$. As before, we observe good estimation of the transition coefficient $a$ throughout all cases, but reduced accuracy on the bias term $b$ and process noise $σ_1$, as expected.
We observe that the estimate of the normal distribution variance, $σ_1^2$ is inflated for the heavy-tailed cases, which is consistent with our earlier discussion and observations of the heavy-tailed nature of the process.
Additionally, the variance of the estimate of $σ_1^2$ increases as the degrees of freedom $ν$ decreases.
Lastly, we find that the mean estimate of the bias term $b$ remains accurate, but the uncertainty in this
estimate increases rapidly as $ν$ decreases.

To aide in our understanding of the predictive performance of the Kalman Filter, we can also inspect the residuals
of the maximum likelihood estimation, which are shown in @fig:1.5_residuals.

#figure(
  image(
    "output/1_5_residuals.png",
  ),
  caption: [
    Distribution of the negative log-likelihood residuals from the parameter estimation
    of the processes given by @eq:1.5 with different process noise distributions.
  ]
) <fig:1.5_residuals>

This confirms our earlier observations, where we see that the mean residual of the maximum likelihood estimation
increases as the degrees of freedom $ν$ decreases, corresponding to a more heavy-tailed distribution.
Additionally, we observe a broadening of the distribution of the residuals, which
additionally would make parameter estimation based on a small number of realisations more challenging.
This is consistent with the broadening of the underlying distributions as shown in @fig:1.5_distributions.

In conclusion, we find that the Kalman Filter still performs well on non-Gaussian processes, though
for particularly heavy-tailed processes, the performance impairment may become unacceptable.
In most real-world scenarios, we would not expect the process noise to be this heavy-tailed,
which explains why the Kalman Filter finds such widespread use in industry and research.

While it is often reasonable to assume that the observation noise is Gaussian due to the Central Limit Theorem,
it may not be reasonable to assume the same of the process noise, $e_{1,t}$.
We have demonstrated, that in most cases non-Gaussian process noise can be handled adequately
by the Kalman Filter, though we note that the performance may be impaired for particularly heavy-tailed processes.

#pagebreak()
= Modelling a Transformer Station <sec:2>

In this section, we will be using simple state-space models (SSMs) to gain insights into the temperature of an electrical transformer station.
We are given a dataset with 168 observations as dependent variable $Y_t$ and 3 exogenous variables $T_(a \, t), Phi_(s \, t), Phi_(I \, t)$ which describe the outdoor temperature for the transformer in degrees Celsius, the horizontal global solar radiation at the station and the electrical load on the transformer.

== Exploratory Analysis <sec:2_1>

#figure(
  image(
    "output/2_1_exo_y_plot.png",
    width: 103%
  ),
  caption: [Given observation $Y_t$ and exogenous variables $T_(a \, t), Phi_(s \, t), Phi_(I \, t)$],
) <fig:2.1_exo_y_plot>

In @fig:2.1_exo_y_plot, there is very clearly a seasonality for day/night-cycles in all variables. The solar radiation $Phi_(s \, t)$ drops to $0$ during the night, which the load $Phi_(I \, t)$ mirrors almost exactly. It has a slightly quicker drop, once the sun is setting and during the peaks it
displays a wiggle, which suggests some sort of load controller or a load maximum with excess being discharged. Lower peaks or crumples in the radiation curve could be explained by cloud cover. A curious thing to notice, is that the load on the transformer $Phi_(I \, t)$ appears to have a quicker attack-time to rise, than the solar radiation $Phi_(s \, t)$.

Intuitively, we would expect the solar radiation to lead and load to lag. The $Y_t$ temperature follows $Phi_(s \, t)$ showing some cool-down period, once solar radiation dropped, hence a slower decay in temperature. Outdoor temperature $T_(a \, t)$ not only follows the solar radiation, hence daily 24h seasonality, but also exhibits a longer period seasonality, which could be climate and wheather effects.

/*
Overall, we can actually deduce a lot from just outdoor temperature and solar radiation cycles, especially the uninterrupted (unclouded) ones. When inspecting the graph, we can deduce about 17h of daylight, which excludes locations betwee $approx plus.minus 54$ degrees N/S.
In the southern-hemisphere there is only \'Tierra de Fuego\' the southern cape of Latin America that is still land-mass, but it does not match the temperature profile (as even in summer, for the long daylight hours, it has max. temperatures of about 8 degrees Celsius). One could possible match the outdoor temperature with weather data to deduce a more
accurate location.
*/


== 1D state-space model <sec:2_2_1D_SSM>

The goal here is to fit (estimate the parameters) of the following model via the Kalman-Filter and MLE:

$
  X_(t+1) = a X_t + B u_t + e_(1,t) \
  Y_t = c X_t + e_(2,t)
$<eq:2.2_1d_ssm>

with 
$
  u_t = [T_(a \, t) \, Phi_(s \, t) \, Phi_(I \, t)]^tack.b in bb(R)^(1 times 3); quad
  a in bb(R) \, B in bb(R)^(1 times 3) \, c in bb(R); quad
  e_(1 \, t) in bb(R) \, e_(2 \, t) in bb(R); quad
  X_t in bb(R)
$

For the fitting, the following contraints were set: The inital value for the hidden-state value was chosen at $X_0 = 20$ and $Sigma_(t+1 | t)^(x x)=1.0$, the parameters were initialized as $a=0.8 \, B=[0.05, 0.1, 0.1]^T \, c=1, sigma_1^2=log(2), sigma_2^2=log(2)$.
All entries of $a, B, c$ were constrained in the interval $[-2,2]$, while $sigma_1^2 \, sigma_2^2$, the variances of $e_(1 \, t) \, e_(2 \, t)$ were constrained in $[1e-3, log(10)]$.

These values were chosen somewhat arbitrarily (based on what worked well) within the boundry of the hints given in the exercise.

We fitted the model analogously to the framework introduced in @eq:1.3_Kalman_filter and @eq:1.4_implicit_Kalman_likelihood (for which the code can be found in the attached `2.ipynb` file). The Kalman-Filter provided predictions for the hidden-state, which where the basis for our negative log-likelihood to minimize.
The resulting estimated parameters rounded to the 4th decimal digit are:

$
  a = 0.7906; quad B=[0.1313 \, 0.0031 \, 0.2524]^T; quad c=0.8863; quad sigma_1 = 1e-3; quad sigma_2 = 1e-3
$ <eq:2.2_parameter_estimates>

It was notable, that the variances $sigma_1^2 \, sigma_2^2$ were always pushed to the lower boundary of the given contraint, no matter the initialization. This could imply, that the model can explain the observations $Y_t$ very well without added noise, or that the noise in the system is not of additive nature.

#figure(
  image(
    "output/2_2_predicted_yt_1d.png",
  ),
  caption: [given observation $Y_t$ compared to the output prediction of our @eq:2.2_1d_ssm model, based on the parameters @eq:2.2_parameter_estimates],
) <fig:2.2_predicted_observations>

@fig:2.2_predicted_observations shows that our simple model captures the dynamics of $Y_t$ relatively well. It does not perform well on the first $50$ hours, most likely because the data is more noisy for cloud-cover wheather conditions. Additionally, we can observe that $B_(1,3)$ is the highest coefficient in the @eq:2.2_1d_ssm model and is the factor to exogenous variable $Phi_(I,t)$, the load. As this variable shows the noisiest behaviour for that period, the 'culprit' is clear.

Looking at the residuals (in this case equivalent to the 'innovation') in @fig:2.2_residual_diagnostics, we can observe that they are approximately normally distributed (with small exceptions in the extreme value quantiles in the QQ plot). Hence, we do not diagnose a systematic error with the model.

#figure(
  image(
    "output/2_2_residual_diagnostics_1d.png",
  ),
  caption: [residual $hat(y)_t - Y_t$ of true temperature and the output prediction of our @eq:2.2_1d_ssm model, based on the parameters @eq:2.2_parameter_estimates],
) <fig:2.2_residual_diagnostics>

Beyond, we report the AIC and BIC as model selection criteria with $text("AIC")=495.15, text("BIC")=517.02$. In this setting, the BIC 'advantage' of penalizing model complexity heavier than AIC has not kicked-in yet, since with $p=7, n=168 arrow.r.double 2p > log(n)p$.

The physical implications and interpretations are extensively discussed in @sec:2_4_2d_state_interpretation for a 2D SSM. As the line of argument is analogous here for the 1D case, we will only briefly discuss the parameters.

The coefficient matrix $B$ can be seen as the weights for a sum-composition of $u_t$. Consequently, we can interpret the magnitude of the coefficients as importance weight. Apparently the last value of $B$ is the largest (@eq:2.2_parameter_estimates), which corresponds to $Phi_(t,I)$ being weighted as the most important predictor. Via the second coefficientin $B$, the importance of $Phi_(t,s)$ is weighted the lowest. Most likely, because from @fig:2.4_correlation_heatmap we can see that those two $Phi_t$ have a high correlation and thus share a lot of informational value for the model. Since one of these is already weighted high, there is no need to weigh the other one highly as well.
$T_(t,a)$ outdoor temperature is given about half as much importance as $Phi_(t,I)$ load.


== 2D state-space model <sec:2_3_2D_SSM>

The goal here is to fit (estimate the parameters) of the following extended 2D (meaning 2 hidden states) model via the Kalman-Filter and MLE:

$
  X_(t+1) = A X_t + B u_t + e_(1,t) \
  Y_t = C X_t + e_(2,t)
$<eq:2.3_2d_ssm>

with
$
  u_t = [T_(a \, t) \, Phi_(s \, t) \, Phi_(I \, t)]^tack.b in bb(R)^(1 times 3); quad
  A in bb(R)^(2 times 2) \, B in bb(R)^(2 times 3) \, C in bb(R)^(1 times 2); quad
  e_(1 \, t) in bb(R)^(2 times 1) \, e_(2 \, t) in bb(R); quad 
  bold(X)_t in bb(R)^(2 times 1)
$
both noise terms $e_t$ following (multivariate-)normal distributions with mean $0$.

For the fitting, the following inital values were chosen:
$
  X_0 = [20, 20]^T quad
  Sigma_(t+1 | t)^(x x)=10.0 \
  A=mat(
    0.8, 0.1;
    0.0, 0.7;
  ) quad
  B=mat(
    -0.1, 0.1, 0.1;
    0.0, -0.1, 0.0;
  ) quad
  c=mat(
    1.0, 0.2;
  )\
  sigma_1^2=log(2) quad sigma_2^2=log(2)
$<eq2.3_initial_values_2d>

subject to the following constraints: All entries of $A, B, C$ were constrained in the interval $[-2,2]$, while $sigma_1^2 \, sigma_2^2$, the variances of $e_(1 \, t) \, e_(2 \, t)$ were constrained in $[1e-3, 2]$.

These values were again chosen somewhat arbitrarily, within the ranges of the hint given in the assignment. It shall be noted, that we could achieve significantly better results, by fitting multiple times and taking the last estimated parameters as new initialisation for the next fitting.
However, we decided to only have one run with intentionally 'worse' initial values, to also get an impression of the performance of the fitting procedure itself.

After fitting, we have the following estimated parameters, rounded to the 4th decimal digit:

$
  A=mat(
    -0.8303, -0.3605;
    0.6865, 0.9543;
  ) quad
  B=mat(
    -1.7612, 1.741, -1.1408;
    0.9138, -0.4893, 0.8132;
  )\
  c=mat(
    0.264, 0.6254;
  ) quad
  sigma_1^2=0.001 quad sigma_2^2=0.001
$<eq:2.3_2d_parameter_estimates>

Again, the variances of the noise are pushed to the lower boundary.
We did initially to estimate the initial state $X_0$, however, even with constraints in the optimiser, the overall model yielded worse results. Yet the initial value estimates always hovered $approx [20, 20] = X_0$, therefore we deemed it a qualified guess.

#figure(
  image(
    "output/2_3_os_pred_2d_vs_1d.png",
  ),
  caption: [one-step predictions $hat(Y)_t$ of both models @eq:2.2_1d_ssm and @eq:2.3_2d_ssm],
) <fig:2.3_os_pred_2d_vs_1d>

In @fig:2.3_os_pred_2d_vs_1d we first look at the one-step predictions and compare the 2D to the 1D model. Visually, the performance seems approximately on par; both having significant difficulties capturing the dynamics for cloud-cover conditions and somewhat difficulties capturing peaks (both models over- and under-shooting at different peaks).

We report the information criteria statistics as $text("AIC")=499.23 \, text("BIC")=542.97$, which indicates that the 2D model does not perform significantly better. Because it is still the case, that AIC penalizes model complexity more heavily here (with $p=14$), we can even argue that the 2D model performs slightly worse than the 1D model, compared with its complexity.

Thus, we turn to the residual analysis in @fig:2.3_residual_diagnostics_2d. This follows the above assumptions, that the model does not perform significantly worse or better. We still cannot diagnose a systematic error of the 2D SSM, as the residuals appear approximately normally distributed.
Despite, one thing to note is that from the top plot we can see a weak trend, that the model tends to change from under-shooting to over-shooting the true $Y_t$.
This is confirmed in the ACF, as $rho$ is not alternating as strongly as in the 1D model. The PACF suggests that there is not enough information to diagnose a systematic trend in residuals (yet).

#figure(
  image(
    "output/2_3_residual_diagnostics_2d.png",
  ),
  caption: [residual $hat(y)_t - Y_t$ of true temperature and the output prediction of our @eq:2.3_2d_ssm model, based on the parameters @eq:2.3_2d_parameter_estimates],
) <fig:2.3_residual_diagnostics_2d>

The lack of performance difference between the 1D hidden-state model and the 2D hidden-state model, may suggest, that there is either no significant value addition in a 2nd hidden-state, hence, $X_(t,0)$ and $X_(t,1)$ would strongly correlate. There could also be a notion of inverse or counter-acting relationship between the two hidden-states, that nulls out or corrects any information gain of the additional second state compared to only one state.
We explore this further in the next section.


== 2D state interpretation & discussion <sec:2_4_2d_state_interpretation>

The goal of this section is to inspect the estimated hidden-states of the 2D model @eq:2.3_2d_ssm and discuss insights, that can be derived from the coefficients $A, B, C, sigma_1, sigma_2$ or the state time-series $bold(X_t)$ themselves.

#figure(
  image(
    "output/2_4_hidden_states_exo.png",
  ),
  caption: [estimated hidden state $X_(t,0)$ (because visually closer to exog.) and exogenous variables $u_t$; standardized for better comparison],
) <fig:2.4_hidden_states_exo>

We start with a visual analysis of the hidden states $bold(X_t)$ in @fig:2.4_hidden_states_exo. We can see that both states are seasonally counter-acting, as assumed in @sec:2_3_2D_SSM. In fact, the Pearson correlation coefficient between $X_(t,0), X_(t,1)$ is $-0.9961$, so almost perfect negative correlation. Thus confirms the assumption, that there is not much value added by a second hidden state.
By construction, the hidden states are leading the observations in signal response (in upper graph). The value ranges for $X_t$ are not quite insightful, as they can easily be absorbed by the magnitude of coefficients in $A$. Yet, the opposing nature suggests that one state is 'buffering' the other.
In the bottom graph of @fig:2.4_hidden_states_exo, we only plot one state, as already deduced that there is not much additional informational value in the second state. 
Overall, we know that, in principle, the hidden state is just a weighted sum (because it is a linear combination) of the exogenous variables. We can see that $X_(t,0)$ is most closely followed by $Phi_(t,s)$ (which can be seen in @fig:2.4_correlation_heatmap), which suggests the highest weight/coefficient assigned. 

#let fig = [#figure(
  image(
    "output/2_4_correlation_heatmap.png",
    width: 100%
  ),
  caption: [correlation coefficient heatmap],
) <fig:2.4_correlation_heatmap>]

#wrap-content(fig, [
  The heatmap @fig:2.4_correlation_heatmap also shows, that even the exogenous variables in $u_t$ have high statistical correlation among themselves. This further supports the conclusion, that estimated parameters will favour one variable, as the others do not add much information to the model.

  We follow with an analysis of the estimated coefficients: As mentioned above, we can interpret the coefficients in $B$ as weights for the sum-composition of $u_t$ for each state. The sign between the first and second row of $B$ are flipped (c.f. @eq:2.3_2d_parameter_estimates), explaining the opposing behaviour and the 'buffering/dampening' nature. We can see in the top graph of @fig:2.4_stepwise_coeff_states, that this weighted sum is somewhere in between the total sum and the average of all exogenous variables in $u_t$ combined. 
])
  
The highest weights are both in column 1 of $B$, the weights for $T_(t,a)$. This is likely due to the extremely different magnitude of the exogenous variables, where $T_(t,a)$ has the smalles values, so it needs to be boosted to even compete with the other exogenous variables in the weighted sum. Additionally, $T_(t,a)$ has the lowest correlation (c.f. @fig:2.4_correlation_heatmap) with any of the other exogenous variables (or states), so it provides the most uncovered information.

$A$ acts as the recursive state factor, also allowing for further dampening of the high magnitude of the states $bold(X)_t$. This matrix parameter is responsible for convergence and stability of the entire SSM (c.f. @sec:2.4_1_stability).
The parameter matrix $C$ is the factor that combines both states into a prediction for the observation $Y_t$. It is therefore import, as it controlls the mixing of states. 

Interestingly, because of the high correlations of variables in the entire SSM, in the bottom plot of @fig:2.4_stepwise_coeff_states, we can observe that skipping the recursive estimation of a hidden-state entirely (by modelling $C (B u_t)$), we still get a decent approximation of our observations $Y_t$. This would suggest and confirm our previous assumptions, that $Y_t$ could possible be predicted directly by a linear model of only $u_t$.

To better illustrate the effects of the parameters, we visualise some time-series with the estimated states $bold(X_t)$ and exogenous variables $u_t$ multiplied by the parameters sequentially.

#figure(
  image(
    "output/2_4_stepwise_coeff_states.png",
  ),
  caption: [sequential application of estimated parameters to the states $bold(X)_t$],
) <fig:2.4_stepwise_coeff_states>

The interpretation of @fig:2.4_stepwise_coeff_states is mostly conducted above. A not yet mentioned interesting insight is that the combination of $A bold(X)_t$ switches sign in the middle of the time axis. This coincides with the slight trend notices in the residuals in @fig:2.2_residual_diagnostics. So after all, it might be that improving the estimates of $A$ could avoid the over- and under-shooting trend.

A physical interpretation of the states is, as mentioned, a combination of the external factors at the measurment station, $u_t$ as predictor, mostly (c.f. @fig:2.4_correlation_heatmap) $Phi_(s,t)$.

One could interpret one state $X_(t,1)$ as a buffer for $X_(t,0)$ (c.f. @fig:2.4_hidden_states_exo). Based on the above conclusions, mainly buffering/adjusting $X_(t,0)$ for the outdoor temperature $T_(t,a)$. Alas, within that line of argument, one state would represent a "cooling" and the other a "solar-radiation-load-temperature" response.

Yet, as already argued, this would not be necessariy as one single state carries almost just as information, hence this buffering effect for $T_(t,a)$ is most likely absorbed into $A$ (not into the noise, as the noise terms both $e_(t,1), e_(t,2) arrow.r 0$ during the MLE).

Overall, the physical interpretation of the model makes sense. Nevertheless, it also became clear that the model is too complex for the informational value it contains/processes.


=== Stability <sec:2.4_1_stability>

While not asked for, this section seems important after we encountered instability in the first version of the assignment in part 1.
Very briefly: For stability of an SSM we focus only on the homogenous part of the system, namely:

$
  X_(t+1) = A X_t + ...\
  Y_t = C X_t + ...
$

Assuming the exogenous part $u_t$ is in itself is bounded and stable, it does not affect the stability of the system, it only gives a 'direction' of movement. We also assume the variance of the state equation $sigma_1^2$ to be reasonably bounded.
More specifically, within the homogenous part of the system, we look at matrix $A$. If the spectrum of $A$ is contractive, in other words, if $forall lambda_i in bold(lambda)(A): |lambda_i|<1$, then the SSM is assumed to be asymptotically stable (for full details, see @Anderson_1979 chapter 4)


=== Outlook <sec:2.4_2_outlook>

Having analysed two simpler SSMs, we can give an outlook on what to improve. As we have seen in this example, increasing the number of states does not significantly improve the performance of the model. Thus, we must think of other routes to improve performance.

One option is to employ Kalman-smootheners, which could help tackle the large prediction error in especially the first 50h of the expirement.

Another option would be to include an auto-regressive term in the observation equation of the model. This could assist in better performance around the peaks, as the hidden-states seem to pre-maturely drop after peeks, in some occasions.

Furthermore, it could be beneficial to include another parameter-exogenous-term ($Y_t = C X_t + D u_t + e$) in the observation equation. As we saw in @fig:2.1_exo_y_plot, the two exogenous variables $Phi_(s,t), Phi_(I,t)$ seemed to be leading in dynamics and seasonality. Re-enforcing the exogenous $u_t$ onto the observation part, could help capture that behaviour better.

#bibliography("report.bib")
