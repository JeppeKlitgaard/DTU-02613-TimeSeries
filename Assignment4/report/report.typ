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
$ <eq:1.1>

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

The 5 independent realisations of the state as given by @eq:1.1 can be seen in @fig:1.1. Note that
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
    "output/1_4_experiment_a_1_b_0.9_sigma_1_1.png",
    width: 50%
  ),
  caption: [Parameter estimates for $A, B, σ_1^2$ using the Kalman Filter and maximum likelihood estimation with $A = 1, B = 0.9, σ_1^2 = 1$],
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

We find that for all three rounds and across all three parameters, we are able to estimate the parameters with a reasonable degree of accuracy, though we note that
the distribution of the estimates is quite large and the true value is only well-approximated when averaging across multiple realisations.

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
  caption: [Observations $Y_t$ of the processes given by @eq:1.5 with t-distributed process noise. Parameters are $A = 0.9, B = 1, σ_1^2 = 1, σ_2^2 = 1$ and $X_0 = 5$],
) <fig:1.5_observations>

= Modelling a Transformer Station <sec:2>

In this section, we will be using simple state-space models (SSMs) to gain insights into the temperature of an electrical transformer station. 
We are given a dataset with 168 observations as dependent variable $Y_t$ and 3 exogenous variables $T_(a \, t), Phi_(s \, t), Phi_(I \, t)$ which describe the outdoor temperature for the transformer in degrees Celsius, the horizontal global solar radiation at the station and the electrical load on the transformer.

== Exploratory Analysis <sec:2_1>

#figure(
  image(
    "output/2_1_exo_y_plot.png",
  ),
  caption: [Given observation $Y_t$ and exogenous variables $T_(a \, t), Phi_(s \, t), Phi_(I \, t)$],
) <fig:2.1_exo_y_plot>

In @fig:2.1_exo_y_plot, there is very clearly a seasonality for day/night-cycles in all variables. The solar radiation $Phi_(s \, t)$ drops to $0$ during the night, which the load $Phi_(I \, t)$ mirrors almost exactly. It has a slightly quicker drop, once the sun is setting and during the peaks it
displays a wiggle, which suggests some sort of load controller or a load maximum with excess being discharged. Lower peaks or crumples in the radiation curve could be explained by cloud cover. A curious thing to notice, is that the load on the transformer $Phi_(I \, t)$ appears to have a quicker attack-time to rise, than the solar radiation $Phi_(s \, t)$. 

Intuitively, we would expect the solar radiation to lead and load to lag. The $Y_t$ temperature follows $Phi_(s \, t)$ showing some cool-down period, once solar radiation dropped, hence a slower decay in temperature. Outdoor temperature $T_(a \, t)$ not only follows the solar radiation, hence daily 24h seasonality, but also exhibits a longer period seasonality, which could be climate and wheather effects.

Overall, we can actually deduce a lot from just outdoor temperature and solar radiation cycles, especially the uninterrupted (unclouded) ones. When inspecting the graph, we can deduce about 17h of daylight, which excludes locations betwee $approx plus.minus 54$ degrees N/S. 
In the southern-hemisphere there is only \'Tierra de Fuego\' the southern cape of Latin America that is still land-mass, but it does not match the temperature profile (as even in summer, for the long daylight hours, it has max. temperatures of about 8 degrees Celsius). One could possible match the outdoor temperature with weather data to deduce a more
accurate location.


== 1D state-space model <sec:2_2_1D_SSM>

The goal here is to fit (estimate the parameters) of the following model via the Kalman-Filter and MLE:

$
  X\_{t+1} = a X\_t + B u\_t + e\_{1,t} \\\\
  Y\_t = c X\_t + e\_{2,t}
$ <eq:2.2_1d_ssm>

with:
- $u_t = [T_(a \, t) \, Phi_(s \, t) \, Phi_(I \, t)]^tack.b in bb(R)^(1 times 3)$
- $a in bb(R) \, B in bb(R)^(1 times 3) \, c in bb(R)$
- $e_(1 \, t) in bb(R) \, e_(2 \, t) in bb(R)$
- $arrow.r.double X_t in bb(R)$

For the fitting, the following contraints were set: The inital value for the hidden-state value was chosen at $X_0 = 20$ and $Sigma_(t+1 | t)^(x x)=1.0$, the parameters were initialized as $a=0.8 \, B=[0.05, 0.1, 0.1]^T \, c=1, sigma_1^2=log(2), sigma_2^2=log(2)$.
All entries of $a, B, c$ were constrained in the interval $[-2,2]$, while $sigma_1^2 \, sigma_2^2$, the variances of $e_(1 \, t) \, e_(2 \, t)$ were constrained in $[1e-3, log(10)]$.

These values were chosen somewhat arbitrarily (based on what worked well) within the boundry of the hints given in the exercise.

We fitted the model analogously to the framework introduced in @eq:1.3_Kalman_filter and @eq:1.4_implicit_Kalman_likelihood (for which the code can be found in the attached `2.ipynb` file). The Kalman-Filter provided predictions for the hidden-state, which where the basis for our negative log-likelihood to minimize.
The resulting estimated parameters rounded to the 4th decimal digit are:

$$
  a = 0.7906 \, B=[0.1313 \, 0.0031 \, 0.2524]^T \, c=0.8863 \, sigma_1 = 1e-3 \, sigma_2 = 1e-3
$$ <eq:2.2_parameter_estimates>

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

Beyond, we report the AIC and BIC as model selection criteria with $text("AIC")=495.15, text("BIC")=517.02$. In this setting, the BIC 'advantage' of penalizing model complexity heavier than AIC already kicks-in, since with $p=7, n=168 arrow.r.double 2p < log(n)p$.


– Are there periods where the model performs poorly (e.g., daytime peaks)?
– What do the estimated parameters tell you about the influence of load, solar radiation,
and outdoor temperature?
• Hint: Consider the physical meaning of the parameters (e.g., effect of load on temperature
rise).

#bibliography("report.bib")
