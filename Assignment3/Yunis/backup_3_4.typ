== Impulse Response<sec:3_4_impulse_response>

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
 & #h(0em) #h(0em) dots.v\
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
  image("../output/3_4_impulse_response_Tt.png"),
  caption: [Impulse Response of the AR(10)-X(1,1) with the unit impulse from $T_d$]
) <fig:3_4_ir_Tt>

If the phrasing "Present it for both variables" is interpreted as presenting the impulse response with a unit impulse from both variables separately, we can give a second plot of just an impulse from $G_v$:

#figure(
  image("../output/3_4_impulse_response_Gv.png"),
  caption: [Impulse Response of the AR(10)-X(1,1) with the unit impulse from $G_v$]
) <fig:3_3_ir_Gv>

As already explained based on @fig:3_1_analysis, the two impulse responses counter-act. The IR (impulse response) for this model decays rather quickly, especially with the impact of an exogenous variable that is reciprocal to the target variable. Most of these effects were already assumed based on available information.

We can interpret, that the model with order $10$ is rather robust to unit shocks from both exogenous variables. As a consequence, with this model, we can accept a certain degree of fluctuation in the exogenous variables, without having to worry about a big impact on the predictive performance of our model. For example, the effect of a measurements error or sensor defect in $T_d$ or $G_v$ that acts as a shock to our system, will not last for very long.
Thus, having a larger order for the AR part of our model pretects against shocks from the other variables. On the other hand, the model is more vulnerable to sudden changes in auto-regressive behaviour.
