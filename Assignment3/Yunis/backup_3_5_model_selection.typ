== Model Selection<sec:3_5_to_3_8_model_selection>

In this section, we will stick to the notation convention for ARX models introduced above and create a selection of models for comparison.

=== Linear Regression <sec:3_5_ar0_x1_1>

For the simple linear regression (OLS) model we have:

$
  P_t = omega_1 T_t + beta_1 G_t + epsilon_t
$<eq:3_5_ols_model>

with the design matrix: $X = [T_t \, G_t]$ for the column-vectors
and with $Theta = [omega_1 \, beta_1]^T$ as parameter vector. The estimation (fitting the model) for $hat(Theta)$ will be analogous to @eq:3_4_ols_parameter_est.

The resulting parameters are $hat(Theta) = [3.8948 \, -0.1099]^T$ with a corresponding RMSE of $approx 5.407$. Thus we have our linear regression predictions as $X dot.op hat(Theta) = hat(P_t)$.
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

For a fair comparison, a burn-in period of 2 observations was granted to the one-step predictions, to deal with reasonable values. Still, a RMSE of $22.6638$ was noted.

#figure(
  image("output/3_5_residual_analysis_ols_os.png"),
  caption: [model and residual analysis or linear model and its one-step predictions]
) <fig:3_5_residual_analysis_ols_os>

In @fig:3_5_residual_analysis_ols_os, we can observe that the one-step prediction model performs slighty worse than the full linear model. Especially in the sharp drops, the one-step approach does not quite capture the movement as fast. The ACF and CCF plots show, that the residuals are far from being white-noise (which is the desirable state; then the error is $epsilon_t$). This indicates a systematic error in the predictions, so there is room for improvement of the model.

TODO
- need for a transfer function?

=== AR(1)-X(1,1) Model<sec:3_6_ar1_x1_1>

Now, we deal with a model that supplements the previous model by one AR component (but already modeled in @eq:3_4_ar1_x1_1):

$
  P_t = -phi.alt_1 P_(t-1) + omega_1 T_t + beta_1 G_t + epsilon_t
$<eq:3_6_ar1_x1_1>

With the OLS estimated parameters $hat(Theta) = [-hat(phi.alt_1) \, hat(omega_1) \, hat(beta_1)]^T = [-0.4127 \, 2.3352 \, -0.0836]^T$ we get a RMSE of $0.353$, which is a significant improvement to before, just by adding an AR component to the model. Since the ACF and CCF only confirmed, what we could already observe in the residual plots themselves, we omit them this time:

#figure(
  image("output/3_6_ar1_x1_1_ols_os.png"),
  caption: [AR(1)-X(1,1) model and residual analysis with OLS estimates and its one-step predictions]
) <fig:3_6_ar1_x1_1_ols_os>

We granted a burn-in period of 4 observations to the one-step estimates for a fairer comparison.
From @fig:3_6_ar1_x1_1_ols_os we can conclude that alas both one-step and full OLS predictions produce significantly smaller residuals, the distribution of the residuals fails in the same ways as before. There is still a a visible seasonality to the magnitude of prediction errors. This can already be seen in the bare plot of the time-series: The steep drops and consecutive ascends are not captured correctly. While indeed a better model already, still a systematic error.


=== AR(2)-X(2,2) Model<sec:3_7_ar2_x1_1>

This time we add one order to the AR and both exogenous parts:

$
  P_t = -phi.alt_1 P_(t-1) - phi.alt_2 P_(t-2) + omega_1 T_t + omega_2 T_(t-1) + beta_1 G_t + beta_2 G_(t-1) + epsilon_t
$<eq:3_7_ar2_x2_2>

With the OLS estimated parameters $hat(Theta) = [-hat(phi.alt_1) \, -hat(phi.alt_2) \, hat(omega_1) \, hat(omega_2) \, hat(beta_1) \, hat(beta_2)]^T = [-0.6274 \, 0.055 \, 0.2704 \, 1.4412 \, -0.0991 \, 0.037]^T$ we get a RMSE of $0.2842$, which is again an improvement to the AR(1)-X(1,1) model.

#figure(
  image("output/3_7_ar2_x2_2_ols_os.png"),
  caption: [AR(2)-X(2,2) model and residual analysis with OLS estimates and its one-step predictions]
) <fig:3_7_ar2_x2_2_ols_os>

While this model fits a lot tighter, the notion of not capturing the drops and steeper ascends persists. One can beautifully see this behaviour around $upright("t_hour")=75$, where there is a jump downwards in the residual and then something like a staircase upwards until the model catched up again. There are several of the 'staircases' in the residuals.

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

Here we only consider full OLS estimates and not the one-step predictions (as those did worse in every model and need a burn-in period).

#figure(
  image("output/3_7_aic_bic_model_comparison.png"),
  caption: [AIC, BIC comparison of different models orders from @sec:3_5_ar0_x1_1, @sec:3_6_ar1_x1_1, @sec:3_7_ar2_x1_1]
) <fig:3_7_aic_bic_model_comparison>

This is a classic Elbow curve. In various modelling fields, information theoretical measures or straight residual measures (such as MSE, etc.) are plotted against different model variations. The general idea is, to deduce at which point, there is a good trade-off between parameters of the model and the model performance.
In this case, we have the x-axis as increasing model order, hence increasing complexity and proneness to overfitting. Hence, we are looking for a trade-off between performance and model complexity. It is an interesting choice of metric to do this kind of plot. AIC and BIC inherently punish model complexity, as long as $log(n) > 2$, the BIC does even more so; this also explains the slight difference between those two curves. Thus they theoretically already account for one of the decision parameter that such Elbow-curve is intended to help with.
If we look at the plot, there is a clear trend towards the more complex AR(2)-X(2,2) model. Traditionally, the decision here would be to select the AR(1)-X(1,1) model, as it already performs relatively well and the performance gain to the more complex model is very small. Hence, favour simplicity before performance, in the hope of better generalisation.

However, as the AIC and BIC both already punish the increased model complexity of AR(2)-X(2,2) model, but still produce a lower score, the interpretation would yield in: "more complex, but worth it". Therefore, the model selection would fall onto the AR(2)-X(2,2) model.

=== One-Step Predictions & RMSE (3.8)<sec:3_8_OSPred_RMSE>

When producing one-step predictions, a critical factor is to mind the burn-in period. As the model parameters are re-fitted at each step of the prediction (refer to @eq:3_5_os_pred). Taking into account the new predicted values of the series, in the first few steps, there is little data to fit on. Thus, the predictions are not very strong. This is what we refer to as the burn-in period.

As we must give one-step predictions on the test-dataset (which is already quite few observations), it makes sense to include some of the last observations of the train-dataset, to counter the burn-in.
This way, there is hopefully not much of an impact on the predictions of the test-dataset.
Additionally, this would also be a more fair comparison to the previous exercises, as those are OLS fitted on the entire train-dataset, thus do not suffer from a burn-in, as do one-step predictions.

To test out, how much of a spill from the train-dataset we need, we can check fitting the one-step prediction on *only* the test-data.

#figure(
  image("output/3_8_rmse_model_comparison_only_testdata.png"),
  caption: [RMSE comparison of different models orders from @sec:3_5_ar0_x1_1, @sec:3_6_ar1_x1_1, @sec:3_7_ar2_x1_1 on only the test dataset]
) <fig:3_8_rmse_model_comparison_only_testdata>

The decision rule on how much burn-in we grant the model is: Increase the burn-in $b$ (meaning, exclude the first $b$ observations and predictions) until there are no more negative values in the one-step predictions of $P_h$. This is a context based rule, as negative wattage for the heater, the variable $P_h$ represents, does not make sense.
The @fig:3_8_rmse_model_comparison_only_testdata shows the results for $b=11$. Therefore, we include the last 11 samples/observations of the train-dataset into the design matrix $X$ for the one-step predictions.

Now including 11 samples from the train-data as spill over, we observe the following:

#figure(
  image("output/3_8_rmse_model_comparison_train_test_spill.png"),
  caption: [RMSE comparison of different models orders from @sec:3_5_ar0_x1_1, @sec:3_6_ar1_x1_1, @sec:3_7_ar2_x1_1 on concatenation of train and test dataset]
) <fig:3_8_rmse_model_comparison_train_test_spill>

On the new concatenated series, we really only need $b ≥ 6$ to not have absurd values in $P_h$ (negative or >1000), but with $b=8$, we now effectively have 64 one-step predictions as asked for in the task.

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

REVISE
- therefore the RMSE formula given in the task is actually incorrect, as it averages over 3 more values than there are one-step predictions available (64 test-data samples; 3 casualties, 2 for model order and 1 for OS-predictions)
- that is not even considering a burn-in!
- hence the effective number of observations fair for comparison to the previous models is much lower than 64
- please, fucking think about what you write Peder, at least 64 fucking times, before you confuse the shit out of people for hours!

Reflecting upon the question: Does this yield the same model selection as via the AIC, BIC metrics / criteria?

Yes and no. Only going via the RMSE, we could conclude from the first approach, where we only used pure test-dataset samples, that we should select the AR(2)-X(2,2) model.
However, after adjusting the setting (to have a spill of samples from the train-dataset) for a technically more fair comparison, the RMSE would yield the conclusion to select the AR(1)-X(1,1).
However, to put this into perspective, the RMSE is overall much lower with the second modelling approach. This makes sense makes sense, since there is much more information available to the model. This lets us conclude that there is a lot of relevant information in the further past.
Yet, there seems to be a notion that a model with shorter lag, an AR(1) component, captures this better than a longer lag. This could be interpreted as that past information is overall valuable, but the series behaviour is ultimately very short-term oriented.


=== k-Step Predictions<sec:3_9_kstep_pred>

We selected the AR(2)-X(2,2) model, mainly because of the AIC, BIC selection criteria. From experience, they produce a very reliable selection process.
We provide a selection of $k$ step-widths, which are usually choices that would make sense in real-life settings: 12h, 24h, 48h.

#figure(
  image("output/3_9_kstep_predictions_selection.png"),
  caption: [k-step predictions for intervals]
) <fig:3_9_kstep_predictions_selectionl>

Beyond that we look at a plot of the RMSE against the step-width k with the intend to draw some conclusions about the error behaviour over longer periods of prediction. As we already pointed out problems with the residual distribution, this seems reasonable. We also acknowledge, that the allowed burn-in period plays a crucial roll in fitting k-step models, thus we compare them as well.

#figure(
  image("output/3_9_rmse_kstep_pred.png"),
  caption: [k-step predictions for intervals]
) <fig:3_9_rmse_kstep_pred>

As expected, the RMSE reacts heavily to the allowed burn-in. Taking the approach of also allowing a spill of data from the train-dataset, we can afford to accept a longer burn-in period, as we can analogously increase the spillage of data from train-dataset.
From $b≥9$ the RMSE stays mostly below $5$, which is on most realistic cases acceptable for this problem set-up.
Nonetheless, there are these devious peaks that display some sort of seasonality for the RMSE over the step-width. How can we explain these?

$P_h$ clearly shows seasonal drops in power (@fig:3_1_analysis). In some intervals however, there are exceptions, where the drop is not as sudden or not as steep.

REWRITE
- there are a few exceptions to these drops:
        - between 2013-02-01 and 2013-02-02
        - between 2013-02-04 and 2013-02-05
        - between 2013-02-07 and 2013-02-08
    - the drops are not as sharp
    - probably the model parameters, when fitted, encode an expectation of seasonality at these exception points as well
    - as this doesn't happen the model predicts wrong, thus the error peaks
    - a 2nd explanation would be the weakness of the predictions around the double (inverse) peaks
        - between 2013-02-02 and 2013-02-03
        - between 2013-02-06 and 2013-02-07

Overall, the model performs well, the deviations are small, however, the residuals do not look like white noise yet. There is still a visible auto-regressive behaviour present. This could mean, that the model order or structure is not yet sufficient, a more complex model could resolve that.
Usually, the last option is to increase model complexity, but another argument for that is: The distribution of residuals is patterned. As indicated earlier, the model systematically overshoots for values of $P_h$ in the 'buckle' right after the steep drop of values in @fig:3_7_ar2_x2_2_ols_os and systematically undershoots values for the bottom of the drops.

Apart from the prediction accuracy (which could potentially be improved), in an operational setting the model would work just fine with multi-step predictions. It is not computationally expensive and would be easily integrated into a prediction pipeline.
Naturally, the further into the future we predict, the higher the uncertainty; that notion is present without calculating the prediction intervals. Beyond that, we saw certain points of the series (towards the drops of $P_h$) where the predictions tend to deviate; FIGURE, the seasonality of the RMSE over $k$ step-width.
This does not improve if we simply choose a $k$ step-width, that produces a good RMSE, since we would work on a continually re-fitting real-time prediction pipeline. Hence, by the nature of the series, this point of drop in $P_h$ will naturally come, regardless of the k-step.

Modelling-wise an RLS model with higher "forgetting" coefficient may be a viable option, to account for these know effects.

Given that knowledge on the historic data, it would be wiser to simply accept an increased uncertainty in specific periods of the series and contextualize consequences: Reduce heating at the drops, but also have power available in case of slight mis-predictions.

It could also be an option to reduce the measurement interval from hourly to 5-min intervals, to be able to react quicker and counter prediction uncertainty.

== Conclusions<sec:3_10_conclusions>

*Some Reflections on the model selection process (@sec:3_6_ar1_x1_1 to @sec:3_8_OSPred_RMSE):*

The construction of the selection process via the exercise objectives, makes the selection process not a fair, comparable process, for a number of reasons:

1. The AIC, BIC are calculated on the training dataset, while the RMSE (without our adjustment) is calculated on the test dataset. This would make for a good approach, if the same metric was used, to verify how well the model behaves on different time-slices, basically how well the model generalises. Yet, this is not the case!
2. two different prediction approaces are used (full OLS fit and step-wise)
3. The burn-in for one-step or multi-step predictions is not considered.
4. Too few options for models are considered, the residuals still show auto-regressive behaviour.

We are comparing 2 different metrics (AIC, BIC together as information criteria) on 2 different time slices with 2 different prediction methods. If anything, it is not surprising to get different conclusions.
Given the drawbacks, this is not a very solid model selection process, comparing apples with pears.

#bibliography("report.bib")
