import numpy as np
import matplotlib.pyplot as plt
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
from statsmodels.graphics.gofplots import qqplot
from matplotlib.gridspec import GridSpec

def compute_aic_bic(log_lik, num_params, num_obs):
    """Compute AIC and BIC given log-likelihood, number of parameters, and observations."""
    aic = 2 * num_params - 2 * log_lik
    bic = np.log(num_obs) * num_params - 2 * log_lik
    return aic, bic

def plot_residual_diagnostics(residuals, series_name:str="Residuals"):
    fig = plt.figure(figsize=(16, 6), dpi=400)
    fig.suptitle("Residual $\hat{y}_t-Y_t$ Diagnostics")
    gs = GridSpec(2, 3, figure=fig)

    ax1 = fig.add_subplot(gs[0, :2])
    plot_acf(residuals, ax=ax1, lags=40, title=f"{series_name} ACF")
    ax1.set_xlabel("$k$ lags")
    ax1.set_ylabel("$\\rho$ auto-corr.")
    ax1.set_ylim((-1.2, 1.2))

    ax2 = fig.add_subplot(gs[1, :2])
    plot_pacf(residuals, ax=ax2, lags=40, title=f"{series_name} PACF", method='ywm')
    ax2.set_xlabel("$k$ lags")
    ax2.set_ylabel("$\\rho$ auto-corr.")
    ax2.set_ylim((-1.2, 1.2))

    ax3 = fig.add_subplot(gs[:, 2])
    qqplot(residuals, line='s', ax=ax3)
    ax3.set_title(f"{series_name} QQ-Plot")

    plt.tight_layout()
    plt.show()

def plot_state_space():
    pass

def plot_obervations(y_pred, df):
    plt.figure(figsize=(10,4), dpi=300)
    plt.title("Transformer Station Temperature")
    plt.plot(y_pred, label="$\hat{Y}_t$")
    plt.plot(df["Y"], label="$Y_t$")
    plt.xlabel("$t$ time (hours)")
    plt.ylabel("$Y_t$ in ◦C")
    plt.legend()