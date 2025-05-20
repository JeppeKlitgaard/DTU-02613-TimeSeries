import pandas as pd
import numpy as np
from scipy.optimize import minimize

# for single exogenous observation ssm / AR(0);
n, p, m = 100, 3, 2 
A = np.random.randn(m,m)
B = np.random.randn(m,p)
u = np.random.randn(n,p)*13     # p exogenous variables
C = np.random.randn(1,m)
D = np.zeros((1,p))
Q = np.eye(m)*0.8
R = np.eye(1)*0.2

def simulate_ssm(A:np.ndarray=A, B:np.ndarray=B, C:np.ndarray=C, D:np.ndarray=D, Q:np.ndarray=Q, R:np.ndarray=R, u:np.ndarray=u, x0=None, seed:int=42):
    """
    Simulate a linear Gaussian state-space model.
    
    Parameters:
        A, B, C, D: system matrices
        Q, R: process and observation noise covariances
        u: exogenous input of shape (T, p)
        x0: initial state (n,), defaults to 0
        seed: random seed for reproducibility
        
    Returns:
        y: simulated observations (T, m)
        x: latent states (T+1, n)
    """
    rng = np.random.default_rng(seed)
    T, p = u.shape
    n = A.shape[0]
    m = C.shape[0]
    
    x = np.zeros((T+1, n))
    y = np.zeros((T, m))
    x[0] = np.ones(n)
    if x0 is not None:
        x[0] = x0

    for t in range(T):
        w_t = rng.multivariate_normal(np.zeros(n), Q)
        v_t = rng.multivariate_normal(np.zeros(m), R)
        x[t+1] = A @ x[t] + B @ u[t] + w_t
        y[t] = C @ x[t+1] + D @ u[t] + v_t

    return y, x

def predict_ssm(A, B, C, D, Q, R, u_future, xT, PT, seed=None, return_samples=False):
    """
    Predict future observations from a linear Gaussian SSM.
    
    Parameters:
        A, B, C, D: system matrices
        Q, R: process and observation noise covariances
        u_future: future exogenous input (H, p)
        xT: current latent state mean (n,)
        PT: current latent state covariance (n, n)
        seed: optional random seed
        return_samples: whether to sample from predictive distribution
        
    Returns:
        y_mean: predictive means (H, m)
        y_std: predictive std deviations (H, m)
        samples (if return_samples=True): (H, m)
    """
    rng = np.random.default_rng(seed)
    H = u_future.shape[0]
    n = A.shape[0]
    m = C.shape[0]
    
    x_pred = xT
    P_pred = PT
    y_mean = np.zeros((H, m))
    y_std = np.zeros((H, m))
    y_samples = np.zeros((H, m)) if return_samples else None

    for t in range(H):
        # Predict next latent state
        x_pred = A @ x_pred + B @ u_future[t]
        P_pred = A @ P_pred @ A.T + Q
        
        # Predict observation
        y_mean[t] = C @ x_pred + D @ u_future[t]
        S = C @ P_pred @ C.T + R
        y_std[t] = np.sqrt(np.diag(S))
        
        if return_samples:
            y_samples[t] = rng.multivariate_normal(y_mean[t], S)
    
    if return_samples:
        return y_mean, y_std, y_samples
    else:
        return y_mean, y_std
    

# SSM - generic State Space Model
def kalman_loglik(Y:np.ndarray, A:np.ndarray=A, B:np.ndarray=B, C:np.ndarray=C, D:np.ndarray=D, Q:np.ndarray=Q, R:np.ndarray=R, u_df:pd.DataFrame=None, x0:int=20):
    U = u_df[["Ta", "S", "I"]].values
    T = len(Y)
    n = A.shape[0]

    # initialization
    x_pred = np.ones(n)*x0
    P_pred = np.eye(n)*20.0

    loglik = 0.0
    for t in range(T):

        # Prediction step
        x_pred = A @ x_pred + B @ U[t]
        P_pred = A @ P_pred @ A.T + Q**2

        # Observation prediction with D @ u_t
        y_pred = (C @ x_pred + D @ U[t]).item()
        S = (C @ P_pred @ C.T).item() + R**2
        innov = Y[t] - y_pred

        # Log-likelihood
        loglik += -0.5 * (np.log(2 * np.pi * S) + (innov**2) / S)

        # Kalman update
        K = (P_pred @ C.T) / S
        x_pred = x_pred + (K.flatten() * innov)
        P_pred = P_pred - K @ C @ P_pred

    return float(-loglik[0][0])

def flatten_params(A, B, C, D, Q, R, x0):
    return np.concatenate([
        A.flatten(),
        B.flatten(),
        C.flatten(),
        D.flatten(),
        Q.flatten(),
        np.atleast_1d(R).flatten(),
        x0.flatten()
    ])

def unpack_params(param_vector, shapes):
    A_shape, B_shape, C_shape, D_shape, Q_shape, R_shape, x0_shape = shapes
    idx = 0

    def next_block(shape):
        nonlocal idx
        size = np.prod(shape)
        block = param_vector[idx:idx + size].reshape(shape)
        idx += size
        return block

    A = next_block(A_shape)
    B = next_block(B_shape)
    C = next_block(C_shape)
    D = next_block(D_shape)
    Q = next_block(Q_shape)
    R = next_block(R_shape)
    x0 = next_block(x0_shape)
    return A, B, C, D, Q, R, x0

def estimate_ssm_parameters(kalman_loglik, Y, u_df,
                            A_shape, B_shape, C_shape, D_shape, Q_shape, R_shape, x0_shape,
                            start, bounds, max_retries=3, tol=1e-6):
    shapes = (A_shape, B_shape, C_shape, D_shape, Q_shape, R_shape, x0_shape)

    def wrapped_loglik(param_vec):
        A, B, C, D, Q, R, x0 = unpack_params(param_vec, shapes)
        return kalman_loglik(Y=Y, A=A, B=B, C=C, D=D, Q=Q, R=R, u_df=u_df, x0=x0)

    attempt = 0
    current_start = start.copy()

    while attempt < max_retries:
        result = minimize(
            fun=wrapped_loglik,
            x0=current_start,
            bounds=bounds,
            method='L-BFGS-B',
            options={"disp": True, "maxiter": 1000}
        )

        if result.success or result.fun < tol:
            break
        else:
            print(f"Retry {attempt + 1} failed, retrying...")
            current_start = result.x
            attempt += 1

    A, B, C, D, Q, R, x0 = unpack_params(result.x, shapes)
    return result, (A, B, C, D, Q, R, x0)


### LEGACY because not generic for all A, B, C, D, Q, R

# SSM - State Space Model 1D
def kalman_loglik_1d(par, df):
    a = par[0]
    B = np.array(par[1:4])
    c = par[4]  # include C!
    sigma1 = np.exp(par[5])  # enforce positivity
    sigma2 = np.exp(par[6])

    Y = df["Y"].values
    U = df[["Ta", "S", "I"]].values

    T = len(Y)
    x_pred = 20.0
    P_pred = 10.0

    loglik = 0.0
    for t in range(T):
        # Prediction step
        x_pred = a * x_pred + np.dot(B, U[t])
        P_pred = a**2 * P_pred + sigma1**2

        # Observation prediction
        y_pred = c * x_pred
        S = c**2 * P_pred + sigma2**2
        innov = Y[t] - y_pred

        # Log-likelihood
        loglik += -0.5 * (np.log(2 * np.pi * S) + (innov**2) / S)

        # Kalman gain
        K = P_pred * c / S
        x_pred = x_pred + K * innov
        P_pred = (1 - K * c) * P_pred

    return -loglik

# SSM - State Space Model 2D
def kalman_loglik_2d(par, df):
    A = par[0:4].reshape(2, 2)
    B = par[4:10].reshape(2, 3)
    C = par[10:12].reshape(1, 2)
    sigma1 = np.exp(par[12])  # system noise (shared scalar for simplicity)
    sigma2 = np.exp(par[13])  # observation noise

    Y = df["Y"].values
    U = df[["Ta", "S", "I"]].values
    T = len(Y)

    x_pred = np.array([20.0, 20.0])
    P_pred = np.eye(2) * 10.0

    loglik = 0.0
    for t in range(T):
        x_pred = A @ x_pred + B @ U[t]
        P_pred = A @ P_pred @ A.T + np.eye(2) * sigma1**2
        
        y_pred = C @ x_pred
        S = C @ P_pred @ C.T + sigma2**2
        innov = Y[t] - y_pred
        loglik += -0.5 * (np.log(2*np.pi*S) + (innov**2)/S)
        
        K = P_pred @ C.T / S
        x_pred = x_pred + K.flatten() * innov
        P_pred = P_pred - K @ C @ P_pred
        
    return -loglik

def estimate_1d_state_space(df):
    # start = [a, B1, B2, B3, c, log(sigma1), log(sigma2)]
    start = np.array([0.8, 0.05, 0.1, 0.1, 1.0, np.log(2), np.log(2)])
    bounds = [(-2, 2), (-5, 5), (-5, 5), (-5, 5), (0.5, 1.5),  # c ~ 1
              (1e-3, np.log(10)), (1e-3, np.log(10))]

    result = minimize(kalman_loglik_1d, start, args=(df,), bounds=bounds, method='L-BFGS-B')
    return result

def estimate_2d_state_space(df):
    start = np.array([0.8, 0.1, 0.0, 0.7,   # A
                      -0.1, 0.1, 0.1,        # B row 1
                      0.0, -0.1, 0.0,        # B row 2
                      1.0, 0.2,             # C
                      np.log(2), np.log(2)])# log sigma1, sigma2
    bounds = [(-2, 2)] * (len(start)-2) + [(1e-3,2)] * 2
    result = minimize(kalman_loglik_2d, start, args=(df,), bounds=bounds, method='L-BFGS-B')
    return result

def simulate_from_model_1d(par, df):
    a = par[0]
    B = np.array(par[1:4])
    c = par[4]
    sigma1 = np.exp(par[5])
    sigma2 = np.exp(par[6])

    U = df[["Ta", "S", "I"]].values
    T = len(U)
    X = np.zeros(T + 1)
    Y = np.zeros(T)
    X[0] = 20.0  # initial condition

    for t in range(T):
        X[t+1] = a * X[t] + np.dot(B, U[t]) + np.random.normal(0, sigma1)
        Y[t] = c * X[t+1] + np.random.normal(0, sigma2)

    return Y, X[1:]

def simulate_from_model_2d(par, df):
    A = par[0:4].reshape(2,2)
    B = np.array(par[4:10]).reshape(2,3)
    C = par[10:12].reshape(1,2)
    sigma1 = np.exp(par[12])
    sigma2 = np.exp(par[13])

    U = df[["Ta", "S", "I"]].values
    T = len(U)
    X = np.zeros((T + 1, 2))
    Y = np.zeros(T)
    X[0] = 20.0  # initial condition

    for t in range(T):
        X[t+1] = A@X[t] + np.dot(B, U[t]) + np.random.multivariate_normal(mean=np.array([0,0]).reshape(-1,), cov=np.eye(2)*sigma1)
        Y[t] = C@X[t+1] + np.random.normal(0, sigma2)

    return Y, X[1:]

def reconstruct_states(par, df):
    A = par[0:4].reshape(2, 2)
    B = par[4:10].reshape(2, 3)
    C = par[10:12].reshape(1, 2)
    sigma1 = np.exp(par[12])
    sigma2 = np.exp(par[13])
    
    Y = df["Y"].values
    U = df[["Ta", "S", "I"]].values
    T = len(Y)
    
    x_pred = np.array([20.0, 20.0])
    P_pred = np.eye(2) * 10.0
    
    states = []
    for t in range(T):
        x_pred = A @ x_pred + B @ U[t]
        P_pred = A @ P_pred @ A.T + np.eye(2) * sigma1**2
        
        y_pred = C @ x_pred
        S = C @ P_pred @ C.T + sigma2**2
        innov = Y[t] - y_pred
        
        K = P_pred @ C.T / S
        x_pred = x_pred + K.flatten() * innov
        P_pred = P_pred - K @ C @ P_pred
        states.append(x_pred.copy())
        
    return np.array(states)
