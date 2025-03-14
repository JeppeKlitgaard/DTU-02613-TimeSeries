


## ---------------------------------------------------------------------------------------
# Write an ARMA simulation
n <- 200
y <- numeric(n)
# AR1: 
# (1 - phi1 * B) * Y_t = eps_t
# 



## ---------------------------------------------------------------------------------------
# Function for plot
plotit <- function(x){
  layout(rbind(1,2:3))
  par(mar=c(3,3,1,1), mgp=c(2, 0.7,0))
  plot(x, ylab="X")
  acf(x, lag.max=50, lwd=2)
  pacf(x, lag.max=50, lwd=2)
}


## ---------------------------------------------------------------------------------------
# AR(1)
n <- 200
plotit( arima.sim(model=list(ar=c(0.8)), n=n) )


## ---------------------------------------------------------------------------------------
# MA(2)
plotit( arima.sim(list(ma=c(0.9, 0.8)), n) )


## ---------------------------------------------------------------------------------------
# arma(1,2)
plotit( arima.sim(list(ar=c(0.8),ma=c(0.9, 0.8)), n) )


## ---------------------------------------------------------------------------------------
# A simulation function for ARMA simulation, use model as arima.sim, i.e. flip sign of phi (into ar) coefficients
sim <- function(model, n, nburnin=100){
  n <- n + nburnin
  # Take the ar and ma part
  ar <- model$ar
  ma <- model$ma
  # The order (i.e. the number of lags)
  p <- length(ar)
  q <- length(ma)
  # The vector for the simulation result
  y <- numeric(n)
  # Generate the random normal values
  eps <- rnorm(n)
  # Run the simulation
  for(i in (max(p,q)+1):n){
    y[i] <- eps[i] + sum(y[i-(1:p)] * ar) + sum(eps[i-(1:q)] * ma)
  }
  # Return without the burn-in period
  return(y[(nburnin+1):n])
}


## ---------------------------------------------------------------------------------------
# Test it by comparing
model <- list(ar=c(0.4), ma=c(0.2))
set.seed(12)
sim(model, 10, nburnin=100)
set.seed(12)
arima.sim(model, 10, n.start=100)

# Non-stationary process
# Do the simulation and plot
n <- 200
model <- list(ar=c(1.01))
#arima.sim(model, n)
x <- sim(model, n)
plot(x, type="l", ylab="x", xlab="t")


## ----include=FALSE----------------------------------------------------------------------
lagvec <- function(x, lag){
    if (lag > 0) {
        ## Lag x, i.e. delay x lag steps
        return(c(rep(NA, lag), x[1:(length(x) - lag)]))
    }else if(lag < 0) {
        ## Lag x, i.e. delay x lag steps
        return(c(x[(abs(lag) + 1):length(x)], rep(NA, abs(lag))))
    }else{
        ## lag = 0, return x
        return(x)
    }
}
lagdf <- function(x, lagseq) {
    ## Return a data.frame
    tmp <- as.data.frame(do.call("cbind", lapply(lagseq, function(lag){
        return(lagvec(x, lag))
    })))
    names(tmp) <- paste0("k",lagseq)
    return(tmp)
}


## ---------------------------------------------------------------------------------------
# Test it by comparing
model <- list(ar=c(0.4))
set.seed(12)
sim(model, 10, nburnin=100)
set.seed(12)
x <- arima.sim(model, 100)

X <- lagdf(x, 0:3)
summary(lm(k0 ~ k1, X))
summary(lm(k0 ~ k1 + k2, X))
summary(lm(k0 ~ k1 + k2 + k3, X))


