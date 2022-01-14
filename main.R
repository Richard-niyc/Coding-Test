# Note: I learned the following backtesting and portfolio construction method in a course named "Portfolio Optimization in R" taught by Professor PALOMAR, P. Daniel last semester at HKUST.

library(portfolioBacktest)
library(riskParityPortfolio)
library(gridExtra)
library(tidyr)
library(xts)

## Firstly, design strategy: the general Risk Parity Portfolio using Tyler estimator with shrinkage

# let's define the Tyler estimator with shrinkage
estimateTylerShrinkage <- function(X, R_target, alpha, N, verbose = FALSE) {
  max_iter <- 100
  error_th_Sigma <- 1e-3
  
  # Gaussian initial point
  Sigma <- 1/(1+alpha)*cov(X) + alpha/(1+alpha)*R_target  # <--- this line is new
  Sigma <- Sigma/sum(diag(Sigma))
  # loop
  obj_value_record <- Sigma_diff_record <- rep(NA, max_iter)
  for (k in 1:max_iter) {
    Sigma_prev <- Sigma
    
    # Tyler update
    weights <- 1/rowSums(X * (X %*% solve(Sigma)))   # 1/diag( X %*% inv(Sigma) %*% t(X) )
    obj_value_record[k] <- - (N/2)*sum(log(weights)) + (T/2)*sum(log(eigen(Sigma)$values))
    Sigma <- (N/T) * crossprod( sqrt(weights)*X )  # (N/T) * t(X) %*% diag(weights) %*% X
    Sigma <- 1/(1+alpha)*Sigma + alpha/(1+alpha)*R_target  # <--- this line is new
    Sigma <- Sigma/sum(diag(Sigma))
    
    # stopping criterion
    Sigma_diff_record[k] <- norm(Sigma - Sigma_prev, "F")/norm(Sigma_prev, "F")
    if (Sigma_diff_record[k] < error_th_Sigma)
      break
  }
  obj_value_record <- obj_value_record[1:k]
  Sigma_diff_record <- Sigma_diff_record[1:k]
  if (verbose)
    plot(obj_value_record, type = "l", col = "blue",
         main = "Convergence of objective value", xlab = "iterations", ylab = "obj value")
  
  # finally, recover missing scaling factor
  sigma2 <- apply(X, 2, var)
  d <- diag(Sigma)
  kappa <- sum(sigma2*d)/sum(d*d)
  Sigma <- kappa*Sigma
  
  return(Sigma)
}

# construct the general Risk Parity Portfolio using Tyler estimator with shrinkage
RPP_mu_TylerShrinkage <- function(data, w_current) {
  X <- diff(log(data$last))[-1]  # compute log returns
  mu_sm <- colMeans(X)
  N <- ncol(X)
  
  # recall the error of the sample covariance matrix with shrinkage
  R_target <- mean(diag(cov(X))) * diag(N)  # shrinkage target
  alpha <- 0.5  #shrinkage factor
  
  R_Tyler_shrinked <- estimateTylerShrinkage(X, R_target, alpha, N)
  rpp_mu <- riskParityPortfolio(R_Tyler_shrinked, mu = mu_sm, lmd_mu = 5e-5, formulation = "rc-double-index")
  return(rpp_mu$w)
}


## Secondly, preprocess data so that it can be input in portfolioBacktest function

# read data
data <- read.csv("data.csv", header = TRUE)

# convert the shape and type of data
last_df <- tidyr::spread(data = data[, -4], key = ticker, value = last)
volume_df <- tidyr::spread(data = data[, -3], key = ticker, value = volume)

last_xts <- xts(last_df[, -1], order.by = as.Date(last_df[, 1]))
volume_xts <- xts(volume_df[, -1], order.by = as.Date(volume_df[, 1]))

# combine two xts data
listdata <- list(last_xts, volume_xts)
names(listdata) <- c("last", "volume")

# resample the data into 100 datasets, each contains 50 stocks' data in 2 years randomly
backtest_data <- financialDataResample(listdata, N_sample = 50, T_sample = 252*2, num_datasets = 100)


## Thirdly, backtest and compare my strategy to the benchmark (uniform portfolio)

portfolios <- list("RPP_mu_TylerShrinkage" = RPP_mu_TylerShrinkage)
backtest <- portfolioBacktest(portfolio_funs = portfolios,
                              dataset_list = backtest_data,
                              price_name = "last",
                              paral_portfolios = 1,
                              paral_datasets = 2,
                              show_progress_bar = TRUE,
                              benchmark = "uniform",
                              shortselling = TRUE,
                              leverage = 1,
                              lookback = 252,
                              optimize_every = 20,
                              rebalance_every = 1,
                              bars_per_year = 252)

# show the performance summary
backtestSummary(backtest)$performance
grid.table(backtestSummary(backtest)$performance)

# check error
# res <- backtestSelector(backtest, portfolio_index = 1)
# 
# # information of 1st error
# error1 <- res$error_message[[1]]
# str(error1)


## Finally, determine if my strategy is better than the benchmark (uniform portfolio)

# design the ranking criteria (Since I can add more strategies to compare if it is needed. But in this case, we just backtest my strategy and the benchmark.)
leaderboard <- backtestLeaderboard(backtest, 
                                   weights = list("Sharpe ratio"  = 7, 
                                                  "max drawdown"  = 1, 
                                                  "annual return" = 1, 
                                                  "ROT (bps)"     = 1))

# show leaderboard
grid.table(leaderboard$leaderboard_scores)

# In my 100 times of backtesting, the strategy works very well! 
# But since I resampled the data randomly, it will perform differently for each backtesting.