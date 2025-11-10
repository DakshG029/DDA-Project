set.seed(101)

# ------------------------------
# 1) Load libraries
# ------------------------------
library(quantmod)
library(fitdistrplus)
library(evir)

# ------------------------------
# 2) Prepare data (NIFTY50)
# ------------------------------
ticker <- "^NSEI"
start_date <- as.Date("2010-01-01")
end_date <- Sys.Date()
dataset_xts <- getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE)
dataset <- data.frame(Date = index(dataset_xts), coredata(dataset_xts))
dataset <- na.omit(dataset)

prices <- dataset$NSEI.Close
returns <- diff(log(prices))
losses <- -as.numeric(returns)    # working in LOSSES
losses <- losses[!is.na(losses)]

pk <- 0.95      # confidence level for VaR/CVaR
pu <- 0.90      # tail threshold
pu_thresh <- quantile(losses, pu)
exceedances <- losses[losses > pu_thresh] - pu_thresh

# ------------------------------
# 3) Helper: GPD log-likelihood + VaR formula
# ------------------------------
gpd_loglik <- function(xi, beta, y) {
  if (beta <= 0) return(-Inf)
  n <- length(y)
  if (abs(xi) < 1e-8) {
    ll <- -n * log(beta) - sum(y)/beta
  } else {
    if (any(1 + xi * y / beta <= 0)) return(-Inf)
    ll <- -n * log(beta) - (1/xi + 1) * sum(log(1 + xi * y / beta))
  }
  return(ll)
}

gpd_VaR_from_params <- function(threshold, xi, beta, pu, pk) {
  q <- (pk - pu)/(1 - pu)
  if (abs(xi) < 1e-8) {
    return(threshold - beta * log(1 - q))
  } else {
    return(threshold + beta/xi * ((1 - q)^(-xi) - 1))
  }
}

# ------------------------------
# 4) Method 1: Empirical
# ------------------------------
VaR_emp <- quantile(losses, pk)
CVaR_emp <- mean(losses[losses >= VaR_emp])

# ------------------------------
# 5) Method 2: Bootstrap
# ------------------------------
B <- 500
VaR_boot <- numeric(B)
CVaR_boot <- numeric(B)
for (i in 1:B) {
  sample_X <- sample(losses, length(losses), replace = TRUE)
  VaR_boot[i] <- quantile(sample_X, pk)
  CVaR_boot[i] <- mean(sample_X[sample_X >= VaR_boot[i]])
}
VaR_bootstrap <- 2 * VaR_emp - mean(VaR_boot)
CVaR_bootstrap <- 2 * CVaR_emp - mean(CVaR_boot)

# ------------------------------
# 6) Method 3: Parametric (t-distribution) [FIXED]
# ------------------------------
library(MASS)

# Fit a Student-t distribution using MASS::fitdistr
fit_t <- fitdistr(losses, densfun = "t", start = list(m = mean(losses), s = sd(losses), df = 5))

m_t <- fit_t$estimate["m"]
s_t <- fit_t$estimate["s"]
df_t <- fit_t$estimate["df"]

VaR_par <- m_t + s_t * qt(pk, df = df_t)
CVaR_par <- m_t + s_t * (dt(qt(pk, df_t), df_t) / (1 - pk)) * ((df_t + qt(pk, df_t)^2) / (df_t - 1))


# ------------------------------
# 7) Method 4: Parametric Bayesian (t-dist priors)
# ------------------------------
n_iter <- 5000
mu_chain <- numeric(n_iter)
sigma_chain <- numeric(n_iter)
df_chain <- numeric(n_iter)
mu_chain[1] <- mean(losses)
sigma_chain[1] <- sd(losses)
df_chain[1] <- 5

for (i in 2:n_iter) {
  mu_prop <- rnorm(1, mu_chain[i-1], 0.02)
  sigma_prop <- abs(rnorm(1, sigma_chain[i-1], 0.02))
  df_prop <- abs(rnorm(1, df_chain[i-1], 0.5))
  
  ll_curr <- sum(log(dt((losses - mu_chain[i-1]) / sigma_chain[i-1], df_chain[i-1]) / sigma_chain[i-1]))
  ll_prop <- sum(log(dt((losses - mu_prop) / sigma_prop, df_prop) / sigma_prop))
  
  prior_curr <- dnorm(mu_chain[i-1], mean(losses), sd(losses), log = TRUE) +
    dgamma(sigma_chain[i-1], 2, 2, log = TRUE)
  prior_prop <- dnorm(mu_prop, mean(losses), sd(losses), log = TRUE) +
    dgamma(sigma_prop, 2, 2, log = TRUE)
  
  log_acc <- (ll_prop + prior_prop) - (ll_curr + prior_curr)
  if (log(runif(1)) < log_acc) {
    mu_chain[i] <- mu_prop
    sigma_chain[i] <- sigma_prop
    df_chain[i] <- df_prop
  } else {
    mu_chain[i] <- mu_chain[i-1]
    sigma_chain[i] <- sigma_chain[i-1]
    df_chain[i] <- df_chain[i-1]
  }
}
mu_post <- mean(mu_chain)
sigma_post <- mean(sigma_chain)
df_post <- mean(df_chain)

VaR_parB <- mu_post + sigma_post * qt(pk, df_post)
CVaR_parB <- mu_post + sigma_post * (dt(qt(pk, df_post), df_post) / (1 - pk)) * ((df_post + qt(pk, df_post)^2) / (df_post - 1))

# ------------------------------
# 8) Method 5: Extreme Value Bayesian (EVB)
# ------------------------------
gpd_fit <- gpd(losses, pu_thresh)
xi_evb <- gpd_fit$par.ests["xi"]
beta_evb <- gpd_fit$par.ests["beta"]

VaR_evb <- gpd_VaR_from_params(pu_thresh, xi_evb, beta_evb, pu, pk)
q <- (pk - pu)/(1 - pu)
CVaR_evb <- pu_thresh + (beta_evb / (1 - xi_evb)) * ((1 - q)^(-xi_evb) - 1)

# ------------------------------
# 9) Method 6: Informative Prior Bayesian (IPB)
# ------------------------------
mu_xi_prior <- 0.1; sd_xi_prior <- 0.06
mean_beta_prior <- max(mean(exceedances, na.rm = TRUE), 1e-6)
sd_beta_prior <- 0.5 * mean_beta_prior

n_mcmc2 <- 8000; burn_mcmc2 <- 3000
xi_chain2 <- numeric(n_mcmc2); beta_chain2 <- numeric(n_mcmc2)
xi_chain2[1] <- mu_xi_prior
beta_chain2[1] <- mean_beta_prior
sd_xi_prop2 <- 0.03
sd_beta_prop2 <- 0.1 * beta_chain2[1]
n_sim_tail <- 5000

for (i in 2:n_mcmc2) {
  xi_prop <- rnorm(1, xi_chain2[i-1], sd_xi_prop2)
  beta_prop <- rnorm(1, beta_chain2[i-1], sd_beta_prop2)
  if (beta_prop <= 0) { xi_chain2[i] <- xi_chain2[i-1]; beta_chain2[i] <- beta_chain2[i-1]; next }
  
  logprior_curr <- dnorm(xi_chain2[i-1], mu_xi_prior, sd_xi_prior, log = TRUE) -
    dgamma(beta_chain2[i-1], shape = (mean_beta_prior/sd_beta_prior)^2, rate = mean_beta_prior/(sd_beta_prior^2), log = TRUE)
  logprior_prop <- dnorm(xi_prop, mu_xi_prior, sd_xi_prior, log = TRUE) -
    dgamma(beta_prop, shape = (mean_beta_prior/sd_beta_prior)^2, rate = mean_beta_prior/(sd_beta_prior^2), log = TRUE)
  
  ll_curr <- gpd_loglik(xi_chain2[i-1], beta_chain2[i-1], exceedances)
  ll_prop <- gpd_loglik(xi_prop, beta_prop, exceedances)
  
  logacc <- (ll_prop + logprior_prop) - (ll_curr + logprior_curr)
  if (log(runif(1)) < logacc) {
    xi_chain2[i] <- xi_prop; beta_chain2[i] <- beta_prop
  } else {
    xi_chain2[i] <- xi_chain2[i-1]; beta_chain2[i] <- beta_chain2[i-1]
  }
}

xi_post2 <- xi_chain2[(burn_mcmc2+1):n_mcmc2]
beta_post2 <- beta_chain2[(burn_mcmc2+1):n_mcmc2]
n_post2 <- length(xi_post2)

VaR_ibd_draws <- numeric(n_post2)
CVaR_ibd_draws <- numeric(n_post2)
for (i in 1:n_post2) {
  xi_i <- xi_post2[i]; beta_i <- beta_post2[i]
  VaR_i <- gpd_VaR_from_params(pu_thresh, xi_i, beta_i, pu, pk)
  U <- runif(n_sim_tail)
  if (abs(xi_i) < 1e-8) {
    Ysim <- -beta_i * log(1 - U)
  } else {
    Ysim <- beta_i/xi_i * ((1 - U)^(-xi_i) - 1)
  }
  Xsim <- pu_thresh + Ysim
  CVaR_i <- mean(Xsim[Xsim >= VaR_i])
  VaR_ibd_draws[i] <- VaR_i
  CVaR_ibd_draws[i] <- CVaR_i
}
VaR_ipb <- mean(VaR_ibd_draws, na.rm = TRUE)
CVaR_ipb <- mean(CVaR_ibd_draws, na.rm = TRUE)

# ------------------------------
# 10) Compile all results
# ------------------------------
results <- data.frame(
  Method = c("Empirical", "Bootstrap", "Parametric (t)", 
             "Parametric Bayesian", "Extreme Value Bayesian", "Informative Prior Bayesian"),
  VaR = c(VaR_emp, VaR_bootstrap, VaR_par, VaR_parB, VaR_evb, VaR_ipb),
  CVaR = c(CVaR_emp, CVaR_bootstrap, CVaR_par, CVaR_parB, CVaR_evb, CVaR_ipb)
)

print(results)
