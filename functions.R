#########################
# 0. Packages & options #
#########################

required_pkgs <- c("bayesforecast", "scoringRules", "parallel")

missing <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "The following packages are not installed: ",
    paste(missing, collapse = ", "), "\n",
    "Please install them manually first, e.g.:\n",
    "install.packages(c('bayesforecast','scoringRules','parallel'))\n"
  )
}

# Load packages
for (pkg in required_pkgs) {
  library(pkg, character.only = TRUE)
}

# rstan options (used internally by bayesforecast)
if (requireNamespace("rstan", quietly = TRUE)) {
  rstan::rstan_options(auto_write = TRUE)
  # We set 1 core per Stan chain and parallelize at the R level
  options(mc.cores = 1L)
}

# Number of OS cores for R-level parallelisation (Monte Carlo replications)
N_CORES <- min(32L, parallel::detectCores())
cat("Detected", parallel::detectCores(), "cores; using up to", N_CORES, "for Monte Carlo.\n")

#########################
# 1. Utility functions  #
#########################

## 1.1 Shifted-gamma estimator of marginal log-likelihood
##     (Raftery et al. style), following the implementation
##     used in the article.

marginal_loglik_shifted_gamma <- function(loglik_samples, d0, lambda = 0.98) {
  loglik_samples <- as.numeric(loglik_samples)
  l_bar <- mean(loglik_samples)
  s2    <- stats::var(loglik_samples)
  
  # Shifted-gamma trick: choose a high reference level l0_hat
  l0_hat <- max(l_bar + s2, max(loglik_samples))
  alpha0 <- d0 / 2
  # Approximate log marginal likelihood
  log_py <- l0_hat + alpha0 * log(1 - lambda)
  log_py
}

## 1.2 Build concentration curve vs uniform reference
##     Given a vector of non-negative weights prob (e.g., predictive
##     densities at forecast points), we normalise it and compute
##     the Lorenz-type curve relative to a uniform measure:
##       Q_k = k / K (uniform cumulative share, K = length(prob))
##       P_k = cumulative sum of prob sorted in ascending order.

build_concentration_uniform <- function(prob) {
  prob <- as.numeric(prob)
  prob <- prob / sum(prob)
  K    <- length(prob)
  
  ord      <- order(prob)      # ascending order
  p_sorted <- prob[ord]
  P        <- cumsum(p_sorted) # cumulative mass
  Q        <- (1:K) / K        # uniform reference (cumulative mass)
  
  list(P = P, Q = Q)
}

## 1.3 Build concentration curve of p relative to q
##     (Cifarelli–Regazzini concentration function in discrete form)
##
##     Inputs:
##       - p: vector of non-negative weights for model M1
##       - q: vector of non-negative weights for model M0 (reference)
##     Both are normalised internally. The ordering is by the
##     predictive likelihood ratio r_k = p_k / q_k as described in
##     the article.

build_concentration_relative <- function(p, q) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  if (length(p) != length(q)) stop("p and q must have the same length.")
  p <- p / sum(p)
  q <- q / sum(q)
  
  if (any(q <= 0)) {
    stop("Reference distribution q has zero or non-positive entries.")
  }
  
  r   <- p / q
  ord <- order(r)         # increasing likelihood ratio
  
  p_ord <- p[ord]
  q_ord <- q[ord]
  
  P <- cumsum(p_ord)
  Q <- cumsum(q_ord)
  
  list(P = P, Q = Q, r_ord = r[ord])
}

## 1.4 Compute Gini coefficient, Pietra index and concentration area
##     from a discrete concentration curve (P vs Q).
##
##     - Gini coefficient G = 2 * area between the equality line
##       and the concentration curve, i.e. G in [0,1].
##     - Pietra index C = max (Q_k - P_k).
##
##     Here, "concentration area" = area between the line y=x and
##     the curve (in [0, 0.5]).

compute_gini_pietra <- function(P, Q) {
  P <- as.numeric(P)
  Q <- as.numeric(Q)
  if (length(P) != length(Q)) stop("P and Q must have the same length.")
  
  K <- length(P)
  
  # Normalise to end at 1 (for numerical stability)
  P <- P / P[K]
  Q <- Q / Q[K]
  
  # Add (0,0) point at the beginning
  P_full <- c(0, P)
  Q_full <- c(0, Q)
  
  dQ <- diff(Q_full)
  
  # Area under the concentration curve (trapezoidal rule)
  area_under <- sum((P_full[-(K + 1)] + P_full[-1]) / 2 * dQ)
  
  # Area under the equality line y=x on [0,1] is 1/2
  area_line <- 0.5
  area_between <- area_line - area_under
  
  # Gini coefficient
  G <- 2 * area_between
  
  # Pietra index (max vertical distance)
  C <- max(Q - P)
  
  list(Gini = G, Pietra = C, Area = area_between)
}

## 1.5 CRPS for a mixture of Normals via Monte Carlo
##     We generate draws from the posterior predictive mixture
##     (y_rep) and use scoringRules::crps_sample.

crps_mixture_normal <- function(y, means, sds, n_draws = 1000L) {
  M <- length(means)
  if (M != length(sds)) stop("means and sds must have the same length.")
  idx   <- sample.int(M, size = n_draws, replace = TRUE)
  draws <- stats::rnorm(n_draws, mean = means[idx], sd = sds[idx])
  scoringRules::crps_sample(y = y, dat = draws)
}

## 1.6 One-step-ahead Bayesian predictive metrics for ARMA(p,q)
##     fitted via bayesforecast::stan_sarima.
##
##     Inputs:
##       - fit      : stan_sarima fitted object
##       - y_train  : numeric vector of training data (length t)
##       - y_future : realised value at horizon 1 (y_{t+1})
##       - p, q     : AR and MA orders
##       - d0       : number of parameters (p + q + 2 typically)
##       - n_crps   : number of predictive draws for CRPS approximation
##
##     Outputs:
##       - pred_density : Monte Carlo estimate of predictive density at y_future
##       - nlpd         : -log(pred_density) (NLPD as in the article)
##       - crps         : CRPS via posterior predictive draws
##       - log_marginal : marginal log-likelihood via shifted-gamma estimator

predictive_metrics_bayes_arma <- function(fit, y_train, y_future, p, q,
                                          d0, n_crps = 1000L) {
  
  # Extract posterior draws from bayesforecast::extract_stan()
  post <- bayesforecast::extract_stan(fit)
  
  mu0   <- as.numeric(post$mu0)     # intercept draws
  sigma <- as.numeric(post$sigma0)  # innovation sd draws
  M     <- length(mu0)
  
  # AR and MA coefficients (if present)
  phi <- if (!is.null(post$ar)) post$ar else NULL
  theta <- if (!is.null(post$ma)) post$ma else NULL
  
  # Residuals (innovations) from the fitted model
  resids <- as.numeric(residuals(fit))
  t      <- length(y_train)
  
  # Last p observed values: y_t, y_{t-1},...,y_{t-p+1}
  if (p > 0) {
    if (t < p) stop("Not enough observations for AR lags.")
    y_lags <- y_train[t:(t - p + 1)]
  } else {
    y_lags <- numeric(0)
  }
  
  # Last q residuals: e_t, e_{t-1},...,e_{t-q+1}
  if (q > 0) {
    if (length(resids) < q) stop("Not enough residuals for MA lags.")
    e_lags <- resids[t:(t - q + 1)]
  } else {
    e_lags <- numeric(0)
  }
  
  # Conditional mean for each posterior draw:
  #   y_{t+1} = mu0 + sum_j phi_j y_{t+1-j} - sum_j theta_j e_{t+1-j} + u_{t+1}
  # The minus sign in the MA part matches the sign convention used
  # in the article's loglik_arma_pq implementation.
  mu_draw <- mu0
  
  if (p > 0) {
    phi_mat <- as.matrix(phi)
    if (nrow(phi_mat) != M) phi_mat <- t(phi_mat)
    mu_draw <- mu_draw + as.numeric(phi_mat %*% y_lags)
  }
  
  if (q > 0) {
    theta_mat <- as.matrix(theta)
    if (nrow(theta_mat) != M) theta_mat <- t(theta_mat)
    mu_draw <- mu_draw - as.numeric(theta_mat %*% e_lags)
  }
  
  # Predictive density at y_future as mixture of normals
  dens_draw     <- stats::dnorm(y_future, mean = mu_draw, sd = sigma)
  pred_density  <- mean(dens_draw)
  nlpd          <- -log(pred_density)  # NLPD as -log predictive density
  
  # CRPS via Monte Carlo from the posterior predictive mixture
  if (n_crps > 0L) {
    crps_val <- crps_mixture_normal(y_future, mu_draw, sigma, n_draws = n_crps)
  } else {
    crps_val <- NA_real_
  }
  
  # Marginal log-likelihood via shifted-gamma estimator
  loglik_samples <- as.numeric(post$loglik)
  log_marginal   <- marginal_loglik_shifted_gamma(loglik_samples, d0 = d0)
  
  list(
    pred_density = pred_density,
    nlpd         = nlpd,
    crps         = crps_val,
    log_marginal = log_marginal
  )
}
