##############################################################
# 3. Real data: US industrial natural gas prices (EIA)       #
##############################################################

cat("\n### Part 2: Real data application (US natural gas prices) ###\n")

## 3.1 Hard-coded series 2010M01–2025M06
##     IMPORTANT: If you want exact match with the article, make sure
##     these values correspond exactly to the EIA series used there.

gas_prices <- c(
  # 2010
  6.93, 6.76, 6.01, 5.12, 5.08, 5.04, 5.49, 5.37, 4.61, 4.73, 4.60, 5.50,
  # 2011
  5.66, 5.77, 5.21, 5.34, 5.21, 5.21, 5.05, 5.21, 4.84, 4.71, 4.64, 4.59,
  # 2012
  4.58, 4.19, 3.71, 3.21, 3.02, 3.34, 3.60, 3.83, 3.56, 3.94, 4.46, 4.73,
  # 2013
  4.58, 4.54, 4.59, 4.95, 5.00, 4.90, 4.47, 4.31, 4.36, 4.36, 4.62, 4.97,
  # 2014
  5.69, 6.63, 6.47, 5.85, 5.74, 5.46, 5.43, 4.96, 5.02, 5.03, 5.02, 5.62,
  # 2015
  4.90, 4.74, 4.46, 3.96, 3.58, 3.76, 3.74, 3.79, 3.65, 3.54, 3.28, 3.48,
  # 2016
  3.62, 3.58, 3.02, 3.00, 2.90, 2.89, 3.57, 3.59, 3.74, 3.87, 3.86, 4.27,
  # 2017
  4.85, 4.53, 3.92, 4.11, 4.02, 4.05, 3.92, 3.78, 3.83, 3.78, 3.84, 4.19,
  # 2018
  4.46, 4.85, 4.00, 3.89, 3.80, 3.77, 3.75, 3.67, 3.75, 4.03, 4.51, 5.47,
  # 2019
  5.02, 4.62, 4.31, 3.99, 3.64, 3.55, 3.33, 3.18, 3.35, 3.43, 3.86, 3.84,
  # 2020
  3.71, 3.58, 3.39, 3.00, 2.91, 2.72, 2.58, 2.85, 3.30, 3.29, 3.98, 4.11,
  # 2021
  4.04, 9.32, 4.41, 4.00, 4.11, 4.16, 4.69, 4.95, 5.42, 6.61, 6.90, 6.77,
  # 2022
  6.49, 7.34, 6.20, 6.70, 8.11, 9.34, 7.89, 9.44, 9.62, 7.18, 6.76, 8.08,
  # 2023
  7.18, 5.95, 5.00, 4.04, 3.54, 3.52, 3.84, 3.80, 3.81, 4.05, 4.35, 4.48,
  # 2024
  5.05, 4.80, 3.76, 3.35, 3.18, 3.70, 3.61, 3.10, 3.28, 3.81, 3.92, 5.05,
  # 2025 (Jan–Jun)
  5.83, 5.74, 5.48, 5.10, 4.51, 4.46
)

n_real <- length(gas_prices)
cat("Length of gas price series:", n_real, "observations.\n")

T0_real <- floor(2 * n_real / 3)
K_real  <- n_real - T0_real

# Stan iterations for real-data models (same as simulation for fidelity)
N_ITER_REAL <- N_ITER_MCMC

# 3.2 Candidate models (as in the article)
candidate_models <- list(
  list(label = "M0: AR(1)",      p = 1, q = 0),
  list(label = "M1: AR(2)",      p = 2, q = 0),
  list(label = "M2: ARMA(1,1)",  p = 1, q = 1),
  list(label = "M3: ARMA(1,2)",  p = 1, q = 2),
  list(label = "M4: ARMA(2,1)",  p = 2, q = 1),
  list(label = "M5: MA(1)",      p = 0, q = 1),
  list(label = "M6: AR(3)",      p = 3, q = 0)
)

J <- length(candidate_models)

dens_mat <- matrix(NA_real_, nrow = K_real, ncol = J)
nlpd_mat <- matrix(NA_real_, nrow = K_real, ncol = J)
crps_mat <- matrix(NA_real_, nrow = K_real, ncol = J)
ml_mat   <- matrix(NA_real_, nrow = K_real, ncol = J)

cat("Rolling-origin Bayesian ARMA fitting for real data...\n")

for (t in T0_real:(n_real - 1L)) {
  k <- t - T0_real + 1L
  y_train  <- gas_prices[1:t]
  y_future <- gas_prices[t + 1L]
  
  cat(sprintf("  - Forecast origin t = %d (k = %d of %d)\n", t, k, K_real))
  
  for (j in seq_len(J)) {
    spec <- candidate_models[[j]]
    p <- spec$p
    q <- spec$q
    d0 <- p + q + 2
    
    fit_j <- bayesforecast::stan_sarima(
      ts     = y_train,
      order  = c(p, 0, q),
      chains = 1,
      iter   = N_ITER_REAL
    )
    
    m_j <- predictive_metrics_bayes_arma(
      fit_j, y_train, y_future,
      p = p, q = q, d0 = d0,
      n_crps = 1000L  # CRPS via mixture draws
    )
    
    dens_mat[k, j] <- m_j$pred_density
    nlpd_mat[k, j] <- m_j$nlpd
    crps_mat[k, j] <- m_j$crps
    ml_mat[k, j]   <- m_j$log_marginal
  }
}

## 3.3 Concentration vs equality line (uniform reference)

concentration_unif <- vector("list", J)
gini_unif   <- numeric(J)
pietra_unif <- numeric(J)

for (j in seq_len(J)) {
  p_j <- dens_mat[, j]
  conc_j <- build_concentration_uniform(p_j)
  concentration_unif[[j]] <- conc_j
  gp_j <- compute_gini_pietra(conc_j$P, conc_j$Q)
  gini_unif[j]   <- gp_j$Gini
  pietra_unif[j] <- gp_j$Pietra
}

## 3.4 Relative concentration of each model vs M0 (AR(1))

p_M0 <- dens_mat[, 1] / sum(dens_mat[, 1])  # reference: AR(1)

gini_rel   <- rep(NA_real_, J)
pietra_rel <- rep(NA_real_, J)

for (j in 2:J) {
  p_j <- dens_mat[, j]
  conc_rel_j <- build_concentration_relative(p_j, p_M0)
  gp_rel_j   <- compute_gini_pietra(conc_rel_j$P, conc_rel_j$Q)
  gini_rel[j]   <- gp_rel_j$Gini
  pietra_rel[j] <- gp_rel_j$Pietra
}

## 3.5 Average NLPD, CRPS and marginal log-likelihood for each model

avg_nlpd <- colMeans(nlpd_mat)
avg_crps <- colMeans(crps_mat)
avg_ml   <- colMeans(ml_mat)

## 3.6 Summary tables analogous to the article

real_summary <- data.frame(
  Model   = vapply(candidate_models, `[[`, character(1), "label"),
  Pietra_index = pietra_unif,
  Gini_concentration_coefficient = gini_unif,
  NLPD   = avg_nlpd,
  CRPS   = avg_crps,
  Marginal_log_likelihood = avg_ml,
  row.names = NULL
)

cat("\nReal-data summary (concentration vs equality line):\n")
print(real_summary)

relative_summary <- data.frame(
  Model   = vapply(candidate_models, `[[`, character(1), "label"),
  Pietra_relative = pietra_rel,
  Gini_relative   = gini_rel,
  row.names = NULL
)

cat("\nPietra and Gini indices of models M1–M6 relative to M0 (AR(1)):\n")
print(relative_summary)
