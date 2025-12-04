###########################################
# 2. Monte Carlo study: AR(1) vs MA(1)    #
###########################################

cat("\n### Part 1: Monte Carlo experiment (DGP = ARMA(1,1), M0 = AR(1), M1 = MA(1)) ###\n")

set.seed(1234)

# 2.1 Settings (as in the article)
N_REP    <- 100L         # number of Monte Carlo replications
N        <- 100L         # length of each time series
PHI_DGP   <- 0.2         # AR coefficient in DGP
THETA_DGP <- 0.5         # MA coefficient in DGP
SIGMA_DGP <- sqrt(0.16)  # innovation sd in DGP

T0 <- floor(2 * N / 3)   # initial training window (~2/3)
K  <- N - T0             # number of forecast points

# Stan settings (as in the article: 15000 iterations, 7500 warm-up)
N_ITER_MCMC <- 15000L
# For quick testing, you may temporarily reduce this number, e.g. 2000

# 2.2 Worker function for a single Monte Carlo replication

simulate_single_rep_bayes <- function(rep_id) {
  cat(sprintf("  - Replication %d\n", rep_id))
  
  set.seed(1000 + rep_id)
  
  # Generate ARMA(1,1) time series
  y <- as.numeric(stats::arima.sim(
    n    = N,
    list(ar = PHI_DGP, ma = THETA_DGP),
    sd   = SIGMA_DGP
  ))
  
  # M0: AR(1), M1: MA(1)
  p0 <- 1; q0 <- 0
  p1 <- 0; q1 <- 1
  
  d0_M0 <- p0 + q0 + 2  # number of parameters: ARs + MAs + mu + sigma
  d0_M1 <- p1 + q1 + 2
  
  dens_M0 <- numeric(K)  # predictive density for M0 at each forecast
  dens_M1 <- numeric(K)  # predictive density for M1
  nlpd_M0 <- numeric(K)  # NLPD for M0
  nlpd_M1 <- numeric(K)  # NLPD for M1
  ml_M0   <- numeric(K)  # marginal log-likelihood (approx) for M0
  ml_M1   <- numeric(K)  # marginal log-likelihood (approx) for M1
  
  # Rolling-origin forecasting
  for (t in T0:(N - 1L)) {
    k <- t - T0 + 1L
    y_train  <- y[1:t]
    y_future <- y[t + 1L]
    
    # Fit Bayesian AR(1)
    fit0 <- bayesforecast::stan_sarima(
      ts     = y_train,
      order  = c(p0, 0, q0),
      chains = 1,
      iter   = N_ITER_MCMC
    )
    
    # Fit Bayesian MA(1)
    fit1 <- bayesforecast::stan_sarima(
      ts     = y_train,
      order  = c(p1, 0, q1),
      chains = 1,
      iter   = N_ITER_MCMC
    )
    
    m0 <- predictive_metrics_bayes_arma(
      fit0, y_train, y_future, p = p0, q = q0,
      d0 = d0_M0, n_crps = 0L  # CRPS not needed in Monte Carlo study
    )
    
    m1 <- predictive_metrics_bayes_arma(
      fit1, y_train, y_future, p = p1, q = q1,
      d0 = d0_M1, n_crps = 0L
    )
    
    dens_M0[k] <- m0$pred_density
    dens_M1[k] <- m1$pred_density
    nlpd_M0[k] <- m0$nlpd
    nlpd_M1[k] <- m1$nlpd
    ml_M0[k]   <- m0$log_marginal
    ml_M1[k]   <- m1$log_marginal
  }
  
  # Build predictive distributions over the forecast window
  p <- dens_M1 / sum(dens_M1)  # M1
  q <- dens_M0 / sum(dens_M0)  # M0
  
  # Concentration vs uniform for each model
  conc_M0_unif <- build_concentration_uniform(q)
  conc_M1_unif <- build_concentration_uniform(p)
  
  gp_M0 <- compute_gini_pietra(conc_M0_unif$P, conc_M0_unif$Q)
  gp_M1 <- compute_gini_pietra(conc_M1_unif$P, conc_M1_unif$Q)
  
  # Relative concentration of M1 relative to M0
  conc_rel <- build_concentration_relative(p, q)
  gp_rel   <- compute_gini_pietra(conc_rel$P, conc_rel$Q)
  
  list(
    conc_M0_P = conc_M0_unif$P,
    conc_M0_Q = conc_M0_unif$Q,
    conc_M1_P = conc_M1_unif$P,
    conc_M1_Q = conc_M1_unif$Q,
    conc_rel_P = conc_rel$P,
    conc_rel_Q = conc_rel$Q,
    gini_M0   = gp_M0$Gini,
    pietra_M0 = gp_M0$Pietra,
    gini_M1   = gp_M1$Gini,
    pietra_M1 = gp_M1$Pietra,
    gini_rel   = gp_rel$Gini,
    pietra_rel = gp_rel$Pietra,
    mean_nlpd_M0 = mean(nlpd_M0),
    mean_nlpd_M1 = mean(nlpd_M1),
    mean_ml_M0   = mean(ml_M0),
    mean_ml_M1   = mean(ml_M1)
  )
}

# 2.3 Run all Monte Carlo replications (parallel on Linux/macOS, sequential on Windows)

use_parallel <- (.Platform$OS.type != "windows" && N_CORES > 1L)

if (use_parallel) {
  sim_list <- parallel::mclapply(
    X        = 1:N_REP,
    FUN      = simulate_single_rep_bayes,
    mc.cores = N_CORES
  )
} else {
  cat("Running Monte Carlo study sequentially (no mclapply on this OS).\n")
  sim_list <- lapply(1:N_REP, simulate_single_rep_bayes)
}

# 2.4 Collect Monte Carlo results

P_M0_mat <- do.call(cbind, lapply(sim_list, `[[`, "conc_M0_P"))  # K x N_REP
P_M1_mat <- do.call(cbind, lapply(sim_list, `[[`, "conc_M1_P"))  # K x N_REP

P_rel_mat <- do.call(cbind, lapply(sim_list, `[[`, "conc_rel_P"))
Q_rel_mat <- do.call(cbind, lapply(sim_list, `[[`, "conc_rel_Q"))

gini_M0_vec   <- sapply(sim_list, `[[`, "gini_M0")
gini_M1_vec   <- sapply(sim_list, `[[`, "gini_M1")
pietra_M0_vec <- sapply(sim_list, `[[`, "pietra_M0")
pietra_M1_vec <- sapply(sim_list, `[[`, "pietra_M1")

nlpd_M0_vec <- sapply(sim_list, `[[`, "mean_nlpd_M0")
nlpd_M1_vec <- sapply(sim_list, `[[`, "mean_nlpd_M1")

ml_M0_vec   <- sapply(sim_list, `[[`, "mean_ml_M0")
ml_M1_vec   <- sapply(sim_list, `[[`, "mean_ml_M1")

# 2.5 Monte Carlo summary (analogous to Table 1)

summary_sim <- data.frame(
  Model = c("M0: AR(1)", "M1: MA(1)"),
  Pietra_index = c(mean(pietra_M0_vec), mean(pietra_M1_vec)),
  Gini_concentration_coefficient = c(mean(gini_M0_vec), mean(gini_M1_vec)),
  NLPD = c(mean(nlpd_M0_vec), mean(nlpd_M1_vec)),
  #Marginal_log_likelihood = c(mean(ml_M0_vec), mean(ml_M1_vec)),
  row.names = NULL
)

cat("\nMonte Carlo summary (averages over replications):\n")
print(summary_sim)

# 2.6 Agreement rates (analogous to Table 2)

M0_better_pietra <- pietra_M0_vec < pietra_M1_vec
M0_better_gini   <- gini_M0_vec   < gini_M1_vec
M0_better_nlpd   <- nlpd_M0_vec   < nlpd_M1_vec
#M0_better_ml     <- ml_M0_vec     > ml_M1_vec

agreement_rates <- data.frame(
  Criterion = c("Pietra index", "Gini concentration coefficient",
                "NLPD", "Marginal log-likelihood"),
  `M0 better than M1` = paste0(
    round(100 * c(mean(M0_better_pietra),
                  mean(M0_better_gini),
                  mean(M0_better_nlpd),
                  #mean(M0_better_ml)), 
          1),
    "%"
  )
)

cat("\nAgreement rates (M0 better than M1):\n")
print(agreement_rates, row.names = FALSE)
