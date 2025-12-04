## --- SIMULATION: concentration functions for AR(1) vs MA(1) ---

# x-axis for the uniform reference
x_unif <- (1:K) / K

# Mean concentration curves vs uniform
mean_P_M0 <- rowMeans(P_M0_mat)  # M0: AR(1)
mean_P_M1 <- rowMeans(P_M1_mat)  # M1: MA(1)

# 95% bands
ci_M0_low  <- apply(P_M0_mat, 1, quantile, probs = 0.025)
ci_M0_high <- apply(P_M0_mat, 1, quantile, probs = 0.975)

ci_M1_low  <- apply(P_M1_mat, 1, quantile, probs = 0.025)
ci_M1_high <- apply(P_M1_mat, 1, quantile, probs = 0.975)

par(mfrow = c(1, 2))

## Panel (a): concentration vs uniform for both models
plot(x_unif, mean_P_M0, type = "l", lwd = 2, col = "blue",
     ylim = c(0, 1), xlab = "Cumulative share of forecast points",
     ylab = "Cumulative predictive mass")
abline(0, 1, lwd = 2, lty = 2)  # equality line

polygon(c(x_unif, rev(x_unif)),
        c(ci_M0_low, rev(ci_M0_high)),
        border = NA, col = rgb(0, 0, 1, 0.15))
lines(x_unif, mean_P_M0, lwd = 2, col = "blue")

polygon(c(x_unif, rev(x_unif)),
        c(ci_M1_low, rev(ci_M1_high)),
        border = NA, col = rgb(1, 0, 0, 0.15))
lines(x_unif, mean_P_M1, lwd = 2, col = "red")

legend("topleft",
       legend = c("M0: AR(1)", "M1: MA(1)", "Equality line"),
       col = c("blue", "red", "black"),
       lty = c(1, 1, 2), lwd = 2, bty = "n")

## Panel (b): relative concentration of M1 relative to M0
mean_P_rel <- rowMeans(P_rel_mat)
mean_Q_rel <- rowMeans(Q_rel_mat)

plot(mean_Q_rel, mean_P_rel, type = "l", lwd = 2, col = "darkorange",
     xlab = "Cumulative mass under M0: AR(1)",
     ylab = "Cumulative mass under M1: MA(1)",
     xlim = c(0,1), ylim = c(0,1))
abline(0, 1, lwd = 2, lty = 2)
legend("topleft",
       legend = c("Concentration M1 vs M0", "Equality line"),
       col = c("darkorange", "black"),
       lty = c(1, 2), lwd = 2, bty = "n")

par(mfrow = c(1, 1))

## --- SIMULATION: histograms of Pietra and Gini indices ---

par(mfrow = c(1, 2))

# Pietra index
hist(pietra_M0_vec, breaks = 10, col = rgb(0, 0, 1, 0.4),
     main = "Pietra index (simulation)",
     xlab = "Pietra", xlim = range(c(pietra_M0_vec, pietra_M1_vec)))
hist(pietra_M1_vec, breaks = 10, col = rgb(1, 0, 0, 0.4), add = TRUE)
legend("topright",
       legend = c("M0: AR(1)", "M1: MA(1)"),
       fill = c(rgb(0, 0, 1, 0.4), rgb(1, 0, 0, 0.4)),
       bty = "n")

# Gini coefficient
hist(gini_M0_vec, breaks = 10, col = rgb(0, 0, 1, 0.4),
     main = "Gini coefficient (simulation)",
     xlab = "Gini", xlim = range(c(gini_M0_vec, gini_M1_vec)))
hist(gini_M1_vec, breaks = 10, col = rgb(1, 0, 0, 0.4), add = TRUE)
legend("topright",
       legend = c("M0: AR(1)", "M1: MA(1)"),
       fill = c(rgb(0, 0, 1, 0.4), rgb(1, 0, 0, 0.4)),
       bty = "n")

par(mfrow = c(1, 1))

#--------------

## --- REAL DATA: concentration functions vs equality line ---

# candidate_models: list with labels and (p,q) as defined in the script
model_labels <- vapply(candidate_models, `[[`, character(1), "label")
plot_indices <- 2:length(candidate_models)  # M1,...,M6

par(mfrow = c(2, 3))
for (j in plot_indices) {
  conc_M0 <- concentration_unif_list[[1]]   # M0: AR(1)
  conc_Mj <- concentration_unif_list[[j]]   # modelo alternativo
  
  plot(conc_M0$Q, conc_M0$P, type = "l", lwd = 2, col = "blue",
       xlab = "Cumulative share of forecast points",
       ylab = "Cumulative predictive mass",
       main = paste0("M0: AR(1) vs ", model_labels[j]),
       xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, lty = 2, lwd = 2)
  lines(conc_Mj$Q, conc_Mj$P, lwd = 2, col = "red")
  
  legend("topleft",
         legend = c("M0: AR(1)", model_labels[j], "Equality line"),
         col = c("blue", "red", "black"),
         lty = c(1, 1, 2), lwd = 2, bty = "n")
}
par(mfrow = c(1, 1))



## --- REAL DATA: relative concentration (each model vs M0) ---

p_M0 <- dens_mat[, 1] / sum(dens_mat[, 1])  # reference distribution (M0)

par(mfrow = c(2, 3))

for (j in plot_indices) {
  p_j <- dens_mat[, j]
  conc_rel_j <- build_concentration_relative(p_j, p_M0)
  
  plot(conc_rel_j$Q, conc_rel_j$P, type = "l", lwd = 2, col = "darkorange",
       xlab = "Cumulative mass under M0: AR(1)",
       ylab = "Cumulative mass under model",
       main = paste0("Relative: ", model_labels[j]),
       xlim = c(0, 1), ylim = c(0, 1))
  abline(0, 1, lty = 2, lwd = 2)
}

par(mfrow = c(1, 1))


#--------------------
summary_sim #Table 1
agreement_rates #Table 2
real_summary #Table 3
relative_summary #Table 4
colMeans(rbind(gini_M0_vec, gini_M1_vec))
colMeans(rbind(pietra_M0_vec, pietra_M1_vec))
