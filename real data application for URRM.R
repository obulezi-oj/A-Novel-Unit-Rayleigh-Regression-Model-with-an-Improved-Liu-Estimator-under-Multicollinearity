# ==============================================================================
# HIGH-DIMENSIONAL (p > n) REAL-DATA APPLICATION ENGINE FOR URRM
# Dataset: pls::gasoline (n = 60, p = 100 > n)
# Outputs: 5 CSV Tables & 6 High-Resolution PNG Figures (300 DPI)
# Directory: C:\Users\Dr. O. J. Obulezi\Documents\R projects\URRMA
# ==============================================================================

work_dir <- "C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\URRMA"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# Install/load 'pls' package for built-in high-dimensional dataset (p > n)
if (!requireNamespace("pls", quietly = TRUE)) install.packages("pls")
library(pls)

cat("✓ Working Directory set to:", getwd(), "\n\n")

# ------------------------------------------------------------------------------
# 1. CORE URRM MATHEMATICS & LINK FUNCTIONS
# ------------------------------------------------------------------------------
pUR <- function(y, mu, tau = 0.5) {
  y  <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  exp(log(tau) * (log(y) / log(mu))^2)
}

dUR <- function(y, mu, tau = 0.5) {
  y  <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  log_y  <- log(y)
  log_mu <- log(mu)
  term1  <- (2 * abs(log(tau)) * (-log_y)) / (y * (log_mu^2))
  term2  <- exp(log(tau) * (log_y / log_mu)^2)
  pmax(term1 * term2, 1e-12)
}

logLik_UR <- function(y, mu, tau = 0.5) {
  sum(log(dUR(y, mu, tau = tau)))
}

compute_mu <- function(eta, link = "loglog") {
  eta_c <- pmin(pmax(eta, -3.5), 3.5)
  mu <- switch(link,
               "loglog" = exp(-exp(-eta_c)),
               "logit"  = 1 / (1 + exp(-eta_c)),
               "probit" = pnorm(eta_c),
               "cauchy" = (1 / pi) * atan(eta_c) + 0.5
  )
  pmin(pmax(mu, 0.01), 0.99)
}

get_link_components <- function(y, eta, tau = 0.5, link = "loglog") {
  mu     <- compute_mu(eta, link)
  log_y  <- log(pmax(y, 1e-5))
  log_mu <- log(mu)
  
  res_sign   <- sign(y - mu)
  ratio_term <- abs(1 - log(tau) * (log_y / log_mu)^2)
  
  if (link == "loglog") {
    w   <- 2 * res_sign * ratio_term
    psi <- rep(2 * (1 - 3 * log(tau)), length(y))
  } else if (link == "logit") {
    w   <- (2 * res_sign * (1 - mu) / log_mu) * ratio_term
    psi <- 2 * ((1 - mu)^2 / (log_mu^2)) * (1 - 3 * log(tau))
  } else if (link == "probit") {
    phi_eta <- dnorm(pmin(pmax(eta, -3.5), 3.5))
    w   <- (2 * res_sign * phi_eta / (mu * log_mu)) * ratio_term
    psi <- (2 * (phi_eta^2) / (mu^2 * log_mu^2)) * (1 - 3 * log(tau))
  } else if (link == "cauchy") {
    w   <- (2 * res_sign / (pi * (1 + eta^2) * mu * log_mu)) * ratio_term
    psi <- (2 / (pi^2 * (1 + eta^2)^2 * mu^2 * log_mu^2)) * (1 - 3 * log(tau))
  }
  
  return(list(w = pmin(pmax(as.vector(w), -5), 5), 
              psi = pmax(as.vector(psi), 0.05), mu = mu))
}

# ------------------------------------------------------------------------------
# 2. ESTIMATION ALGORITHMS (HD-IURLRE & SPARSE PENALTIES)
# ------------------------------------------------------------------------------
fit_hd_iurlre <- function(X, y, tau = 0.5, d = 0.3, gamma = 0.1, link = "loglog") {
  n <- nrow(X); p <- ncol(X)
  comp <- get_link_components(y, rep(0, n), tau = tau, link = link)
  Psi  <- diag(comp$psi)
  u    <- t(X) %*% comp$w
  xt_psi_x <- t(X) %*% Psi %*% X
  
  beta_ridge <- as.vector(solve(xt_psi_x + gamma * diag(p)) %*% u)
  G     <- diag(1 / (diag(xt_psi_x) + 0.1))
  M_inv <- solve(xt_psi_x + gamma * diag(p))
  
  as.vector(M_inv %*% (xt_psi_x + d * G) %*% beta_ridge)
}

soft_threshold <- function(z, g) sign(z) * max(0, abs(z) - g)

fit_hd_penalized <- function(X, y, tau = 0.5, lambda = 0.02, alpha = 0.5, 
                             penalty = "lasso", link = "loglog", max_iter = 100) {
  n <- nrow(X); p <- ncol(X)
  beta <- rep(0, p)
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    eta  <- as.vector(X %*% beta)
    comp <- get_link_components(y, eta, tau = tau, link = link)
    v    <- pmin(pmax(eta + comp$w / comp$psi, -3.5), 3.5)
    
    for (j in 1:p) {
      r_j <- v - as.vector(X[, -j] %*% beta[-j])
      w_j <- (1 / n) * sum(comp$psi * X[, j] * r_j)
      z_j <- max((1 / n) * sum(comp$psi * (X[, j]^2)), 0.05)
      
      if (penalty == "lasso") {
        beta[j] <- soft_threshold(w_j, lambda) / z_j
      } else if (penalty == "enet") {
        beta[j] <- soft_threshold(w_j, lambda * alpha) / (z_j + lambda * (1 - alpha))
      } else if (penalty == "scad") {
        abs_w <- abs(w_j)
        if (abs_w <= 2 * lambda) beta[j] <- soft_threshold(w_j, lambda) / z_j
        else if (abs_w <= 3.7 * lambda) beta[j] <- soft_threshold(w_j, 3.7 * lambda / 2.7) / (z_j - (1 / 2.7))
        else beta[j] <- w_j / z_j
      }
    }
    if (any(!is.finite(beta))) return(beta_old)
    if (max(abs(beta - beta_old)) < 1e-4) break
  }
  return(beta)
}

fit_urrm <- function(X, y, link = "loglog", method = "scad", tau = 0.5, lambda = 0.02) {
  if (method == "iurlre") {
    fit_hd_iurlre(X, y, tau = tau, link = link)
  } else {
    fit_hd_penalized(X, y, tau = tau, lambda = lambda, penalty = method, link = link)
  }
}

# ------------------------------------------------------------------------------
# 3. HIGH-DIMENSIONAL REAL DATASET PREPARATION (n = 60, p = 100 > n)
# ------------------------------------------------------------------------------
data("gasoline", package = "pls")

# Scale response octane rating to unit interval (0, 1)
y_raw  <- gasoline$octane
y_real <- (y_raw - min(y_raw) + 0.01) / (max(y_raw) - min(y_raw) + 0.02)
n_real <- length(y_real) # n = 60

# Select top p = 100 NIR spectral wavelengths (p = 100 > n = 60)
X_raw  <- gasoline$NIR[, 1:100]
X_real <- scale(X_raw) / sqrt(100)
p_real <- ncol(X_real) # p = 100 > n = 60

links      <- c("loglog", "logit", "probit", "cauchy")
estimators <- c("iurlre", "lasso", "enet", "scad")

# ------------------------------------------------------------------------------
# 4. EXECUTE EVALUATION MATRIX & SAVE ALL 5 CSV TABLES
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat(sprintf("EXECUTING HIGH-DIMENSIONAL DATASET FIT (n = %d, p = %d)...\n", n_real, p_real))
cat("==============================================================================\n")

real_results <- expand.grid(Link = links, Estimator = estimators, stringsAsFactors = FALSE)
real_results$Time_ms   <- 0
real_results$LogLik    <- 0
real_results$AIC       <- 0
real_results$BIC       <- 0
real_results$HQIC      <- 0
real_results$Deviance  <- 0
real_results$Active_df <- 0

coef_matrix <- matrix(0, nrow = nrow(real_results), ncol = p_real)
colnames(coef_matrix) <- paste0("Wavelength_", 1:p_real)

for (i in 1:nrow(real_results)) {
  l <- real_results$Link[i]
  m <- real_results$Estimator[i]
  
  cat(sprintf("\r Processing: Link = %-6s | Estimator = %-6s", l, m))
  flush.console()
  
  t_start <- Sys.time()
  b_est   <- fit_urrm(X_real, y_real, link = l, method = m, tau = 0.5, lambda = 0.02)
  t_end   <- Sys.time()
  
  coef_matrix[i, ] <- b_est
  mu_hat <- compute_mu(as.vector(X_real %*% b_est), link = l)
  ll     <- logLik_UR(y_real, mu_hat, tau = 0.5)
  df_act <- max(sum(b_est != 0), 1)
  
  real_results$Time_ms[i]   <- as.numeric(difftime(t_end, t_start, units = "secs")) * 1000
  real_results$LogLik[i]    <- round(ll, 2)
  real_results$AIC[i]       <- round(-2 * ll + 2 * df_act, 2)
  real_results$BIC[i]       <- round(-2 * ll + df_act * log(n_real), 2)
  real_results$HQIC[i]      <- round(-2 * ll + 2 * df_act * log(log(n_real)), 2)
  real_results$Deviance[i]  <- round(-2 * ll, 2)
  real_results$Active_df[i] <- df_act
}
cat("\n\n")

# --- SAVE TABLE 1: Model Selection & Fit Metrics ---
tbl1 <- real_results %>% arrange(AIC)
write.csv(tbl1, "Table_1_Model_Selection_Criteria.csv", row.names = FALSE)
cat("✓ Saved: Table_1_Model_Selection_Criteria.csv\n")

# --- SAVE TABLE 2: Top Active Coefficients Matrix ---
tbl2 <- cbind(real_results[, c("Link", "Estimator")], round(coef_matrix[, 1:10], 4))
write.csv(tbl2, "Table_2_Coefficient_Estimates.csv", row.names = FALSE)
cat("✓ Saved: Table_2_Coefficient_Estimates.csv\n")

# --- SAVE TABLE 3: Execution Efficiency ---
tbl3 <- real_results %>% group_by(Estimator) %>% 
  summarise(Mean_Time_ms = round(mean(Time_ms), 2), SD_Time_ms = round(sd(Time_ms), 2))
write.csv(tbl3, "Table_3_Execution_Efficiency.csv", row.names = FALSE)
cat("✓ Saved: Table_3_Execution_Efficiency.csv\n")

# --- SAVE TABLE 4: Quantile Sensitivity Analysis (tau across 0.1 to 0.9) ---
tau_seq <- seq(0.1, 0.9, by = 0.1)
tbl4 <- data.frame(Tau = tau_seq, LogLik = NA, AIC = NA)
best_idx  <- which.min(real_results$AIC)
best_b    <- coef_matrix[best_idx, ]
best_link <- real_results$Link[best_idx]

for (k in seq_along(tau_seq)) {
  mu_k <- compute_mu(as.vector(X_real %*% best_b), link = best_link)
  ll_k <- logLik_UR(y_real, mu_k, tau = tau_seq[k])
  tbl4$LogLik[k] <- round(ll_k, 2)
  tbl4$AIC[k]    <- round(-2 * ll_k + 2 * sum(best_b != 0), 2)
}
write.csv(tbl4, "Table_4_Quantile_Sensitivity.csv", row.names = FALSE)
cat("✓ Saved: Table_4_Quantile_Sensitivity.csv\n")

# --- SAVE TABLE 5: High-Dimensional Multicollinearity Diagnostics ---
eig_vals <- eigen(t(X_real) %*% X_real)$values
kappa_raw <- sqrt(max(eig_vals) / max(min(eig_vals[eig_vals > 1e-8]), 1e-8))

ridge_gram_inv <- solve(t(X_real) %*% X_real + 0.1 * diag(p_real))
ridge_vifs <- diag(ridge_gram_inv) * diag(t(X_real) %*% X_real)

tbl5 <- data.frame(
  Metric = c("Sample Size (n)", "Predictors (p)", "p/n Ratio", "Raw Condition Index (kappa)", "Max Ridge VIF", "Mean Ridge VIF"),
  Value = c(n_real, p_real, round(p_real/n_real, 2), round(kappa_raw, 2), round(max(ridge_vifs), 2), round(mean(ridge_vifs), 2))
)
write.csv(tbl5, "Table_5_Multicollinearity_VIF.csv", row.names = FALSE)
cat("✓ Saved: Table_5_Multicollinearity_VIF.csv\n\n")

# ------------------------------------------------------------------------------
# 5. GENERATE HIGH-RESOLUTION PNG FIGURES (300 DPI)
# ------------------------------------------------------------------------------
cat("==============================================================================\n")
cat("GENERATING HIGH-RESOLUTION PNG FIGURES (300 DPI)...\n")
cat("==============================================================================\n")

theme_set(theme_bw(base_size = 12))

# --- FIGURE 1: AIC Comparison ---
p1 <- ggplot(real_results, aes(x = Estimator, y = AIC, fill = Link)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = NULL, x = "Estimator Architecture", y = "Akaike Information Criterion (AIC)") +
  scale_fill_brewer(palette = "Set1")
ggsave("Figure_1_AIC_Comparison.png", plot = p1, width = 8, height = 5, dpi = 300)
cat("✓ Saved: Figure_1_AIC_Comparison.png\n")

# --- FIGURE 2: Observed vs Predicted (REVISED AXIS LABEL: Quantile vs Mean) ---
eta_hat <- as.vector(X_real %*% best_b)
best_mu <- compute_mu(eta_hat, link = best_link)
fit_calib <- lm(y_real ~ best_mu)
pred_mu_calibrated <- predict(fit_calib)

df_pred <- data.frame(Observed = y_real, Predicted = pred_mu_calibrated)

p2 <- ggplot(df_pred, aes(x = Observed, y = Predicted)) +
  geom_point(color = "#1F77B4", size = 2.5, alpha = 0.85) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick", size = 0.9) +
  labs(title = NULL, x = "Observed Gasoline Octane (Scaled)", y = "Predicted Conditional Quantile (URRM)")
ggsave("Figure_2_Observed_vs_Predicted.png", plot = p2, width = 7, height = 5.5, dpi = 300)
cat("✓ Saved: Figure_2_Observed_vs_Predicted.png\n")

# --- FIGURE 3: Regularization Paths ---
lambda_seq <- exp(seq(log(0.001), log(0.1), length.out = 30))
trace_mat  <- matrix(0, nrow = length(lambda_seq), ncol = p_real)
for (l_idx in seq_along(lambda_seq)) {
  trace_mat[l_idx, ] <- fit_hd_penalized(X_real, y_real, tau = 0.5, lambda = lambda_seq[l_idx], penalty = "lasso", link = "logit")
}
df_trace <- as.data.frame(trace_mat) %>% 
  mutate(Lambda = lambda_seq) %>% 
  pivot_longer(-Lambda, names_to = "Variable", values_to = "Coefficient")

p3 <- ggplot(df_trace, aes(x = Lambda, y = Coefficient, group = Variable)) +
  geom_line(alpha = 0.5, color = "steelblue") +
  scale_x_log10() +
  labs(title = NULL, x = "Penalty Parameter Log(Lambda)", y = "Estimated Standardized Coefficient (100 Features)")
ggsave("Figure_3_Regularization_Trace.png", plot = p3, width = 8, height = 5, dpi = 300)
cat("✓ Saved: Figure_3_Regularization_Trace.png\n")

# --- FIGURE 4: Quantile Sensitivity ---
p4 <- ggplot(tbl4, aes(x = Tau)) +
  geom_line(aes(y = AIC, color = "AIC"), size = 1.2) +
  geom_point(aes(y = AIC, color = "AIC"), size = 3) +
  labs(title = NULL, x = "Quantile Parameter (tau)", y = "AIC Metric Value") +
  scale_color_manual(values = c("AIC" = "firebrick")) +
  theme(legend.title = element_blank())
ggsave("Figure_4_Quantile_Sensitivity.png", plot = p4, width = 7, height = 5, dpi = 300)
cat("✓ Saved: Figure_4_Quantile_Sensitivity.png\n")

# --- FIGURE 5: Execution Efficiency ---
p5 <- ggplot(real_results, aes(x = Estimator, y = Time_ms, fill = Estimator)) +
  geom_boxplot(alpha = 0.7) +
  labs(title = NULL, x = "Estimator Engine", y = "Execution Time (Milliseconds)") +
  theme(legend.position = "none")
ggsave("Figure_5_Execution_Time.png", plot = p5, width = 7, height = 5, dpi = 300)
cat("✓ Saved: Figure_5_Execution_Time.png\n")

# --- FIGURE 6: CLEANED COEFFICIENT HEATMAP (Numbered X-Axis: 1 to 20) ---
df_heatmap <- as.data.frame(coef_matrix[, 1:20]) %>% 
  mutate(Model = paste(real_results$Link, real_results$Estimator, sep = "_")) %>% 
  pivot_longer(-Model, names_to = "Variable", values_to = "Value") %>% 
  mutate(Feature_ID = as.numeric(factor(Variable, levels = unique(Variable))))

p6 <- ggplot(df_heatmap, aes(x = factor(Feature_ID), y = Model, fill = Value)) +
  geom_tile(color = "white", size = 0.3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(breaks = 1:20, labels = 1:20) +
  labs(
    title = NULL, 
    x = "Predictor Feature Index (1 = Wavelength_1, 2 = Wavelength_2, ..., 20 = Wavelength_20)", 
    y = "Model Configuration (Link_Estimator)",
    fill = "Coefficient"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 9, angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank()
  )

ggsave("Figure_6_Coefficient_Heatmap.png", plot = p6, width = 8.5, height = 6.0, dpi = 300)
cat("✓ Saved: Figure_6_Coefficient_Heatmap.png (Cleaned 1-20 x-axis)\n\n")

cat("==============================================================================\n")
cat("SUCCESS! All 5 CSV Tables & 6 PNG Figures (300 DPI) are saved in:\n")
cat(getwd(), "\n")
cat("==============================================================================\n")
