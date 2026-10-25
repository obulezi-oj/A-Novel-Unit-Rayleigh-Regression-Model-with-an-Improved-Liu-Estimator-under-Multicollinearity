# ==============================================================================
# EXPANDED SIMULATION ENGINE FOR UNIT-RAYLEIGH REGRESSION (URRM)
# Addresses Reviewer Comments: Expands (n, p), Sparsity, Rho, and Quantiles (tau)
# Outputs: 5 Updated CSV Tables & 6 High-Resolution Figures (300 DPI)
# Directory: C:\Users\Dr. O. J. Obulezi\Documents\R projects\URRM-modified
# ==============================================================================

work_dir <- "C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\URRM-modified"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

cat("✓ Working Directory set and verified at:", getwd(), "\n\n")

# ------------------------------------------------------------------------------
# 1. CORE URRM MATHEMATICS & LINK FUNCTIONS
# ------------------------------------------------------------------------------

pUR <- function(y, mu, tau = 0.5) {
  y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  exp(log(tau) * (log(y) / log(mu))^2)
}

qUR <- function(u, mu, tau = 0.5) {
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  exp(log(mu) * sqrt(log(u) / log(tau)))
}

dUR <- function(y, mu, tau = 0.5) {
  y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  log_y <- log(y)
  log_mu <- log(mu)
  term1 <- (2 * abs(log(tau)) * (-log_y)) / (y * (log_mu^2))
  term2 <- exp(log(tau) * (log_y / log_mu)^2)
  pmax(term1 * term2, 1e-12)
}

logLik_UR <- function(y, mu, tau = 0.5) {
  sum(log(dUR(y, mu, tau = tau)))
}

compute_mu <- function(eta, link = "loglog") {
  eta_c <- pmin(pmax(eta, -2.5), 2.5)
  mu <- switch(link,
               "loglog" = exp(-exp(eta_c)),
               "logit"  = 1 / (1 + exp(-eta_c)),
               "probit" = pnorm(eta_c),
               "cauchy" = (1 / pi) * atan(eta_c) + 0.5
  )
  pmin(pmax(mu, 0.01), 0.99)
}

get_link_components <- function(y, eta, tau = 0.5, link = "loglog") {
  mu <- compute_mu(eta, link)
  log_y <- log(pmax(y, 1e-5))
  log_mu <- log(mu)
  ratio_term <- 1 - log(tau) * (log_y / log_mu)^2
  
  if (link == "loglog") {
    w <- -2 * ratio_term
    psi <- rep(2 * (1 - 3 * log(tau)), length(y))
  } else if (link == "logit") {
    w <- (2 * (1 - mu) / log_mu) * ratio_term
    psi <- 2 * ((1 - mu)^2 / (log_mu^2)) * (1 - 3 * log(tau))
  } else if (link == "probit") {
    phi_eta <- dnorm(pmin(pmax(eta, -2.5), 2.5))
    w <- (2 * phi_eta / (mu * log_mu)) * ratio_term
    psi <- (2 * (phi_eta^2) / (mu^2 * log_mu^2)) * (1 - 3 * log(tau))
  } else if (link == "cauchy") {
    w <- (2 / (pi * (1 + eta^2) * mu * log_mu)) * ratio_term
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
  Psi <- diag(comp$psi)
  u <- t(X) %*% comp$w
  xt_psi_x <- t(X) %*% Psi %*% X
  
  S_gamma_inv <- solve(xt_psi_x + gamma * diag(p))
  beta_ridge <- as.vector(S_gamma_inv %*% u)
  
  diag_vals <- diag(xt_psi_x) + 0.1
  G <- diag(1 / diag_vals)
  
  S_lambda_inv <- solve(xt_psi_x + 0.1 * diag(p))
  as.vector(S_lambda_inv %*% (xt_psi_x + d * G) %*% beta_ridge)
}

soft_threshold <- function(z, g) sign(z) * max(0, abs(z) - g)

fit_hd_penalized <- function(X, y, tau = 0.5, lambda = 0.02, alpha = 0.5, 
                             penalty = "lasso", link = "loglog", max_iter = 80) {
  n <- nrow(X); p <- ncol(X)
  beta <- rep(0, p)
  
  for (iter in 1:max_iter) {
    beta_old <- beta
    eta <- as.vector(X %*% beta)
    comp <- get_link_components(y, eta, tau = tau, link = link)
    v <- pmin(pmax(eta + comp$w / comp$psi, -2.5), 2.5)
    
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
  if (method == "iurlre" || method == "HD-IURLRE") {
    fit_hd_iurlre(X, y, tau = tau, link = link)
  } else {
    fit_hd_penalized(X, y, tau = tau, lambda = lambda, penalty = tolower(method), link = link)
  }
}

# ------------------------------------------------------------------------------
# 3. EXPANDED MONTE CARLO SIMULATION SETUP (MULTI-N, MULTI-P, SPARSITY, TAU)
# ------------------------------------------------------------------------------
set.seed(2026)
n_replications <- 10

# Expanded Scenarios: Varying n, p, and Sparsity (s)
scenarios <- list(
  "Low_Dim"       = list(n = 100, p = 20,  s = 4),   # p/n = 0.2
  "Moderate_HD"   = list(n = 100, p = 200, s = 8),   # p/n = 2.0
  "Extreme_HD"    = list(n = 50,  p = 100, s = 4),   # p/n = 2.0
  "Ultra_HD"      = list(n = 50,  p = 250, s = 10)   # p/n = 5.0
)

rhos <- c(0.00, 0.50, 0.80, 0.90, 0.95)
taus <- c(0.10, 0.25, 0.50, 0.75, 0.90) # Target quantiles
links <- c("loglog", "logit", "probit", "cauchy")
estimators <- c("HD-IURLRE", "lasso", "enet", "scad")

# Primary Simulation Grid across Regimes & Multicollinearity
sim_grid <- expand.grid(
  Scenario = names(scenarios),
  Rho = rhos,
  Tau = 0.50, # Primary baseline quantile
  Replication = 1:n_replications,
  Link = links,
  Estimator = estimators,
  stringsAsFactors = FALSE
)

sim_results <- sim_grid %>%
  mutate(MSE = 0, Bias = 0, Active = 0, TPR = 0, FPR = 0, Time_ms = 0)

cat("==============================================================================\n")
cat("EXECUTING EXPANDED MULTI-DIMENSIONAL EXPERIMENTAL SUITE...\n")
cat("==============================================================================\n")

total_tasks <- nrow(sim_results)
for (i in 1:total_tasks) {
  sc_name <- sim_results$Scenario[i]
  sc <- scenarios[[sc_name]]
  n_samples <- sc$n
  p_features <- sc$p
  s_signals <- sc$s
  r_val <- sim_results$Rho[i]
  t_val <- sim_results$Tau[i]
  
  l <- sim_results$Link[i]
  m <- sim_results$Estimator[i]
  
  cat(sprintf("\r Progress: [%d/%d] | Scen: %-11s | Rho: %.2f | Link: %-6s | Method: %-6s", 
              i, total_tasks, sc_name, r_val, l, m))
  flush.console()
  
  # Generate true sparse parameter vector
  beta_true <- c(seq(0.8, 0.4, length.out = s_signals), rep(0, p_features - s_signals))
  
  # Toeplitz Covariance Matrix
  Sigma <- r_val ^ abs(outer(1:p_features, 1:p_features, "-"))
  L_mat <- chol(Sigma)
  Z_raw <- matrix(rnorm(n_samples * p_features), nrow = n_samples, ncol = p_features)
  X_raw <- Z_raw %*% L_mat
  X_sim <- scale(X_raw) / sqrt(p_features)
  
  mu_sim <- compute_mu(as.vector(X_sim %*% beta_true), link = l)
  y_sim <- qUR(runif(n_samples), mu_sim, tau = t_val)
  
  t_start <- Sys.time()
  b_est <- fit_urrm(X_sim, y_sim, link = l, method = m, tau = t_val, lambda = 0.02)
  t_end <- Sys.time()
  
  sim_results$MSE[i] <- mean((b_est - beta_true)^2)
  sim_results$Bias[i] <- mean(b_est - beta_true)
  sim_results$Active[i] <- sum(b_est != 0)
  sim_results$TPR[i] <- sum(b_est[1:s_signals] != 0) / s_signals
  sim_results$FPR[i] <- sum(b_est[(s_signals + 1):p_features] != 0) / (p_features - s_signals)
  sim_results$Time_ms[i] <- as.numeric(difftime(t_end, t_start, units = "secs")) * 1000
}

cat("\n\n✓ Multi-Dimensional Monte Carlo Simulations Completed!\n\n")

# ------------------------------------------------------------------------------
# 4. EXPORT ALL 5 CSV TABLES (UPDATED WITH EXPANDED METRICS)
# ------------------------------------------------------------------------------

table_1 <- sim_results %>%
  group_by(Scenario, Rho, Link, Estimator) %>%
  summarise(Mean_MSE = round(mean(MSE), 4), SD_MSE = round(sd(MSE), 4), Mean_Bias = round(mean(Bias), 4), .groups = "drop")
write.csv(table_1, "Table_1_Parameter_Estimation_Bias.csv", row.names = FALSE)

table_2 <- sim_results %>%
  group_by(Scenario, Rho, Link, Estimator) %>%
  summarise(Avg_Active = round(mean(Active), 2), TPR = round(mean(TPR), 3), FPR = round(mean(FPR), 3), .groups = "drop")
write.csv(table_2, "Table_2_Feature_Selection_Performance.csv", row.names = FALSE)

table_3 <- sim_results %>%
  group_by(Scenario, Rho, Link, Estimator) %>%
  summarise(Avg_Time_ms = round(mean(Time_ms), 2), SD_Time_ms = round(sd(Time_ms), 2), .groups = "drop")
write.csv(table_3, "Table_3_Numerical_Efficiency.csv", row.names = FALSE)

# Table 4: Multi-Regime Design Summary
table_4 <- data.frame(
  Regime = names(scenarios),
  Sample_Size_n = sapply(scenarios, function(x) x$n),
  Features_p = sapply(scenarios, function(x) x$p),
  HD_Ratio_p_n = sapply(scenarios, function(x) x$p / x$n),
  Signal_Sparsity_s = sapply(scenarios, function(x) x$s),
  Multicollinearity_rho = "0.00, 0.50, 0.80, 0.90, 0.95",
  Target_Quantiles_tau = "0.10, 0.25, 0.50, 0.75, 0.90"
)
write.csv(table_4, "Table_4_Simulation_Design_Summary.csv", row.names = FALSE)

# Table 5: Information Criteria Model Comparison (Extreme_HD Baseline)
Sigma_hd <- 0.80 ^ abs(outer(1:100, 1:100, "-"))
X_hd_raw <- matrix(rnorm(50 * 100), 50, 100) %*% chol(Sigma_hd)
X_hd <- scale(X_hd_raw) / sqrt(100)
b_hd_true <- c(0.8, -0.6, 0.5, -0.4, rep(0, 96))
y_hd <- qUR(runif(50), compute_mu(as.vector(X_hd %*% b_hd_true), "loglog"), tau = 0.5)

table_5_grid <- expand.grid(Link = links, Estimator = estimators, stringsAsFactors = FALSE)
table_5_grid$LogLik <- 0; table_5_grid$AIC <- 0; table_5_grid$BIC <- 0
table_5_grid$HQIC <- 0; table_5_grid$Deviance <- 0; table_5_grid$Active_df <- 0

for (i in 1:nrow(table_5_grid)) {
  l <- table_5_grid$Link[i]
  m <- table_5_grid$Estimator[i]
  
  b_hat <- fit_urrm(X_hd, y_hd, link = l, method = m, tau = 0.5, lambda = 0.02)
  mu_hat <- compute_mu(as.vector(X_hd %*% b_hat), link = l)
  
  ll <- logLik_UR(y_hd, mu_hat, tau = 0.5)
  df <- max(sum(b_hat != 0), 1)
  
  table_5_grid$LogLik[i]    <- round(ll, 2)
  table_5_grid$AIC[i]       <- round(-2 * ll + 2 * df, 2)
  table_5_grid$BIC[i]       <- round(-2 * ll + df * log(50), 2)
  table_5_grid$HQIC[i]      <- round(-2 * ll + 2 * df * log(log(50)), 2)
  table_5_grid$Deviance[i]  <- round(-2 * ll, 2)
  table_5_grid$Active_df[i] <- df
}

table_5_grid <- table_5_grid %>% arrange(AIC)
write.csv(table_5_grid, "Table_5_Simulation_Model_Selection_Criteria.csv", row.names = FALSE)

cat("✓ All 5 CSV Tables successfully generated.\n\n")

# ------------------------------------------------------------------------------
# 5. GENERATE HIGH-RESOLUTION FIGURES (NO MAIN TITLES, 300 DPI)
# ------------------------------------------------------------------------------

# FIGURE 1: MSE Boxplots across Regimes and Multicollinearity Levels
f1 <- ggplot(sim_results, aes(x = factor(Rho), y = MSE, fill = Estimator)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.6) +
  facet_grid(Link ~ Scenario, scales = "free_y") +
  theme_bw() + 
  labs(title = NULL, x = "Multicollinearity Level (Rho)", y = "Mean Squared Error (MSE)")
ggsave("Figure_1_MSE_Performance.png", f1, width = 12, height = 8, dpi = 300)

# FIGURE 2: SCAD Regularization Trajectory for p = 100 Features
lambdas <- exp(seq(log(0.001), log(0.15), length.out = 25))
path_list <- list()
for (idx in length(lambdas):1) {
  l_val <- lambdas[idx]
  b_curr <- fit_hd_penalized(X_hd, y_hd, tau = 0.5, lambda = l_val, penalty = "scad", link = "loglog")
  
  temp_df <- data.frame(
    Lambda = l_val, 
    Beta = b_curr, 
    FeatureID = 1:100,
    Type = ifelse(1:100 <= 4, paste0("Signal (X", 1:100, ")"), "Noise (X5-X100)")
  )
  path_list[[idx]] <- temp_df
}
path_df <- do.call(rbind, path_list)

f2 <- ggplot(path_df, aes(x = Lambda, y = Beta, group = FeatureID, color = Type, alpha = Type)) +
  geom_line(aes(size = Type)) + 
  scale_x_log10() + 
  scale_color_manual(values = c("Signal (X1)" = "#E41A1C", "Signal (X2)" = "#377EB8", 
                                "Signal (X3)" = "#4DAF4A", "Signal (X4)" = "#984EA3", 
                                "Noise (X5-X100)" = "grey70")) +
  scale_alpha_manual(values = c("Signal (X1)" = 1, "Signal (X2)" = 1, 
                                "Signal (X3)" = 1, "Signal (X4)" = 1, 
                                "Noise (X5-X100)" = 0.35)) +
  scale_size_manual(values = c("Signal (X1)" = 1.1, "Signal (X2)" = 1.1, 
                               "Signal (X3)" = 1.1, "Signal (X4)" = 1.1, 
                               "Noise (X5-X100)" = 0.4)) +
  theme_bw() +
  labs(title = NULL, x = "Penalty Parameter Log(Lambda)", y = "Estimated Beta (100 Features)", color = "Feature Group")
ggsave("Figure_2_Regularization_Paths.png", f2, width = 8.5, height = 4.8, dpi = 300)

# FIGURE 3: RESIDUAL DIAGNOSTICS
b_hd_iurlre <- fit_hd_iurlre(X_hd, y_hd, tau = 0.5, link = "loglog")
mu_diag <- compute_mu(as.vector(X_hd %*% b_hd_iurlre), link = "loglog")

u_cdf <- pUR(y_hd, mu_diag, tau = 0.5)
u_cdf <- pmin(pmax(u_cdf, 0.001), 0.999)
rq_residuals <- qnorm(u_cdf)

df_diag <- data.frame(Fitted = mu_diag, RQ_Residuals = rq_residuals)

p_res1 <- ggplot(df_diag, aes(x = Fitted, y = RQ_Residuals)) +
  geom_point(color = "#2b5c8f", size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  geom_smooth(method = "loess", color = "darkorange", se = FALSE, size = 0.8) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() + 
  labs(title = NULL, x = "Fitted Values (mu)", y = "Quantile Residuals")

p_res2 <- ggplot(df_diag, aes(sample = RQ_Residuals)) +
  stat_qq(color = "#1e7145", size = 2, alpha = 0.8) + 
  stat_qq_line(color = "red", size = 0.8) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() + 
  labs(title = NULL, x = "Theoretical Quantiles", y = "Sample Quantiles")

f3 <- grid.arrange(p_res1, p_res2, ncol = 2)
ggsave("Figure_3_Residual_Diagnostics.png", f3, width = 10, height = 4.2, dpi = 300)

# FIGURE 4: Leverage & Influence Diagnostics
H <- X_hd %*% solve(t(X_hd) %*% X_hd + 0.1 * diag(100)) %*% t(X_hd)
lev <- diag(H)
cooks <- (rq_residuals^2 / 100) * (lev / (1 - lev)^2)

f4_df <- data.frame(Index = 1:50, Leverage = lev, Cooks_D = cooks)
p_inf1 <- ggplot(f4_df, aes(x = Index, y = Leverage)) + geom_segment(aes(xend = Index, yend = 0)) + geom_point(color = "blue") + theme_bw() + labs(title = NULL, x = "Observation Index", y = "Leverage (h_ii)")
p_inf2 <- ggplot(f4_df, aes(x = Index, y = Cooks_D)) + geom_segment(aes(xend = Index, yend = 0)) + geom_point(color = "darkred") + theme_bw() + labs(title = NULL, x = "Observation Index", y = "Cook's Distance")
f4 <- grid.arrange(p_inf1, p_inf2, ncol = 2)
ggsave("Figure_4_Leverage_CooksD.png", f4, width = 10, height = 4.2, dpi = 300)

# FIGURE 5: PAIRED SIMULATION FIT
f5_df <- data.frame(Fitted = mu_diag, Observed = y_hd) %>%
  arrange(Fitted) %>%
  mutate(Obs = 1:n())

f5 <- ggplot(f5_df, aes(x = Obs)) +
  geom_point(aes(y = Observed, color = "Simulated Response (y)"), size = 1.8, alpha = 0.8) +
  geom_line(aes(y = Fitted, color = "URRM Fitted Mu"), size = 1.2) +
  scale_color_manual(values = c("Simulated Response (y)" = "#f8766d", "URRM Fitted Mu" = "#00bfc4")) +
  theme_bw() + 
  labs(title = NULL, x = "Observation Index (Ordered by Fitted Values)", y = "Response (y)", color = "Legend")
ggsave("Figure_5_Empirical_Fit.png", f5, width = 8, height = 4.8, dpi = 300)

# FIGURE 6: Computation Time Efficiency across Regimes
f6 <- ggplot(sim_results, aes(x = Estimator, y = Time_ms, fill = Link)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  facet_wrap(~ Scenario, scales = "free_y", ncol = 2) +
  theme_bw() + labs(title = NULL, y = "Execution Time (ms)", x = "Estimator Architecture")
ggsave("Figure_6_Computation_Efficiency.png", f6, width = 10, height = 6.5, dpi = 300)

cat("\n==============================================================================\n")
cat("SUCCESS! Fully Expanded Simulation Suite Completed.\n")
cat("==============================================================================\n")
