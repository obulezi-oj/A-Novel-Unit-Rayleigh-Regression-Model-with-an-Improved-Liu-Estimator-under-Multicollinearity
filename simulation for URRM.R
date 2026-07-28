# ==============================================================================
# UNIFIED SIMULATION PUBLICATION ENGINE FOR UNIT-RAYLEIGH REGRESSION (URRM)
# Outputs: 5 CSV Tables & 6 Corrected High-Resolution Figures (300 DPI)
# Directory: C:\Users\Dr. O. J. Obulezi\Documents\R projects\URRM
# ==============================================================================

work_dir <- "C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\URRM"
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

cat("✓ Working Directory set to:", getwd(), "\n\n")

# ------------------------------------------------------------------------------
# 1. CORE URRM MATHEMATICS & LINK FUNCTIONS
# ------------------------------------------------------------------------------

# Exact Unit-Rayleigh CDF: F(y; mu, tau) = exp( log(tau) * (log(y)/log(mu))^2 )
pUR <- function(y, mu, tau = 0.5) {
  y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  exp(log(tau) * (log(y) / log(mu))^2)
}

# Quantile Function
qUR <- function(u, mu, tau = 0.5) {
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  exp(log(mu) * sqrt(log(u) / log(tau)))
}

# Probability Density Function
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
# 2. ESTIMATION ALGORITHMS
# ------------------------------------------------------------------------------

fit_hd_iurlre <- function(X, y, tau = 0.5, d = 0.3, gamma = 0.1, link = "loglog") {
  n <- nrow(X); p <- ncol(X)
  comp <- get_link_components(y, rep(0, n), tau = tau, link = link)
  Psi <- diag(comp$psi)
  u <- t(X) %*% comp$w
  xt_psi_x <- t(X) %*% Psi %*% X
  
  beta_ridge <- as.vector(solve(xt_psi_x + gamma * diag(p)) %*% u)
  G <- diag(1 / (diag(xt_psi_x) + 0.1))
  M_inv <- solve(xt_psi_x + gamma * diag(p))
  
  as.vector(M_inv %*% (xt_psi_x + d * G) %*% beta_ridge)
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
  if (method == "iurlre") {
    fit_hd_iurlre(X, y, tau = tau, link = link)
  } else {
    fit_hd_penalized(X, y, tau = tau, lambda = lambda, penalty = method, link = link)
  }
}

# ------------------------------------------------------------------------------
# 3. MONTE CARLO SIMULATION (GENERATING TABLES 1, 2, 3)
# ------------------------------------------------------------------------------
set.seed(2026)
n_samples <- 100
p_features <- 20
tau_level <- 0.5
n_replications <- 10

beta_true <- c(0.8, -0.6, 0.5, -0.4, rep(0, p_features - 4))
links <- c("loglog", "logit", "probit", "cauchy")
estimators <- c("iurlre", "lasso", "enet", "scad")

sim_results <- expand.grid(Replication = 1:n_replications, Link = links, Estimator = estimators, stringsAsFactors = FALSE)
sim_results$MSE <- 0
sim_results$Bias <- 0
sim_results$Active <- 0
sim_results$TPR <- 0
sim_results$FPR <- 0
sim_results$Time_ms <- 0

cat("==============================================================================\n")
cat("EXECUTING MONTE CARLO SIMULATION SUITE...\n")
cat("==============================================================================\n")

total_tasks <- nrow(sim_results)
for (i in 1:total_tasks) {
  r <- sim_results$Replication[i]
  l <- sim_results$Link[i]
  m <- sim_results$Estimator[i]
  
  cat(sprintf("\r Progress: [%d/%d] | Rep: %d | Link: %-6s | Method: %-6s", i, total_tasks, r, l, m))
  flush.console()
  
  X_raw <- matrix(rnorm(n_samples * p_features), nrow = n_samples, ncol = p_features)
  X_sim <- scale(X_raw) / sqrt(p_features)
  
  mu_sim <- compute_mu(as.vector(X_sim %*% beta_true), link = l)
  y_sim <- qUR(runif(n_samples), mu_sim, tau = tau_level)
  
  t_start <- Sys.time()
  b_est <- fit_urrm(X_sim, y_sim, link = l, method = m, tau = tau_level, lambda = 0.02)
  t_end <- Sys.time()
  
  sim_results$MSE[i] <- mean((b_est - beta_true)^2)
  sim_results$Bias[i] <- mean(b_est - beta_true)
  sim_results$Active[i] <- sum(b_est != 0)
  sim_results$TPR[i] <- sum(b_est[1:4] != 0) / 4
  sim_results$FPR[i] <- sum(b_est[5:p_features] != 0) / (p_features - 4)
  sim_results$Time_ms[i] <- as.numeric(difftime(t_end, t_start, units = "secs")) * 1000
}

cat("\n\n✓ Monte Carlo Simulations Completed!\n\n")

# ------------------------------------------------------------------------------
# 4. EXPORT ALL 5 TABLES
# ------------------------------------------------------------------------------

# Table 1: Parameter Estimation Accuracy & Bias
table_1 <- sim_results %>%
  group_by(Link, Estimator) %>%
  summarise(Mean_MSE = round(mean(MSE), 4), SD_MSE = round(sd(MSE), 4), Mean_Bias = round(mean(Bias), 4), .groups = "drop")
write.csv(table_1, "Table_1_Parameter_Estimation_Bias.csv", row.names = FALSE)

# Table 2: Feature Selection Metrics
table_2 <- sim_results %>%
  group_by(Link, Estimator) %>%
  summarise(Avg_Active = round(mean(Active), 2), TPR = round(mean(TPR), 3), FPR = round(mean(FPR), 3), .groups = "drop")
write.csv(table_2, "Table_2_Feature_Selection_Performance.csv", row.names = FALSE)

# Table 3: Numerical Efficiency
table_3 <- sim_results %>%
  group_by(Link, Estimator) %>%
  summarise(Avg_Time_ms = round(mean(Time_ms), 2), SD_Time_ms = round(sd(Time_ms), 2), .groups = "drop")
write.csv(table_3, "Table_3_Numerical_Efficiency.csv", row.names = FALSE)

# Table 4: Simulation Design & Response Summary
set.seed(2026)
X_bench_raw <- matrix(rnorm(100 * 20), 100, 20)
X_bench <- scale(X_bench_raw) / sqrt(20)
b_bench_true <- c(0.8, -0.6, 0.5, -0.4, rep(0, 16))
y_bench <- qUR(runif(100), compute_mu(as.vector(X_bench %*% b_bench_true), "loglog"), tau = 0.5)

table_4 <- data.frame(
  Metric = c("Sample Size (n)", "Features (p)", "Signal Predictors", "Response Mean", "Response SD", "Response Median", "Response Min", "Response Max"),
  Value = c(100, 20, 4, round(mean(y_bench), 4), round(sd(y_bench), 4), round(median(y_bench), 4), round(min(y_bench), 4), round(max(y_bench), 4))
)
write.csv(table_4, "Table_4_Simulation_Design_Summary.csv", row.names = FALSE)

# Table 5: 16-Model Selection Grid
table_5_grid <- expand.grid(Link = links, Estimator = estimators, stringsAsFactors = FALSE)
table_5_grid$LogLik <- 0
table_5_grid$AIC <- 0
table_5_grid$BIC <- 0
table_5_grid$HQIC <- 0
table_5_grid$Deviance <- 0
table_5_grid$Active_df <- 0

for (i in 1:nrow(table_5_grid)) {
  l <- table_5_grid$Link[i]
  m <- table_5_grid$Estimator[i]
  
  b_hat <- fit_urrm(X_bench, y_bench, link = l, method = m, tau = 0.5, lambda = 0.02)
  mu_hat <- compute_mu(as.vector(X_bench %*% b_hat), link = l)
  
  ll <- logLik_UR(y_bench, mu_hat, tau = 0.5)
  df <- max(sum(b_hat != 0), 1)
  
  table_5_grid$LogLik[i]    <- round(ll, 2)
  table_5_grid$AIC[i]       <- round(-2 * ll + 2 * df, 2)
  table_5_grid$BIC[i]       <- round(-2 * ll + df * log(100), 2)
  table_5_grid$HQIC[i]      <- round(-2 * ll + 2 * df * log(log(100)), 2)
  table_5_grid$Deviance[i]  <- round(-2 * ll, 2)
  table_5_grid$Active_df[i] <- df
}

table_5_grid <- table_5_grid %>% arrange(AIC)
write.csv(table_5_grid, "Table_5_Simulation_Model_Selection_Criteria.csv", row.names = FALSE)

cat("✓ All 5 CSV Tables successfully generated.\n\n")

# ------------------------------------------------------------------------------
# 5. GENERATE ALL 6 PUBLICATION FIGURES
# ------------------------------------------------------------------------------

cat("Generating All 6 High-Resolution PNG Figures...\n")

# FIGURE 1: MSE Boxplots
f1 <- ggplot(sim_results, aes(x = Estimator, y = MSE, fill = Estimator)) +
  geom_boxplot(alpha = 0.8, outlier.size = 1) +
  facet_wrap(~ Link, scales = "free_y") +
  theme_bw() + labs(title = "Figure 1: Simulation MSE Performance Across Link Functions", y = "Mean Squared Error")
ggsave("Figure_1_MSE_Performance.png", f1, width = 9, height = 6, dpi = 300)

# FIGURE 2: SCAD Regularization Path
lambdas <- exp(seq(log(0.001), log(0.15), length.out = 25))
path_list <- list()
for (idx in length(lambdas):1) {
  l_val <- lambdas[idx]
  b_curr <- fit_hd_penalized(X_bench, y_bench, tau = 0.5, lambda = l_val, penalty = "scad", link = "loglog")
  path_list[[idx]] <- data.frame(Lambda = l_val, Beta = b_curr[1:5], Feature = paste0("X", 1:5))
}
path_df <- do.call(rbind, path_list)

f2 <- ggplot(path_df, aes(x = Lambda, y = Beta, color = Feature)) +
  geom_line(size = 1.1) + scale_x_log10() + theme_bw() +
  labs(title = "Figure 2: SCAD Regularization Paths for Top Predictors", x = "Penalty Parameter Log(Lambda)", y = "Estimated Beta")
ggsave("Figure_2_Regularization_Paths.png", f2, width = 8, height = 5, dpi = 300)

# FIGURE 3: PERFECTLY CONDITIONED RESIDUAL DIAGNOSTICS
nll_urrm <- function(beta, X, y, link = "loglog", tau = 0.5) {
  eta <- as.vector(X %*% beta)
  mu <- compute_mu(eta, link)
  -sum(log(dUR(y, mu, tau = tau)))
}

fit_mle <- optim(par = rep(0, 20), fn = nll_urrm, X = X_bench, y = y_bench, 
                 link = "loglog", tau = 0.5, method = "BFGS")
mu_diag <- compute_mu(as.vector(X_bench %*% fit_mle$par), link = "loglog")

u_cdf <- pUR(y_bench, mu_diag, tau = 0.5)
u_cdf <- pmin(pmax(u_cdf, 0.001), 0.999)
rq_residuals <- qnorm(u_cdf)

df_diag <- data.frame(Fitted = mu_diag, RQ_Residuals = rq_residuals)

p_res1 <- ggplot(df_diag, aes(x = Fitted, y = RQ_Residuals)) +
  geom_point(color = "#2b5c8f", size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
  geom_smooth(method = "loess", color = "darkorange", se = FALSE, size = 0.8) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() + 
  labs(title = "Randomized Quantile Residuals vs Fitted", x = "Fitted Values (mu)", y = "Quantile Residuals")

p_res2 <- ggplot(df_diag, aes(sample = RQ_Residuals)) +
  stat_qq(color = "#1e7145", size = 2, alpha = 0.8) + 
  stat_qq_line(color = "red", size = 0.8) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() + 
  labs(title = "Normal Q-Q Plot", x = "Theoretical Quantiles", y = "Sample Quantiles")

f3 <- grid.arrange(p_res1, p_res2, ncol = 2)
ggsave("Figure_3_Residual_Diagnostics.png", f3, width = 10, height = 4.5, dpi = 300)

# FIGURE 4: Leverage & Influence
H <- X_bench %*% solve(t(X_bench) %*% X_bench + 0.1 * diag(20)) %*% t(X_bench)
lev <- diag(H)
cooks <- (rq_residuals^2 / 20) * (lev / (1 - lev)^2)

f4_df <- data.frame(Index = 1:100, Leverage = lev, Cooks_D = cooks)
p_inf1 <- ggplot(f4_df, aes(x = Index, y = Leverage)) + geom_segment(aes(xend = Index, yend = 0)) + geom_point(color = "blue") + theme_bw() + labs(title = "Leverage (h_ii)")
p_inf2 <- ggplot(f4_df, aes(x = Index, y = Cooks_D)) + geom_segment(aes(xend = Index, yend = 0)) + geom_point(color = "darkred") + theme_bw() + labs(title = "Cook's Distance")
f4 <- grid.arrange(p_inf1, p_inf2, ncol = 2)
ggsave("Figure_4_Leverage_CooksD.png", f4, width = 10, height = 4.5, dpi = 300)

# FIGURE 5: CORRECTED PAIRED SIMULATION FIT
f5_df <- data.frame(Fitted = mu_diag, Observed = y_bench) %>%
  arrange(Fitted) %>%
  mutate(Obs = 1:n())

f5 <- ggplot(f5_df, aes(x = Obs)) +
  geom_point(aes(y = Observed, color = "Simulated Response (y)"), size = 1.8, alpha = 0.8) +
  geom_line(aes(y = Fitted, color = "URRM Fitted Mu"), size = 1.2) +
  scale_color_manual(values = c("Simulated Response (y)" = "#f8766d", "URRM Fitted Mu" = "#00bfc4")) +
  theme_bw() + 
  labs(title = "Figure 5: Model Fit vs Generated Bounded Response", 
       x = "Observation Index (Ordered by Fitted Values)", y = "Response (y)", color = "Legend")
ggsave("Figure_5_Empirical_Fit.png", f5, width = 8, height = 5, dpi = 300)

# FIGURE 6: Computation Time Efficiency
f6 <- ggplot(sim_results, aes(x = Estimator, y = Time_ms, fill = Link)) +
  geom_bar(stat = "summary", fun = "mean", position = "dodge") +
  theme_bw() + labs(title = "Figure 6: Execution Time Efficiency (Milliseconds)", y = "Time (ms)")
ggsave("Figure_6_Computation_Efficiency.png", f6, width = 8, height = 5, dpi = 300)

cat("\n==============================================================================\n")
cat("SUCCESS! All 5 CSV Tables & 6 Figures Saved in Directory:\n", getwd(), "\n")
cat("==============================================================================\n")
