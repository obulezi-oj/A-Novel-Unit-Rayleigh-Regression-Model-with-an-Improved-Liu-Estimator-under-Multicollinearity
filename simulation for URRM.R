# ==============================================================================
# Comprehensive Simulation & Diagnostics Script for Unit-Rayleigh Regression (URRM)
# STABLE VERSION: Monte Carlo Simulation, Manuscript Tables & Model Diagnostics
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. ENVIRONMENT & DIRECTORY SETUP
# ------------------------------------------------------------------------------

# Set working directory safely
work_dir <- "C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\URRM"
if (dir.exists(work_dir)) {
  setwd(work_dir)
} else {
  warning("Specified working directory does not exist. Using current directory: ", getwd())
}

# Define output directories using absolute pathing
fig_dir <- file.path(getwd(), "figures")
res_dir <- file.path(getwd(), "results")

# Create output subdirectories recursively to prevent file connection errors
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(res_dir)) dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(MASS)
  library(gridExtra)
})

# ==============================================================================
# 1. CORE DISTRIBUTION & REGRESSION FUNCTIONS
# ==============================================================================

#' Quantile Function Q(u; mu, tau)
qUR <- function(u, mu, tau = 0.5) {
  exp(log(mu) * sqrt(log(u) / log(tau)))
}

#' Cumulative Distribution Function F(x; mu, tau)
pUR <- function(x, mu, tau = 0.5) {
  exp(log(tau) * (log(x) / log(mu))^2)
}

#' Probability Density Function f(x; mu, tau)
dUR <- function(x, mu, tau = 0.5) {
  term1 <- 2 * log(tau) * log(x) / (x * (log(mu))^2)
  term2 <- exp(log(tau) * (log(x) / log(mu))^2)
  return(term1 * term2)
}

#' Hazard Rate Function h(x; mu, tau) = f(x) / (1 - F(x))
hUR <- function(x, mu, tau = 0.5) {
  fx <- dUR(x, mu, tau)
  Fx <- pUR(x, mu, tau)
  return(fx / pmax(1 - Fx, 1e-10))
}

#' Negative Log-Likelihood Function for URRM
nll_urrm <- function(beta, X, y, tau = 0.5) {
  eta <- as.vector(X %*% beta)
  mu <- exp(-exp(eta))
  mu <- pmin(pmax(mu, 1e-10), 1 - 1e-10)
  
  ll <- log(-2 * log(tau)) + log(-log(y)) - log(y) - 2 * eta + log(tau) * (log(y))^2 * exp(-2 * eta)
  return(-sum(ll))
}

#' Maximum Likelihood Estimation (MLE)
fit_mle <- function(X, y, tau = 0.5) {
  p <- ncol(X)
  init_beta <- rep(0, p)
  opt <- optim(init_beta, nll_urrm, X = X, y = y, tau = tau, method = "BFGS")
  return(opt$par)
}

#' Compute 6 Regularized Shrinkage Parameters d (eps = 0.01)
get_d_values <- function(beta_mle, X_design, eps = 0.01) {
  Info <- 4 * t(X_design) %*% X_design
  eig <- eigen(Info)
  Q <- eig$vectors
  lambda <- eig$values
  
  alpha <- as.vector(t(Q) %*% beta_mle)
  
  d1_vec <- (alpha^2 - 1) / ((1 / lambda) + alpha^2)
  d1_raw <- min(d1_vec)
  
  d2_raw <- sum((alpha^2 - 1) / (lambda + 1)^2) / sum(((1 / lambda) + alpha^2) / (lambda + 1)^2)
  d3_raw <- median(d1_vec)
  d4_raw <- (max(alpha^2) - 1) / ((1 / min(lambda)) + max(alpha^2))
  d5_raw <- 1 - (sum(1 / (lambda + 1)) / sum(alpha^2 / (lambda + 1)))
  d6_raw <- sum(lambda * (alpha^2 - 1)) / sum(1 + lambda * alpha^2)
  
  d1 <- min(1 - eps, max(eps, d1_raw))
  d2 <- min(1 - eps, max(eps, d2_raw))
  d3 <- min(1 - eps, max(eps, d3_raw))
  d4 <- min(1 - eps, max(eps, d4_raw))
  d5 <- min(1 - eps, max(eps, d5_raw))
  d6 <- min(1 - eps, max(eps, d6_raw))
  
  return(c(d1 = d1, d2 = d2, d3 = d3, d4 = d4, d5 = d5, d6 = d6))
}

#' Compute Improved Liu Estimator (URLRE)
fit_urlre <- function(beta_mle, X_design, d) {
  p <- length(beta_mle)
  Info <- 4 * t(X_design) %*% X_design
  Fd <- solve(Info + diag(p)) %*% (Info + d * diag(p))
  beta_urlre <- as.vector(Fd %*% beta_mle)
  return(beta_urlre)
}

# ==============================================================================
# 2. MONTE CARLO SIMULATION STUDY (MSE, MAE, RE, d ESTIMATES)
# ==============================================================================
set.seed(2026)

sample_sizes <- c(20, 25, 50, 100, 250, 500)
rhos <- c(0.90, 0.99, 0.999)
tau <- 0.5
true_beta <- c(0.5, -0.8, 1.2, -0.5)
k <- length(true_beta) - 1
n_sim <- 500

mse_results_list <- list()
mae_results_list <- list()
re_results_list  <- list()
d_vals_results_list <- list()

cat("Executing Monte Carlo Simulations...\n")

for (n in sample_sizes) {
  for (rho in rhos) {
    cat(sprintf(" -> Processing n = %d, rho = %.3f\n", n, rho))
    
    mse_matrix <- matrix(0, nrow = n_sim, ncol = 7)
    colnames(mse_matrix) <- c("MLE", "d1", "d2", "d3", "d4", "d5", "d6")
    
    mae_matrix <- matrix(0, nrow = n_sim, ncol = 7)
    colnames(mae_matrix) <- c("MLE", "d1", "d2", "d3", "d4", "d5", "d6")
    
    d_matrix <- matrix(0, nrow = n_sim, ncol = 6)
    colnames(d_matrix) <- paste0("d", 1:6)
    
    for (m in 1:n_sim) {
      W <- matrix(rnorm(n * (k + 1)), nrow = n, ncol = k + 1)
      Z_preds <- matrix(0, nrow = n, ncol = k)
      for (j in 1:k) {
        Z_preds[, j] <- sqrt(1 - rho^2) * W[, j] + rho * W[, k + 1]
      }
      X_design <- cbind(1, Z_preds)
      
      eta <- as.vector(X_design %*% true_beta)
      mu_true <- exp(-exp(eta))
      
      u <- runif(n)
      y <- qUR(u, mu_true, tau = tau)
      
      b_mle <- fit_mle(X_design, y, tau = tau)
      d_hats <- get_d_values(b_mle, X_design, eps = 0.01)
      d_matrix[m, ] <- d_hats
      
      # MLE Metrics
      mse_matrix[m, "MLE"] <- sum((b_mle - true_beta)^2)
      mae_matrix[m, "MLE"] <- sum(abs(b_mle - true_beta))
      
      # URLRE Metrics
      for (idx in 1:6) {
        b_urlre <- fit_urlre(b_mle, X_design, d_hats[idx])
        col_name <- paste0("d", idx)
        mse_matrix[m, col_name] <- sum((b_urlre - true_beta)^2)
        mae_matrix[m, col_name] <- sum(abs(b_urlre - true_beta))
      }
    }
    
    avg_mse <- colMeans(mse_matrix)
    avg_mae <- colMeans(mae_matrix)
    avg_d   <- colMeans(d_matrix)
    
    # Relative Efficiency (RE) = MSE_MLE / MSE_URLRE
    re_vals <- avg_mse["MLE"] / avg_mse
    
    mse_results_list[[length(mse_results_list) + 1]] <- data.frame(
      Sample_Size = n, Rho = rho, t(avg_mse)
    )
    
    mae_results_list[[length(mae_results_list) + 1]] <- data.frame(
      Sample_Size = n, Rho = rho, t(avg_mae)
    )
    
    re_results_list[[length(re_results_list) + 1]] <- data.frame(
      Sample_Size = n, Rho = rho, t(re_vals)
    )
    
    d_vals_results_list[[length(d_vals_results_list) + 1]] <- data.frame(
      Sample_Size = n, Rho = rho, t(avg_d)
    )
  }
}

# Compile DataFrames
sim_mse_df <- do.call(rbind, mse_results_list)
sim_mae_df <- do.call(rbind, mae_results_list)
sim_re_df  <- do.call(rbind, re_results_list)
sim_d_df   <- do.call(rbind, d_vals_results_list)

# ==============================================================================
# 3. EXPORT MANUSCRIPT TABLES (CSVs)
# ==============================================================================
cat("\nExporting 4 Comprehensive Tables to /results directory...\n")

# Table 1: Mean Squared Errors (MSE)
write.csv(sim_mse_df, file.path(res_dir, "Table1_Mean_Squared_Error_MSE.csv"), row.names = FALSE)

# Table 2: Relative Efficiency (RE)
write.csv(sim_re_df, file.path(res_dir, "Table2_Relative_Efficiency_RE.csv"), row.names = FALSE)

# Table 3: Mean Absolute Error (MAE)
write.csv(sim_mae_df, file.path(res_dir, "Table3_Mean_Absolute_Error_MAE.csv"), row.names = FALSE)

# Table 4: Optimal d Parameter Estimates
write.csv(sim_d_df, file.path(res_dir, "Table4_Optimal_d_Parameter_Estimates.csv"), row.names = FALSE)

# ==============================================================================
# 4. MSE CONVERGENCE & PARAMETER D PLOTS (FIGURES 3 & 4)
# ==============================================================================

# Figure 3: MSE Convergence Line Plots
df_line <- sim_mse_df %>%
  pivot_longer(cols = c("MLE", paste0("d", 1:6)), names_to = "Estimator", values_to = "MSE") %>%
  mutate(Rho_Fct = factor(paste0("rho = ", Rho)))

fig3_line <- ggplot(df_line, aes(x = Sample_Size, y = MSE, color = Estimator, group = Estimator)) +
  geom_line(size = 0.9) +
  geom_point(size = 2) +
  scale_y_log10() +
  facet_wrap(~ Rho_Fct, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "MSE Convergence as Sample Size (n) Increases",
       subtitle = "Log-scale MSE trajectory across collinearity levels",
       x = "Sample Size (n)", y = "Mean Squared Error (Log Scale)") +
  theme_minimal(base_size = 11)

ggsave(file.path(fig_dir, "Fig3_MSE_Convergence_Lines.png"), fig3_line, width = 10, height = 6, dpi = 300)

# Figure 4: Average Estimated Shrinkage Parameters (d1 - d6)
df_d_long <- sim_d_df %>%
  pivot_longer(cols = starts_with("d"), names_to = "d_Parameter", values_to = "Value") %>%
  mutate(Rho_Fct = factor(paste0("rho = ", Rho)),
         Sample_Size_Fct = factor(paste0("n = ", Sample_Size), levels = paste0("n = ", sample_sizes)))

fig4_d_vals <- ggplot(df_d_long, aes(x = d_Parameter, y = Value, fill = d_Parameter)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Rho_Fct ~ Sample_Size_Fct) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Average Bounded Shrinkage Parameter Estimates (d1 - d6)",
       subtitle = "Behaviors under boundary regularization (eps = 0.01)",
       x = "Shrinkage Estimator", y = "Average Estimated Value") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Fig4_Optimal_d_Estimates.png"), fig4_d_vals, width = 12, height = 7, dpi = 300)

# ==============================================================================
# 5. FIGURE 5: COMPLETE MULTI-PANEL BAR COMPARISON
# ==============================================================================

df_bar_all <- sim_mse_df %>%
  pivot_longer(cols = c("MLE", paste0("d", 1:6)), names_to = "Estimator", values_to = "MSE") %>%
  mutate(Sample_Size_Fct = factor(paste0("n = ", Sample_Size), levels = paste0("n = ", sample_sizes)),
         Rho_Fct = factor(paste0("rho = ", Rho)))

fig5_bar_all <- ggplot(df_bar_all, aes(x = Estimator, y = MSE, fill = Estimator)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  facet_grid(Rho_Fct ~ Sample_Size_Fct, scales = "free_y") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "MSE Comparison Across All Sample Sizes and Multicollinearity Degrees",
       subtitle = "Evaluating MLE vs. URLRE (d1-d6) strategies",
       x = "Estimator Variant", y = "Mean Squared Error (MSE)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "#f0f0f0", color = NA),
        panel.grid.minor = element_blank())

ggsave(file.path(fig_dir, "Fig5_Estimators_MSE_All_Sample_Sizes.png"), fig5_bar_all, width = 14, height = 8, dpi = 300)

# ==============================================================================
# 6. FIGURE 6: DIAGNOSTIC RESIDUAL PLOTS
# ==============================================================================
cat("Generating Diagnostic Residual Plots (Figure 6)...\n")

n_diag <- 250
rho_diag <- 0.99
W <- matrix(rnorm(n_diag * (k + 1)), nrow = n_diag, ncol = k + 1)
Z_preds <- matrix(0, nrow = n_diag, ncol = k)
for (j in 1:k) Z_preds[, j] <- sqrt(1 - rho_diag^2) * W[, j] + rho_diag * W[, k + 1]
X_diag <- cbind(1, Z_preds)
eta_diag <- as.vector(X_diag %*% true_beta)
mu_diag <- exp(-exp(eta_diag))
y_diag <- qUR(runif(n_diag), mu_diag, tau = tau)

b_mle_diag <- fit_mle(X_diag, y_diag, tau = tau)
mu_hat <- exp(-exp(as.vector(X_diag %*% b_mle_diag)))

rq_residuals <- qnorm(pUR(y_diag, mu_hat, tau = tau))
cs_residuals <- -log(1 - pUR(y_diag, mu_hat, tau = tau))

df_diag <- data.frame(Fitted = mu_hat, RQ_Res = rq_residuals, CS_Res = cs_residuals)

# Panel A: Normal Q-Q Plot
p_diag1 <- ggplot(df_diag, aes(sample = RQ_Res)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Panel A: Q-Q Plot (Quantile Residuals)", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

# Panel B: Residuals vs Fitted Values
p_diag2 <- ggplot(df_diag, aes(x = Fitted, y = RQ_Res)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Panel B: Residuals vs Fitted", x = expression(hat(mu)), y = "Randomized Quantile Residuals") +
  theme_minimal()

# Panel C: Cox-Snell Residual Plot
df_cs <- data.frame(Observed = sort(cs_residuals), Theoretical = qexp(ppoints(n_diag)))
p_diag3 <- ggplot(df_cs, aes(x = Theoretical, y = Observed)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "Panel C: Cox-Snell Residual Plot", x = "Theoretical Exponential Quantiles", y = "Observed Residuals") +
  theme_minimal()

fig6_diag <- grid.arrange(p_diag1, p_diag2, p_diag3, ncol = 3)
ggsave(file.path(fig_dir, "Fig6_Residual_Diagnostics.png"), fig6_diag, width = 14, height = 4.5, dpi = 300)

cat("\nScript execution complete. Output tables (1-4) and Figures (3-6) saved successfully.\n")