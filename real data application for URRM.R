# ==============================================================================
# Empirical Application Script: URRM vs Unit-Lindley QR (ULRM)
# Path: C:/Users/Dr. O. J. Obulezi/Documents/R projects/URRM
# Output: Single Master Table for Quantiles + Summary Tables & Plots
# ==============================================================================

# 1. Set Directory Hierarchy
base_dir <- "C:/Users/Dr. O. J. Obulezi/Documents/R projects/URRM"
if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
setwd(base_dir)

dir.create("output_results", showWarnings = FALSE)
dir.create("output_results/plots", showWarnings = FALSE)
dir.create("output_results/tables", showWarnings = FALSE)

# 2. Package Setup
required_pkgs <- c("betareg", "car", "ggplot2", "gridExtra")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

library(betareg)
library(car)
library(ggplot2)
library(gridExtra)

# 3. Data Preparation & Diagnostics
data("GasolineYield", package = "betareg")
df <- GasolineYield

y <- df$yield
n <- length(y)

Z_raw <- as.matrix(df[, c("gravity", "pressure", "temp10", "temp")])
Z <- cbind(Intercept = 1, Z_raw)
k <- ncol(Z_raw)
p <- ncol(Z) # p = 5

ols_fit <- lm(yield ~ gravity + pressure + temp10 + temp, data = df)
vif_vals <- car::vif(ols_fit)

xtx <- t(Z) %*% Z
eigen_vals <- eigen(xtx)$values
condition_index <- sqrt(max(eigen_vals) / min(eigen_vals))

mc_table <- data.frame(
  Variable = c(names(vif_vals), "Condition Index (CI)"),
  Value = c(round(vif_vals, 4), round(condition_index, 4))
)
write.csv(mc_table, "output_results/tables/Table1_Multicollinearity_Diagnostics.csv", row.names = FALSE)

# 4. Model Functions
ur_loglik <- function(beta, y, Z, tau = 0.5) {
  eta <- Z %*% beta
  mu <- 1 / (1 + exp(-eta))
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  theta <- -log(tau) / ((log(mu))^2)
  log_density <- log(2) + log(theta) - log(y) - log(-log(y)) - (theta / (log(y))^2)
  return(sum(log_density))
}

fit_urrm_mle <- function(y, Z, tau = 0.5) {
  init_beta <- coef(lm(qlogis(y) ~ Z[, -1]))
  opt <- optim(par = init_beta, 
               fn = function(b) -ur_loglik(b, y = y, Z = Z, tau = tau), 
               method = "BFGS", hessian = TRUE)
  beta_hat <- opt$par
  vcov_mle <- solve(opt$hessian)
  se_mle <- sqrt(diag(vcov_mle))
  list(beta = beta_hat, se = se_mle, vcov = vcov_mle, loglik = -opt$value)
}

ul_loglik <- function(beta, y, Z, tau = 0.5) {
  eta <- Z %*% beta
  mu <- 1 / (1 + exp(-eta))
  mu <- pmin(pmax(mu, 1e-6), 1 - 1e-6)
  theta <- (1 - mu) / mu
  log_density <- 2 * log(theta) - log(1 + theta) - 3 * log(1 - y) - (theta * y / (1 - y))
  return(sum(log_density))
}

fit_ulrm_mle <- function(y, Z, tau = 0.5) {
  init_beta <- coef(lm(qlogis(y) ~ Z[, -1]))
  opt <- optim(par = init_beta, 
               fn = function(b) -ul_loglik(b, y = y, Z = Z, tau = tau), 
               method = "BFGS", hessian = TRUE)
  list(beta = opt$par, loglik = -opt$value)
}

fit_urlre <- function(beta_mle, Z, d) {
  FIM <- t(Z) %*% Z 
  I_p <- diag(p)
  W_d <- solve(FIM + I_p) %*% (FIM + d * I_p)
  return(as.vector(W_d %*% beta_mle))
}

get_optimal_d <- function(beta_mle, vcov_mle, Z) {
  sigma2_hat <- mean((y - 1/(1 + exp(-Z %*% beta_mle)))^2)
  alpha <- beta_mle
  lambda <- eigen(t(Z) %*% Z)$values
  
  d1 <- max(0, min(1, median((alpha^2 - sigma2_hat) / (alpha^2 + sigma2_hat / lambda))))
  d2 <- max(0, min(1, mean((alpha^2) / (alpha^2 + sigma2_hat / lambda))))
  d3 <- max(0, min(1, min((alpha^2) / (alpha^2 + 2 * sigma2_hat / lambda))))
  d4 <- max(0, min(1, 1 - sigma2_hat * mean(1 / (lambda * alpha^2 + sigma2_hat))))
  d5 <- max(0, min(1, median(1 - (sigma2_hat / (lambda * alpha^2)))))
  d6 <- max(0, min(1, min(lambda * alpha^2 / (sigma2_hat + lambda * alpha^2))))
  
  return(c(d1 = d1, d2 = d2, d3 = d3, d4 = d4, d5 = d5, d6 = d6))
}

calc_ic <- function(loglik, p, n) {
  aic <- -2 * loglik + 2 * p
  aicc <- aic + (2 * p * (p + 1)) / (n - p - 1)
  bic <- -2 * loglik + p * log(n)
  hqic <- -2 * loglik + 2 * p * log(log(n))
  return(c(LogLik = loglik, AIC = aic, AICC = aicc, BIC = bic, HQIC = hqic))
}

# 5. Fit Comparison (URRM vs ULRM at tau = 0.50)
mle_urrm_med <- fit_urrm_mle(y, Z, tau = 0.5)
mle_ulrm_med <- fit_ulrm_mle(y, Z, tau = 0.5)

comp_table <- data.frame(
  Model = c("Unit-Rayleigh QR (URRM)", "Unit-Lindley QR (ULRM)"),
  rbind(
    round(calc_ic(mle_urrm_med$loglik, p = 5, n = n), 3),
    round(calc_ic(mle_ulrm_med$loglik, p = 5, n = n), 3)
  )
)
write.csv(comp_table, "output_results/tables/Table2_Model_Fit_Comparison.csv", row.names = FALSE)

# 6. SINGLE CONSOLIDATED TABLE FOR ALL QUANTILES (tau = 0.10 to 0.90)
quantiles_all <- seq(0.10, 0.90, by = 0.10)
master_quantile_list <- list()
summary_process_list <- list()

for (tau in quantiles_all) {
  mle_fit <- fit_urrm_mle(y, Z, tau = tau)
  beta_mle <- mle_fit$beta
  se_mle <- mle_fit$se
  
  d_vec <- get_optimal_d(beta_mle, mle_fit$vcov, Z)
  
  liu_matrix <- matrix(NA, nrow = p, ncol = 6)
  colnames(liu_matrix) <- paste0("URLRE_d", 1:6)
  for (j in 1:6) {
    liu_matrix[, j] <- fit_urlre(beta_mle, Z, d = d_vec[j])
  }
  
  # Consolidated block for current quantile
  tau_block <- data.frame(
    Tau = sprintf("%.2f", tau),
    Parameter = c("Intercept", "Gravity", "Pressure", "Temp10", "Temp"),
    MLE = round(beta_mle, 4),
    MLE_SE = round(se_mle, 4),
    round(liu_matrix, 4)
  )
  
  master_quantile_list[[as.character(tau)]] <- tau_block
  
  # Process summary row
  summary_process_list[[as.character(tau)]] <- data.frame(
    Quantile = tau,
    Intercept_MLE = beta_mle[1],
    Gravity_MLE = beta_mle[2],
    Pressure_MLE = beta_mle[3],
    Temp10_MLE = beta_mle[4],
    Temp_MLE = beta_mle[5],
    LogLik = mle_fit$loglik
  )
}

# Bind into ONE SINGLE MASTER TABLE
consolidated_quantile_table <- do.call(rbind, master_quantile_list)

# Export the single big table
write.csv(consolidated_quantile_table, "output_results/tables/Table3_Consolidated_Quantile_Estimates.csv", row.names = FALSE)

summary_quant_df <- do.call(rbind, summary_process_list)

# 7. Plots
# Figure 1: Model Diagnostics
mu_hat <- 1 / (1 + exp(-Z %*% mle_urrm_med$beta))
theta_hat <- -log(0.5) / ((log(mu_hat))^2)
F_y <- exp(-theta_hat / ((log(y))^2))
quant_residuals <- qnorm(pmin(pmax(F_y, 1e-6), 1 - 1e-6))

diag_df <- data.frame(
  Index = 1:n,
  Fitted = mu_hat,
  QuantileResiduals = quant_residuals
)

p1 <- ggplot(diag_df, aes(sample = QuantileResiduals)) +
  stat_qq(color = "#1F78B4", size = 2) +
  stat_qq_line(color = "#E31A1C", linetype = "dashed", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(title = "(a) Q-Q Plot of Quantile Residuals",
       x = "Theoretical Normal Quantiles", y = "Sample Quantiles")

p2 <- ggplot(diag_df, aes(x = Fitted, y = QuantileResiduals)) +
  geom_point(color = "#33A02C", size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, color = "#FF7F00", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(title = "(b) Residuals vs Fitted Values",
       x = "Fitted Values (Predicted Quantile)", y = "Randomized Quantile Residuals")

d_seq <- seq(0, 1, length.out = 100)
coef_path <- t(sapply(d_seq, function(d) fit_urlre(mle_urrm_med$beta, Z, d)[-1]))
path_df <- data.frame(d = d_seq, coef_path)
colnames(path_df) <- c("d", "Gravity", "Pressure", "Temp10", "Temp")

path_long <- reshape(path_df, 
                     varying = c("Gravity", "Pressure", "Temp10", "Temp"),
                     v.names = "Coefficient",
                     timevar = "Predictor",
                     times = c("Gravity", "Pressure", "Temp10", "Temp"),
                     direction = "long")

p3 <- ggplot(path_long, aes(x = d, y = Coefficient, color = Predictor)) +
  geom_line(linewidth = 1.2) +
  theme_minimal(base_size = 12) +
  labs(title = "(c) URLRE Shrinkage Trace Plot (tau = 0.50)",
       x = "Liu Shrinkage Parameter (d)", y = "Coefficient Estimate") +
  theme(legend.position = "bottom")

png("output_results/plots/Figure1_Model_Diagnostics.png", width = 10, height = 8, units = "in", res = 300)
grid.arrange(p1, p2, p3, layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

# Figure 2: Quantile Trajectories Across Predictors
q_process_long <- reshape(summary_quant_df[, c("Quantile", "Gravity_MLE", "Pressure_MLE", "Temp10_MLE", "Temp_MLE")],
                          varying = c("Gravity_MLE", "Pressure_MLE", "Temp10_MLE", "Temp_MLE"),
                          v.names = "Estimate",
                          timevar = "Predictor",
                          times = c("Gravity", "Pressure", "Temp10", "Temp"),
                          direction = "long")

p_quantiles <- ggplot(q_process_long, aes(x = Quantile, y = Estimate)) +
  geom_line(color = "#1F78B4", linewidth = 1.2) +
  geom_point(color = "#E31A1C", size = 2.5) +
  facet_wrap(~ Predictor, scales = "free_y") +
  theme_bw(base_size = 12) +
  labs(title = "Quantile Trajectories of Parameter Estimates across Quantiles (0.10 to 0.90)",
       x = expression(Quantile ~ Level ~ (tau)), y = "MLE Estimate")

png("output_results/plots/Figure2_Quantile_Process_Plots.png", width = 9, height = 7, units = "in", res = 300)
print(p_quantiles)
dev.off()

cat("========================================================\n")
cat("SUCCESS! Cleaned output directory contains ONLY 3 consolidated tables:\n")
cat(" 1. Table1_Multicollinearity_Diagnostics.csv\n")
cat(" 2. Table2_Model_Fit_Comparison.csv\n")
cat(" 3. Table3_Consolidated_Quantile_Estimates.csv (MASTER TABLE)\n")
cat("========================================================\n")