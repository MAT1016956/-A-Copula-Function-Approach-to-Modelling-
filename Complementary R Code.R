#PROBLEM SECTION 2(ADJUSTED CODE)
library(copula)
library(VGAM)
library(ggplot2)
theta <- 2
claytonCop <- claytonCopula(dim = 2, param = theta)
wind_speed_dist <- function(n) qweibull(runif(n), shape = 2, scale = 10)
wind_direction_dist <- function(n) runif(n, min = 0, max = 360)
set.seed(123) 
n_samples <- 1000
copula_sample <- rCopula(n_samples, claytonCop)
u <- qweibull(copula_sample[,1], shape = 2, scale = 10) 
v <- 360 * copula_sample[,2]  
u_cop <- pweibull(u, shape = 2, scale = 10)
v_cop <- (v - min(v)) / (max(v) - min(v))
upper_bound <- pmin(u_cop, v_cop)
copula_df <- data.frame(W = u, D = v)
ggplot(copula_df, aes(x = W, y = D)) +
  geom_point(alpha = 0.4) +
  labs(title = "Joint Distribution of Wind Speed and Direction",
       x = "Wind Speed (m/s)", y = "Wind Direction (degrees)") +
  theme_minimal()
bound_df <- data.frame(W = u, Bounds = upper_bound)
ggplot(bound_df, aes(x = W, y = Bounds)) +
  geom_line(color = "red") +
  labs(title = "Frechet-Hoeffding Upper Bounds",
       x = "Wind Speed (m/s)", y = "Bound Value") +
  theme_minimal()
copula_df <- data.frame(W = u, D = v)
bound_df <- data.frame(W = qweibull(upper_bound, shape = 2, scale = 10), Bounds = upper_bound * 360)
ggplot() +
  geom_point(data = copula_df, aes(x = W, y = D), alpha = 0.4) +
  geom_line(data = bound_df, aes(x = W, y = Bounds), color = "red") +
  labs(
    title = "Joint Distribution of Wind Speed and Direction",
    x = "Wind Speed (m/s)",
    y = "Wind Direction (degrees)"
  ) +
  theme_minimal()




#PROBLEM SECTION 3
library(MASS)  
library(copula)
library(ggplot2)
set.seed(123)
X <- rnorm(1000, mean = 70, sd = 10)
Y <- rnorm(1000, mean = 50, sd = 20)
plot(X, Y, main = "Scatterplot of Signal Strength vs Data Traffic", xlab = "Signal Strength (dBm)", ylab = "Data Traffic (MB/s)")
fit_X <- fitdistr(X, "normal")
fit_Y <- fitdistr(Y, "normal")
correlation <- cor(X, Y, method = "pearson")
print(correlation)
rho_spearman <- cor(X, Y, method = "spearman")
print(rho_spearman)
tau_kendall <- cor(X, Y, method = "kendall")
print(tau_kendall)
ecdf_X <- ecdf(X)
ecdf_Y <- ecdf(Y)
U <- ecdf_X(X)
V <- ecdf_Y(Y)
pseudo_sample <- cbind(U, V)
gauss_cop <- normalCopula(param = correlation, dim = 2)
fit_gauss_cop <- fitCopula(gauss_cop, pseudo_sample, method = "ml")
sample_cop_gauss <- rCopula(1000, fit_gauss_cop@copula)
sample_X_gauss <- qnorm(sample_cop_gauss[,1], mean = fit_X$estimate['mean'], sd = fit_X$estimate['sd'])
sample_Y_gauss <- qnorm(sample_cop_gauss[,2], mean = fit_Y$estimate['mean'], sd = fit_Y$estimate['sd'])
plot(sample_X_gauss, sample_Y_gauss, main = "Sample from Gaussian Copula (Pearson)", xlab = "Signal Strength (dBm)", ylab = "Data Traffic (MB/s)")
gauss_cop_spearman <- normalCopula(param = rho_spearman, dim = 2)
fit_gauss_cop_spearman <- fitCopula(gauss_cop_spearman, pseudo_sample, method = "ml")
sample_cop_spearman <- rCopula(1000, fit_gauss_cop_spearman@copula)
sample_X_spearman <- qnorm(sample_cop_spearman[,1], mean = fit_X$estimate['mean'], sd = fit_X$estimate['sd'])
sample_Y_spearman <- qnorm(sample_cop_spearman[,2], mean = fit_Y$estimate['mean'], sd = fit_Y$estimate['sd'])
plot(sample_X_spearman, sample_Y_spearman, main = "Sample from Gaussian Copula (Spearman)", xlab = "Signal Strength (dBm)", ylab = "Data Traffic (MB/s)")
gumbel_param <- 1 / (1 - tau_kendall)
gumbel_cop <- gumbelCopula(param = gumbel_param, dim = 2)
fit_gumbel_cop <- fitCopula(gumbel_cop, pseudo_sample, method = "ml")
sample_cop_gumbel <- rCopula(1000, fit_gumbel_cop@copula)
sample_X_gumbel <- qnorm(sample_cop_gumbel[,1], mean = fit_X$estimate['mean'], sd = fit_X$estimate['sd'])
sample_Y_gumbel <- qnorm(sample_cop_gumbel[,2], mean = fit_Y$estimate['mean'], sd = fit_Y$estimate['sd'])
plot(sample_X_gumbel, sample_Y_gumbel, main = "Sample from Gumbel Copula", xlab = "Signal Strength (dBm)", ylab = "Data Traffic (MB/s)")
param_kendall <- sin(pi * tau_kendall / 2)
gauss_cop_kendall <- normalCopula(param = param_kendall, dim = 2)
fit_gauss_cop_kendall <- fitCopula(gauss_cop_kendall, pseudo_sample, method = "ml")
sample_cop_kendall <- rCopula(1000, fit_gauss_cop_kendall@copula)
sample_X_kendall <- qnorm(sample_cop_kendall[,1], mean = fit_X$estimate['mean'], sd = fit_X$estimate['sd'])
sample_Y_kendall <- qnorm(sample_cop_kendall[,2], mean = fit_Y$estimate['mean'], sd = fit_Y$estimate['sd'])
plot(sample_X_kendall, sample_Y_kendall, main = "Sample from Gaussian Copula (Kendall)", xlab = "Signal Strength (dBm)", ylab = "Data Traffic (MB/s)")


#PROBLEM SECTION 4
library(quantmod)
library(copula)
library(rugarch)
library(tseries)
library(MASS)
library(ggplot2)
library(corrplot)
symbols <- c("AAPL", "JNJ", "XOM")  
getSymbols(symbols, src = "yahoo", from = "2010-01-01", to = "2024-05-31")
data <- na.omit(merge(Cl(AAPL), Cl(JNJ), Cl(XOM)))
returns <- na.omit(ROC(data, type = "discrete"))
chartSeries(data, theme = chartTheme("white"), name = "Adjusted Closing Prices")
chartSeries(returns, theme = chartTheme("white"), name = "Historical Returns") 
tau <- cor(returns, method = "kendall")
corrplot(tau, method = "color", col = colorRampPalette(c("blue", "white", "red"))(200), tl.col = "black", tl.srt = 45)
t_copula <- tCopula(dim = ncol(returns), df = 4)
fit <- fitCopula(t_copula, data = pobs(as.matrix(returns)), method = "ml")
set.seed(123)
simulated_copula <- rCopula(1000, fit@copula)
u <- seq(0, 1, length = 100)
upper_tau <- rep(NA, length(u))
lower_tau <- rep(NA, length(u))
for (i in 1:length(u)) {
  upper_tail <- simulated_copula[, 1] > u[i] & simulated_copula[, 2] > u[i]
  lower_tail <- simulated_copula[, 1] < u[i] & simulated_copula[, 2] < u[i]
  
  if (sum(upper_tail) > 1) {
    upper_tau[i] <- cor(simulated_copula[upper_tail, ], method = "kendall")[1, 2]
  }
  if (sum(lower_tail) > 1) {
    lower_tau[i] <- cor(simulated_copula[lower_tail, ], method = "kendall")[1, 2]
  }
}

ggplot(data.frame(u = u, upper_tau = upper_tau, lower_tau = lower_tau)) +
  geom_line(aes(x = u, y = upper_tau), color = "blue", linetype = "dashed", size = 1) +
  geom_line(aes(x = u, y = lower_tau), color = "red", linetype = "dashed", size = 1) +
  geom_line(aes(x = u, y = u), color = "black", linetype = "solid") +
  labs(title = "Tail Dependence Function (TDF) and Empirical Kendall's Tau",
       x = "u", y = "Kendall's Tau") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(-1, 1)) +
  theme_bw()


#PROBLEM SECTION 5
library(quantmod)
library(PerformanceAnalytics)
library(copula)
library(quadprog)
library(fPortfolio)
tickers <- c("AAPL", "JNJ", "XOM")
getSymbols(tickers, from = "2020-01-01", to = "2024-05-31", src = 'yahoo')
returns <- na.omit(merge(dailyReturn(AAPL), dailyReturn(JNJ), dailyReturn(XOM)))
colnames(returns) <- tickers
alpha <- 0.05
VaR_values <- apply(returns, 2, function(x) VaR(x, p = alpha, method = "historical"))
CVaR_values <- apply(returns, 2, function(x) CVaR(x, p = alpha, method = "historical"))
returns_matrix <- as.matrix(pobs(returns)) # Convert to matrix
copula_fit <- fitCopula(tCopula(dim = 3), returns_matrix, method = "ml")
n_sim <- 10000
copula_samples <- rCopula(n_sim, copula_fit@copula)
simulated_returns <- t(apply(copula_samples, 2, function(x) quantile(returns[, 1], x)))
mu <- colMeans(returns)
sigma <- cov(returns)
Dmat <- 2 * sigma
dvec <- rep(0, length(mu))
Amat <- cbind(1, diag(length(mu)))
bvec <- c(1, rep(0, length(mu)))
portfolio <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
optimal_weights <- portfolio$solution
cat("VaR at 95% confidence level:\n", VaR_values, "\n")
cat("CVaR at 95% confidence level:\n", CVaR_values, "\n")
cat("Optimal Portfolio Weights:\n", optimal_weights, "\n")