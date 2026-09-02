# Install once if needed:
# install.packages("genlasso")

library(genlasso)

set.seed(123)

# -----------------------------
# 1. Simulate data
# -----------------------------

n <- 200
x <- seq(0, 1, length.out = n)

f_true <- ifelse(
  x < 0.25,
  2 * x,
  ifelse(
    x < 0.50,
    0.5 + 0.5 * (x - 0.25),
    ifelse(
      x < 0.75,
      0.625 + 3 * (x - 0.50),
      1.375 - 1.5 * (x - 0.75)
    )
  )
)

sigma <- 0.30
y <- f_true + rnorm(n, sd = sigma)

# -----------------------------
# 2. Trend filtering
# -----------------------------

# ord = 1 -> piecewise-linear trend filtering
# pos = x makes the design explicit
tf_fit <- trendfilter(y, pos = x, ord = 1)

# Inspect the range of lambda values
range(tf_fit$lambda)

lambda_max <- max(tf_fit$lambda)

# Choose clearly separated regularization levels
lambda_large  <- 0.8 * lambda_max
lambda_medium <- 0.15 * lambda_max
lambda_small  <- 0.01 * lambda_max

# Check ordering explicitly
stopifnot(lambda_small < lambda_medium,
          lambda_medium < lambda_large)

fit_small <- as.vector(
  predict(tf_fit, lambda = lambda_small)$fit
)

fit_medium <- as.vector(
  predict(tf_fit, lambda = lambda_medium)$fit
)

fit_large <- as.vector(
  predict(tf_fit, lambda = lambda_large)$fit
)

# -----------------------------
# 3. Evaluation functions
# -----------------------------

mse <- function(fitted, truth) {
  mean((fitted - truth)^2)
}

# For ord = 1, slope changes correspond to non-zero
# second differences.
#
# Use a relative tolerance instead of 1e-6, because
# numerical solutions are not exactly zero.
count_knots <- function(fitted, tol = 1e-4) {
  d2 <- diff(fitted, differences = 2)
  
  threshold <- tol * max(1, max(abs(fitted)))
  
  sum(abs(d2) > threshold)
}

# -----------------------------
# 4. Numerical summary
# -----------------------------

results <- data.frame(
  Model = c(
    "Small lambda",
    "Medium lambda",
    "Large lambda"
  ),
  Lambda = c(
    lambda_small,
    lambda_medium,
    lambda_large
  ),
  Knots = c(
    count_knots(fit_small),
    count_knots(fit_medium),
    count_knots(fit_large)
  ),
  MSE = c(
    mse(fit_small, f_true),
    mse(fit_medium, f_true),
    mse(fit_large, f_true)
  )
)

print(results)

# -----------------------------
# 5. Plot
# -----------------------------

ylim <- range(c(y, f_true, fit_small, fit_medium, fit_large))

par(mfrow = c(2, 2))

# ============================================================
# 1. Simulated data
# ============================================================

plot(
  x, y,
  pch = 16,
  cex = 0.45,
  ylim = ylim,
  main = "Simulated data",
  xlab = "x",
  ylab = "y"
)

# True function
lines(x, f_true, col = "blue", lwd = 3)

legend(
  "topleft",
  legend = c("Observations", "True function"),
  col = c("black", "blue"),
  pch = c(16, NA),
  lty = c(NA, 1),
  lwd = c(NA, 3),
  bty = "n",
  cex = 0.8
)


# ============================================================
# 2. Small lambda
# ============================================================

plot(
  x, y,
  pch = 16,
  cex = 0.4,
  ylim = ylim,
  main = paste0(
    "Small lambda = ",
    signif(lambda_small, 3)
  ),
  xlab = "x",
  ylab = "y"
)

# True function
lines(x, f_true, col = "blue", lwd = 3)

# Trend filtering estimate
lines(x, fit_small, col = "red", lwd = 2)

legend(
  "topleft",
  legend = c(
    "Observations",
    "True function",
    "Trend filtering"
  ),
  col = c("black", "blue", "red"),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = c(NA, 3, 2),
  bty = "n",
  cex = 0.8
)


# ============================================================
# 3. Medium lambda
# ============================================================

plot(
  x, y,
  pch = 16,
  cex = 0.4,
  ylim = ylim,
  main = paste0(
    "Medium lambda = ",
    signif(lambda_medium, 3)
  ),
  xlab = "x",
  ylab = "y"
)

# Same true function
lines(x, f_true, col = "blue", lwd = 3)

# Different trend filtering estimate
lines(x, fit_medium, col = "red", lwd = 2)

legend(
  "topleft",
  legend = c(
    "Observations",
    "True function",
    "Trend filtering"
  ),
  col = c("black", "blue", "red"),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = c(NA, 3, 2),
  bty = "n",
  cex = 0.8
)


# ============================================================
# 4. Large lambda
# ============================================================

plot(
  x, y,
  pch = 16,
  cex = 0.4,
  ylim = ylim,
  main = paste0(
    "Large lambda = ",
    signif(lambda_large, 3)
  ),
  xlab = "x",
  ylab = "y"
)

# Same true function
lines(x, f_true, col = "blue", lwd = 3)

# Different trend filtering estimate
lines(x, fit_large, col = "red", lwd = 2)

legend(
  "topleft",
  legend = c(
    "Observations",
    "True function",
    "Trend filtering"
  ),
  col = c("black", "blue", "red"),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = c(NA, 3, 2),
  bty = "n",
  cex = 0.8
)

par(mfrow = c(1, 1))
