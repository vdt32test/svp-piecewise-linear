# ============================================================
# OPTIMAL PARTITIONING / PELT EXPERIMENT
# Piecewise-constant Gaussian mean changes
# ============================================================

library(changepoint)

set.seed(123)

# ------------------------------------------------------------
# 1. Simulate data
# ------------------------------------------------------------

n <- 300
x <- 1:n

# True piecewise-constant signal
f_true <- c(
  rep(0.0,  80),
  rep(2.0,  70),
  rep(-1.0, 90),
  rep(1.5,  60)
)

# True changepoints are after observations 80, 150, 240
true_cpts <- c(80, 150, 240)

sigma <- 0.7
y <- f_true + rnorm(n, mean = 0, sd = sigma)


# ------------------------------------------------------------
# 2. Fit PELT
# ------------------------------------------------------------

# Detect changes in the mean
# method = "PELT" uses Pruned Exact Linear Time
#
# penalty = "BIC" chooses a standard complexity penalty automatically

fit_pelt <- cpt.mean(
  y,
  method = "PELT",
  penalty = "BIC",
  class = TRUE
)

# Estimated changepoints
estimated_cpts <- cpts(fit_pelt)

# Remove n if the package includes it as the final endpoint
estimated_cpts <- estimated_cpts[estimated_cpts < n]


# ------------------------------------------------------------
# 3. Recover fitted piecewise-constant signal
# ------------------------------------------------------------

segment_ends <- c(estimated_cpts, n)
segment_starts <- c(1, estimated_cpts + 1)

fit_signal <- numeric(n)

for (j in seq_along(segment_starts)) {
  
  a <- segment_starts[j]
  b <- segment_ends[j]
  
  fit_signal[a:b] <- mean(y[a:b])
}


# ------------------------------------------------------------
# 4. Evaluation
# ------------------------------------------------------------

mse <- function(fitted, truth) {
  mean((fitted - truth)^2)
}

mse_value <- mse(
  fit_signal,
  f_true
)

results <- data.frame(
  True_changepoints = length(true_cpts),
  Estimated_changepoints = length(estimated_cpts),
  MSE = round(mse_value, 5)
)

print(results)

cat(
  "\nTrue changepoint locations:      ",
  paste(true_cpts, collapse = ", "),
  "\n"
)

cat(
  "Estimated changepoint locations: ",
  paste(estimated_cpts, collapse = ", "),
  "\n"
)


# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------

ylim_plot <- range(
  c(y, f_true, fit_signal)
)

plot(
  x,
  y,
  pch = 16,
  cex = 0.45,
  ylim = ylim_plot,
  xlab = "Time index",
  ylab = "y",
  main = "PELT changepoint detection"
)

# True signal
lines(
  x,
  f_true,
  lwd = 3
)

# PELT fitted signal
lines(
  x,
  fit_signal,
  lwd = 2,
  lty = 2
)

# True changepoints
for (cp in true_cpts) {
  abline(
    v = cp,
    lty = 3
  )
}

# Estimated changepoints
for (cp in estimated_cpts) {
  abline(
    v = cp,
    lty = 2
  )
}

legend(
  "topright",
  legend = c(
    "Observations",
    "True signal",
    "PELT estimate",
    "True changepoints",
    "Estimated changepoints"
  ),
  pch = c(16, NA, NA, NA, NA),
  lty = c(NA, 1, 2, 3, 2),
  lwd = c(NA, 3, 2, 1, 1),
  bty = "n",
  cex = 0.8
)