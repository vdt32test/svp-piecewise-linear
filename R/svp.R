# ============================================================
# SIMPLE EXAMPLE: SMALLEST VALID PARTITIONING (SVP)
# ============================================================

set.seed(123)

# ------------------------------------------------------------
# 1. Simulate a piecewise-constant signal
# ------------------------------------------------------------

n <- 300
x <- 1:n

true_cpts <- c(80, 150, 240)

f_true <- c(
  rep(0.0, 80),
  rep(2.0, 70),
  rep(-1.0, 90),
  rep(1.5, 60)
)

sigma <- 0.35
y <- f_true + rnorm(n, mean = 0, sd = sigma)

# ------------------------------------------------------------
# 2. Precompute cumulative sums
# ------------------------------------------------------------

cum_y  <- c(0, cumsum(y))
cum_y2 <- c(0, cumsum(y^2))


# ------------------------------------------------------------
# 3. Gaussian segment cost
#
# Segment (a,b] contains:
# y[a+1], ..., y[b]
#
# C = 1/2 * sum (y_i - mean)^2
# ------------------------------------------------------------

segment_cost <- function(a, b) {
  
  len <- b - a
  
  sum_y  <- cum_y[b + 1]  - cum_y[a + 1]
  sum_y2 <- cum_y2[b + 1] - cum_y2[a + 1]
  
  0.5 * (
    sum_y2 -
      sum_y^2 / len
  )
}


# ------------------------------------------------------------
# 4. SVP validity statistic
#
# f_OP(a,b) =
# max_u [
#     C(a,b) - C(a,u) - C(u,b)
# ]
#
# This measures the largest improvement obtained by
# introducing one changepoint inside the segment.
# ------------------------------------------------------------

validity_statistic <- function(a, b) {
  
  # A segment containing only one observation cannot
  # contain an internal changepoint.
  if (b - a <= 1) {
    return(0)
  }
  
  cost_full <- segment_cost(a, b)
  
  best_improvement <- 0
  
  for (u in (a + 1):(b - 1)) {
    
    improvement <-
      cost_full -
      segment_cost(a, u) -
      segment_cost(u, b)
    
    if (improvement > best_improvement) {
      best_improvement <- improvement
    }
  }
  
  best_improvement
}


# ------------------------------------------------------------
# 5. Choose the validity threshold gamma
# ------------------------------------------------------------

gamma <- 2 * log(n)

# Segment is valid iff:
# validity_statistic <= gamma


# ------------------------------------------------------------
# 6. SVP dynamic programming
#
# For every endpoint t we store:
#
# K[t]    = smallest number of valid segments
# Q[t]    = smallest total cost among solutions with K[t]
# prev[t] = beginning of final segment
#
# Lexicographic comparison:
#
# first minimize K,
# then minimize Q.
# ------------------------------------------------------------

K <- rep(Inf, n + 1)
Q <- rep(Inf, n + 1)
previous <- rep(NA_integer_, n + 1)

# Prefix containing zero observations
K[1] <- 0
Q[1] <- 0


for (t in 1:n) {
  
  for (s in 0:(t - 1)) {
    
    f_value <- validity_statistic(s, t)
    
    # Only consider valid segments
    if (f_value <= gamma) {
      
      candidate_K <- K[s + 1] + 1
      
      candidate_Q <-
        Q[s + 1] +
        segment_cost(s, t)
      
      # Lexicographic minimization
      if (
        candidate_K < K[t + 1] ||
        (
          candidate_K == K[t + 1] &&
          candidate_Q < Q[t + 1]
        )
      ) {
        
        K[t + 1] <- candidate_K
        Q[t + 1] <- candidate_Q
        previous[t + 1] <- s
      }
    }
  }
}


# ------------------------------------------------------------
# 7. Backtracking
# ------------------------------------------------------------

estimated_cpts <- integer(0)

t <- n

while (t > 0) {
  
  s <- previous[t + 1]
  
  if (s > 0) {
    estimated_cpts <- c(s, estimated_cpts)
  }
  
  t <- s
}


# ------------------------------------------------------------
# 8. Construct fitted piecewise-constant signal
# ------------------------------------------------------------

boundaries <- c(
  0,
  estimated_cpts,
  n
)

fit_svp <- numeric(n)

for (j in 1:(length(boundaries) - 1)) {
  
  a <- boundaries[j]
  b <- boundaries[j + 1]
  
  segment_mean <- mean(
    y[(a + 1):b]
  )
  
  fit_svp[(a + 1):b] <- segment_mean
}


# ------------------------------------------------------------
# 9. Evaluation
# ------------------------------------------------------------

mse <- mean(
  (fit_svp - f_true)^2
)

results <- data.frame(
  Quantity = c(
    "Gamma",
    "True number of changepoints",
    "Estimated number of changepoints",
    "MSE"
  ),
  Result = c(
    round(gamma, 4),
    length(true_cpts),
    length(estimated_cpts),
    round(mse, 6)
  )
)

print(
  results,
  row.names = FALSE
)

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
# 10. Plot
# ------------------------------------------------------------

plot(
  x,
  y,
  pch = 16,
  cex = 0.45,
  xlab = "Time",
  ylab = "y",
  main = "Smallest Valid Partitioning"
)

lines(
  x,
  f_true,
  lwd = 3
)

lines(
  x,
  fit_svp,
  lwd = 2,
  lty = 2
)

abline(
  v = true_cpts,
  lty = 3
)

if (length(estimated_cpts) > 0) {
  abline(
    v = estimated_cpts,
    lty = 2
  )
}

legend(
  "topright",
  legend = c(
    "Observations",
    "True signal",
    "SVP estimate",
    "True changepoints",
    "Estimated changepoints"
  ),
  pch = c(16, NA, NA, NA, NA),
  lty = c(NA, 1, 2, 3, 2),
  lwd = c(NA, 3, 2, 1, 1),
  cex = 0.8
)