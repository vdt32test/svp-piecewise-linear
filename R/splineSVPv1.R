# ============================================================
# Spline-SVP v1: exhaustive oracle for tiny examples
# ============================================================


# ------------------------------------------------------------
# 1. Validity statistic T_{s,t}
#
# R indices s,t correspond directly to mathematical indices.
# The candidate segment is (s,t], so observations s+1,...,t
# are used in the local validity statistic.
# ------------------------------------------------------------

segment_validity_stat <- function(x, y, s, t, sigma) {
  
  stopifnot(
    length(x) == length(y),
    1 <= s, s < t, t <= length(x),
    sigma > 0
  )
  
  idx <- (s + 1):t
  
  # No possible internal knot
  if (t - s <= 1) {
    return(0)
  }
  
  xs <- x[idx]
  ys <- y[idx]
  
  # Null model:
  # y_i = a + b x_i + error_i
  X0 <- cbind(1, xs)
  
  fit0 <- lm.fit(x = X0, y = ys)
  rss0 <- sum(fit0$residuals^2)
  
  # Candidate internal knot indices
  candidate_r <- (s + 1):(t - 1)
  
  best_rss1 <- Inf
  
  for (r in candidate_r) {
    
    # Continuous one-knot alternative:
    #
    # a + b*x + c*(x - x_r)_+
    hinge <- pmax(xs - x[r], 0)
    
    X1 <- cbind(1, xs, hinge)
    
    fit1 <- lm.fit(x = X1, y = ys)
    rss1 <- sum(fit1$residuals^2)
    
    if (rss1 < best_rss1) {
      best_rss1 <- rss1
    }
  }
  
  # Numerical protection: theoretically rss1 <= rss0
  improvement <- max(0, rss0 - best_rss1)
  
  improvement / sigma^2
}


# ------------------------------------------------------------
# 2. Precompute all T_{s,t}
# ------------------------------------------------------------

compute_T_matrix <- function(x, y, sigma) {
  
  n <- length(y)
  
  Tmat <- matrix(
    NA_real_,
    nrow = n,
    ncol = n,
    dimnames = list(1:n, 1:n)
  )
  
  for (s in 1:(n - 1)) {
    for (t in (s + 1):n) {
      Tmat[s, t] <- segment_validity_stat(
        x, y, s, t, sigma
      )
    }
  }
  
  Tmat
}


# ------------------------------------------------------------
# 3. Construct the continuous piecewise-linear design matrix
#    for a fixed partition tau.
#
# tau = c(tau_0, ..., tau_K)
#
# Columns correspond to boundary states
# z_0, ..., z_K.
#
# Each observation is represented by linear interpolation
# between the two adjacent boundary states.
# ------------------------------------------------------------

spline_design_matrix <- function(x, tau) {
  
  n <- length(x)
  K <- length(tau) - 1
  
  X <- matrix(0, nrow = n, ncol = K + 1)
  
  # First observation corresponds to z_0
  X[1, 1] <- 1
  
  for (k in 1:K) {
    
    s <- tau[k]
    t <- tau[k + 1]
    
    idx <- (s + 1):t
    
    w <- (x[idx] - x[s]) / (x[t] - x[s])
    
    # l_{s,t}^{u,v}(x_i)
    # = (1-w_i) u + w_i v
    
    X[idx, k]     <- 1 - w
    X[idx, k + 1] <- w
  }
  
  X
}


# ------------------------------------------------------------
# 4. Solve the continuous least-squares problem for one
#    fixed partition.
#
# Returns z*(tau) and C*(tau).
# ------------------------------------------------------------

fit_fixed_partition <- function(x, y, tau) {
  
  X <- spline_design_matrix(x, tau)
  
  # QR-based least squares is preferable to explicitly computing
  # solve(t(X) %*% X) %*% t(X) %*% y.
  fit <- lm.fit(x = X, y = y)
  
  z_star <- as.numeric(fit$coefficients)
  fitted <- as.numeric(X %*% z_star)
  residuals <- y - fitted
  
  C_star <- sum(residuals^2)
  
  list(
    tau = tau,
    K = length(tau) - 1,
    z = z_star,
    cost = C_star,
    fitted = fitted,
    residuals = residuals
  )
}


# ------------------------------------------------------------
# 5. Enumerate all partitions
#
# Endpoints 1 and n are always included.
# Each index 2,...,n-1 is either an interior boundary or not.
#
# Number of partitions = 2^(n-2).
# ------------------------------------------------------------

all_partitions <- function(n) {
  
  stopifnot(n >= 2)
  
  interior <- if (n > 2) 2:(n - 1) else integer(0)
  m <- length(interior)
  
  partitions <- vector("list", 2^m)
  
  for (mask in 0:(2^m - 1)) {
    
    if (m == 0) {
      chosen <- integer(0)
    } else {
      bits <- as.logical(
        intToBits(mask)[seq_len(m)]
      )
      
      chosen <- interior[bits]
    }
    
    partitions[[mask + 1]] <- c(1, chosen, n)
  }
  
  partitions
}


# ------------------------------------------------------------
# 6. Check whether a partition belongs to V_gamma
# ------------------------------------------------------------

partition_is_valid <- function(tau, Tmat, gamma) {
  
  K <- length(tau) - 1
  
  for (k in 1:K) {
    
    s <- tau[k]
    t <- tau[k + 1]
    
    if (Tmat[s, t] > gamma) {
      return(FALSE)
    }
  }
  
  TRUE
}


# ------------------------------------------------------------
# 7. Exhaustive spline-SVP oracle
#
# Primary objective:
#       minimise K
#
# Secondary objective:
#       minimise C*(tau)
# ------------------------------------------------------------

spline_svp_oracle <- function(x, y, sigma, gamma,
                              verbose = FALSE) {
  
  n <- length(y)
  
  stopifnot(
    length(x) == n,
    n >= 2,
    all(diff(x) > 0),
    sigma > 0,
    gamma >= 0
  )
  
  # Precompute state-independent validity
  Tmat <- compute_T_matrix(x, y, sigma)
  
  # Enumerate all candidate partitions
  partitions <- all_partitions(n)
  
  best <- NULL
  number_valid <- 0L
  
  for (tau in partitions) {
    
    # Reject invalid partitions immediately
    if (!partition_is_valid(tau, Tmat, gamma)) {
      next
    }
    
    number_valid <- number_valid + 1L
    
    candidate <- fit_fixed_partition(x, y, tau)
    
    if (verbose) {
      cat(
        "tau =", paste(tau, collapse = ","),
        " | K =", candidate$K,
        " | cost =", candidate$cost,
        "\n"
      )
    }
    
    # Lexicographic comparison:
    # first K, then cost
    if (
      is.null(best) ||
      candidate$K < best$K ||
      (
        candidate$K == best$K &&
        candidate$cost < best$cost
      )
    ) {
      best <- candidate
    }
  }
  
  if (is.null(best)) {
    stop("No valid partition found.")
  }
  
  best$T <- Tmat
  best$gamma <- gamma
  best$sigma <- sigma
  best$number_partitions <- length(partitions)
  best$number_valid_partitions <- number_valid
  
  best
}


# ------------------------------------------------------------
# 8. Evaluate fitted spline at arbitrary x values
# ------------------------------------------------------------

predict_spline <- function(object, x, x_new) {
  
  tau <- object$tau
  z <- object$z
  K <- object$K
  
  out <- numeric(length(x_new))
  
  for (j in seq_along(x_new)) {
    
    xx <- x_new[j]
    
    if (xx < x[1] || xx > x[length(x)]) {
      out[j] <- NA_real_
      next
    }
    
    # Find segment containing xx
    k <- which(
      xx >= x[tau[1:K]] &
        xx <= x[tau[2:(K + 1)]]
    )[1]
    
    s <- tau[k]
    t <- tau[k + 1]
    
    w <- (xx - x[s]) / (x[t] - x[s])
    
    out[j] <- (1 - w) * z[k] + w * z[k + 1]
  }
  
  out
}


# ============================================================
# Tiny example
# ============================================================

set.seed(1)

n <- 10
x <- 1:n
sigma <- 0.25

# Continuous PWL true signal:
# slope changes around x = 5
f0 <- ifelse(
  x <= 5,
  0.5 * x,
  2.5 - 0.7 * (x - 5)
)

y <- f0 + rnorm(n, sd = sigma)

# Try a threshold
gamma <- 4

result <- spline_svp_oracle(
  x = x,
  y = y,
  sigma = sigma,
  gamma = gamma,
  verbose = TRUE
)

cat("\nSelected partition:\n")
print(result$tau)

cat("\nNumber of segments K:\n")
print(result$K)

cat("\nOptimal boundary states:\n")
print(result$z)

cat("\nOptimal cost C*(tau):\n")
print(result$cost)

cat("\nPartitions checked:\n")
print(result$number_partitions)

cat("\nValid partitions:\n")
print(result$number_valid_partitions)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

grid <- seq(min(x), max(x), length.out = 500)

fhat_grid <- predict_spline(
  result,
  x = x,
  x_new = grid
)

plot(
  x, y,
  pch = 19,
  xlab = "x",
  ylab = "y",
  main = "Exhaustive spline-SVP v1 oracle"
)

lines(
  x, f0,
  lty = 2,
  lwd = 2
)

lines(
  grid, fhat_grid,
  lwd = 2
)

points(
  x[result$tau],
  result$z,
  pch = 4,
  cex = 1.4,
  lwd = 2
)

legend(
  "topright",
  legend = c(
    "observations",
    "true signal",
    "spline-SVP fit",
    "fitted boundary states"
  ),
  pch = c(19, NA, NA, 4),
  lty = c(NA, 2, 1, NA),
  lwd = c(NA, 2, 2, 2)
)