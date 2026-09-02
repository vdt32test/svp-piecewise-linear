# ============================================================
# L0 TREND FILTERING
# Sequential AMIAS implementation
# Based on Wen, Wang and Zhang (2023)
#
# q = 1 -> piecewise-linear trend filtering
# ============================================================

set.seed(123)


# ============================================================
# 1. SIMULATE THE SAME DATA AS IN THE L1 EXPERIMENT
# ============================================================

n <- 200
x <- seq(0, 1, length.out = n)

# True continuous piecewise-linear signal
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

# True slope-change locations
true_knots <- c(0.25, 0.50, 0.75)

# Add Gaussian noise
sigma <- 0.30
y <- f_true + rnorm(n, mean = 0, sd = sigma)


# ============================================================
# 2. BUILD THE DISCRETE DIFFERENCE MATRIX
# ============================================================

# Construct D^(q+1)
#
# q = 0 -> first differences
# q = 1 -> second differences
# q = 2 -> third differences
# etc.

make_difference_matrix <- function(n, q) {
  
  # First difference matrix
  D <- matrix(0, nrow = n - 1, ncol = n)
  
  for (i in 1:(n - 1)) {
    D[i, i]     <- -1
    D[i, i + 1] <-  1
  }
  
  if (q == 0) {
    return(D)
  }
  
  # Repeated differences
  for (order in 2:(q + 1)) {
    
    nr <- nrow(D)
    
    D_new <- matrix(
      0,
      nrow = nr - 1,
      ncol = n
    )
    
    for (i in 1:(nr - 1)) {
      D_new[i, ] <- D[i + 1, ] - D[i, ]
    }
    
    D <- D_new
  }
  
  return(D)
}


# q = 1 -> piecewise-linear model
q <- 1

D <- make_difference_matrix(n, q)

dim(D)

# For q = 1, this should be:
# 198 x 200


# ============================================================
# 3. HELPER FUNCTION:
# SOLVE FOR alpha, v, u GIVEN AN ACTIVE SET
# ============================================================

solve_given_active_set <- function(y, D, A) {
  
  m <- nrow(D)
  
  S <- seq_len(m)
  
  # Inactive set
  I <- setdiff(S, A)
  
  # Dual vector
  u <- rep(0, m)
  
  # If inactive set is nonempty:
  #
  # u_I = (D_I D_I^T)^(-1) D_I y
  if (length(I) > 0) {
    
    
    D_I <- D[I, , drop = FALSE]
    
    G <- D_I %*% t(D_I)
    
    rhs <- as.vector(D_I %*% y)
    
    u[I] <- as.vector(
      solve(G, rhs)
    )
  }
  
  # L0-TF estimate:
  #
  # alpha_hat = y - D_I^T u_I
  alpha_hat <- as.vector(
    y - t(D) %*% u
  )
  
  # Primal variable
  v <- as.vector(
    D %*% alpha_hat
  )
  
  # In theory:
  # v_I = 0
  #
  # Set these explicitly to zero to avoid numerical noise.
  if (length(I) > 0) {
    v[I] <- 0
  }
  
  # Active components of u are zero by construction
  if (length(A) > 0) {
    u[A] <- 0
  }
  
  list(
    alpha = alpha_hat,
    v = v,
    u = u,
    A = A,
    I = I
  )
}


# ============================================================
# 4. NAIVE AMIAS ALGORITHM
# ============================================================

amias_fixed_k <- function(
    y,
    D,
    k,
    A0,
    rho,
    max_iter = 20
) {
  
  m <- nrow(D)
  
  S <- seq_len(m)
  
  # Make sure active set has exactly k entries
  A <- sort(unique(A0))
  
  if (length(A) != k) {
    stop("Initial active set A0 must contain exactly k indices.")
  }
  
  converged <- FALSE
  
  for (iter in seq_len(max_iter)) {
    
    # ------------------------------------------
    # Step 1:
    # solve primal/dual variables
    # ------------------------------------------
    
    sol <- solve_given_active_set(
      y = y,
      D = D,
      A = A
    )
    
    v <- sol$v
    u <- sol$u
    
    # ------------------------------------------
    # Step 2:
    # primal-dual mixing variable
    #
    # xi = v + u/rho
    # ------------------------------------------
    
    xi <- v + u / rho
    
    # ------------------------------------------
    # Step 3:
    # hard thresholding:
    # choose top k components by magnitude
    # ------------------------------------------
    
    A_new <- sort(
      order(
        abs(xi),
        decreasing = TRUE
      )[seq_len(k)]
    )
    
    # ------------------------------------------
    # Step 4:
    # convergence check
    # ------------------------------------------
    
    if (identical(A_new, A)) {
      converged <- TRUE
      break
    }
    
    A <- A_new
  }
  
  # Recompute final solution using final active set
  final_sol <- solve_given_active_set(
    y = y,
    D = D,
    A = A
  )
  
  final_sol$converged <- converged
  final_sol$iterations <- iter
  final_sol$k <- k
  
  return(final_sol)
}


# ============================================================
# 5. SEQUENTIAL AMIAS
# ============================================================

sequential_amias <- function(
    y,
    D,
    q,
    kmax = 10,
    rho = NULL,
    max_iter = 20
) {
  
  n <- length(y)
  
  m <- nrow(D)
  
  S <- seq_len(m)
  
  # The paper recommends rho related to a power of n.
  #
  # For q = 1 this default is n^(q+1).
  if (is.null(rho)) {
    rho <- n^(q + 1)
  }
  
  results <- vector(
    "list",
    kmax
  )
  
  # ----------------------------------------------------------
  # k = 0 initial solution
  # ----------------------------------------------------------
  
  A_prev <- integer(0)
  
  sol0 <- solve_given_active_set(
    y = y,
    D = D,
    A = A_prev
  )
  
  u_prev <- sol0$u
  I_prev <- sol0$I
  
  # ----------------------------------------------------------
  # Sequentially fit k = 1,...,kmax
  # ----------------------------------------------------------
  
  for (k in 1:kmax) {
    
    # Warm-start rule:
    #
    # add the inactive index with the largest
    # magnitude dual variable
    
    candidate <- I_prev[
      which.max(abs(u_prev[I_prev]))
    ]
    
    A0 <- sort(
      c(A_prev, candidate)
    )
    
    # Run AMIAS for this fixed cardinality k
    fit <- amias_fixed_k(
      y = y,
      D = D,
      k = k,
      A0 = A0,
      rho = rho,
      max_iter = max_iter
    )
    
    # ------------------------------------------
    # Model-selection quantities
    # ------------------------------------------
    
    alpha_hat <- fit$alpha
    
    mse_data <- mean(
      (y - alpha_hat)^2
    )
    
    # Degrees of freedom from the paper:
    #
    # df = q + 1 + k
    
    df <- q + 1 + k
    
    # BIC used in the paper
    bic <- n * log(mse_data) +
      2 * df * log(n)
    
    fit$mse_data <- mse_data
    fit$df <- df
    fit$bic <- bic
    
    results[[k]] <- fit
    
    # Warm start for next cardinality
    A_prev <- fit$A
    I_prev <- fit$I
    u_prev <- fit$u
  }
  
  # ----------------------------------------------------------
  # Select k using minimum BIC
  # ----------------------------------------------------------
  
  bic_values <- sapply(
    results,
    function(z) z$bic
  )
  
  k_selected <- which.min(
    bic_values
  )
  
  list(
    fits = results,
    k_selected = k_selected,
    best_fit = results[[k_selected]],
    bic = bic_values,
    rho = rho
  )
}


# ============================================================
# 6. RUN L0 TREND FILTERING
# ============================================================

# Search models with between 1 and 10 knots
kmax <- 10

l0_result <- sequential_amias(
  y = y,
  D = D,
  q = q,
  kmax = kmax,
  max_iter = 20
)

# BIC-selected model
best_l0 <- l0_result$best_fit

fit_l0 <- best_l0$alpha

selected_k <- l0_result$k_selected


# ============================================================
# 7. CONVERT ACTIVE INDICES TO APPROXIMATE x LOCATIONS
# ============================================================

# For second differences, an active index j corresponds
# approximately to the middle observation x_(j+1).

active_to_x <- function(A, x, q) {
  
  if (length(A) == 0) {
    return(numeric(0))
  }
  
  if (q == 1) {
    return(x[A + 1])
  }
  
  # Generic approximate mapping for other q
  shift <- floor((q + 1) / 2)
  
  idx <- pmin(
    length(x),
    A + shift
  )
  
  x[idx]
}


estimated_knot_locations <- active_to_x(
  best_l0$A,
  x,
  q
)


# ============================================================
# 8. EVALUATION AGAINST THE KNOWN TRUE SIGNAL
# ============================================================

mse_true <- function(fitted, truth) {
  mean((fitted - truth)^2)
}

l0_true_mse <- mse_true(
  fit_l0,
  f_true
)

# Clean summary table
results <- data.frame(
  Selected_knots = selected_k,
  True_knots = length(true_knots),
  MSE = round(l0_true_mse, 5),
  BIC = round(best_l0$bic, 2)
)

print(results)

# ============================================================
# 9. TABLE FOR ALL CARDINALITIES k
# ============================================================

summary_table <- data.frame(
  k = 1:kmax,
  MSE_data = sapply(
    l0_result$fits,
    function(z) z$mse_data
  ),
  MSE_true = sapply(
    l0_result$fits,
    function(z)
      mse_true(z$alpha, f_true)
  ),
  BIC = sapply(
    l0_result$fits,
    function(z) z$bic
  ),
  Iterations = sapply(
    l0_result$fits,
    function(z) z$iterations
  ),
  Converged = sapply(
    l0_result$fits,
    function(z) z$converged
  )
)

print(summary_table)


# ============================================================
# 10. PLOT: FINAL BIC-SELECTED L0 FIT
# ============================================================

ylim_plot <- range(
  c(y, f_true, fit_l0)
)

plot(
  x,
  y,
  pch = 16,
  cex = 0.45,
  ylim = ylim_plot,
  xlab = "x",
  ylab = "y",
  main = paste0(
    "L0 trend filtering: selected k = ",
    selected_k
  )
)

# True function
lines(
  x,
  f_true,
  col = "blue",
  lwd = 3
)

# L0 trend-filtering estimate
lines(
  x,
  fit_l0,
  col = "red",
  lwd = 2
)

# Estimated knots
for (z in estimated_knot_locations) {
  
  abline(
    v = z,
    col = "red",
    lty = 3
  )
}

# True knots
for (z in true_knots) {
  
  abline(
    v = z,
    col = "blue",
    lty = 3
  )
}

legend(
  "topleft",
  legend = c(
    "Observations",
    "True function",
    "L0 trend filtering",
    "True knots",
    "Estimated knots"
  ),
  col = c(
    "black",
    "blue",
    "red",
    "blue",
    "red"
  ),
  pch = c(
    16,
    NA,
    NA,
    NA,
    NA
  ),
  lty = c(
    NA,
    1,
    1,
    3,
    3
  ),
  lwd = c(
    NA,
    3,
    2,
    1,
    1
  ),
  bty = "n",
  cex = 0.75
)


# ============================================================
# 11. PLOT BIC AGAINST NUMBER OF KNOTS
# ============================================================

plot(
  1:kmax,
  l0_result$bic,
  type = "b",
  pch = 16,
  xlab = "Number of knots k",
  ylab = "BIC",
  main = "BIC model selection"
)

abline(
  v = selected_k,
  lty = 2
)


# ============================================================
# 12. OPTIONAL:
# COMPARE DIFFERENT VALUES OF k
# ============================================================

k_values <- c(1, 3, 6)

par(
  mfrow = c(1, 3)
)

for (k in k_values) {
  
  current_fit <-
    l0_result$fits[[k]]$alpha
  
  current_knots <-
    active_to_x(
      l0_result$fits[[k]]$A,
      x,
      q
    )
  
  plot(
    x,
    y,
    pch = 16,
    cex = 0.35,
    ylim = ylim_plot,
    xlab = "x",
    ylab = "y",
    main = paste("k =", k)
  )
  
  lines(
    x,
    f_true,
    col = "blue",
    lwd = 3
  )
  
  lines(
    x,
    current_fit,
    col = "red",
    lwd = 2
  )
  
  for (z in current_knots) {
    
    abline(
      v = z,
      col = "red",
      lty = 3
    )
  }
}

par(
  mfrow = c(1, 1)
)