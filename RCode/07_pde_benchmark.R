#-------------------------------------------------------------------------------

##  07_pde_benchmark.R
##  A deterministic benchmark for the reduced XVA BSDE.
##
##  Pang and Sturm report no benchmark for their deep BSDE figures, so the
##  accuracy of Figures 4 and 5 cannot be assessed from the paper alone.  This
##  module supplies one.  The Markov structure of the reduced BSDE (17) is
##  explicit: the driver depends on the path only through (t, S_t, beta_t), so
##  the nonlinear Feynman-Kac correspondence gives.

#------------------------------------------------------------------------------

## First derivative by central differences, one-sided at the two end points.
d_dx <- function(u, dx) {
  n <- length(u)
  g <- numeric(n)
  g[2:(n - 1)] <- (u[3:n] - u[1:(n - 2)]) / (2 * dx)
  g[1] <- (u[2] - u[1]) / dx
  g[n] <- (u[n] - u[n - 1]) / dx
  g
}

# -----------------------------------------------------------------------------

## pde_operator(): builds the log-price grid and inverts the implicit operator
## once.

# -----------------------------------------------------------------------------
pde_operator <- function(p, nx = 801L, nt = 1000L, xwidth = 6) {

  sdev <- p$sigma * sqrt(p$Tmat)
  x    <- seq(log(p$S0) - xwidth * sdev, log(p$S0) + xwidth * sdev, length.out = nx)
  dx   <- x[2] - x[1]
  dt   <- p$Tmat / nt

  dif <- 0.5 * p$sigma^2
  adv <- p$rD - 0.5 * p$sigma^2
  lo  <- dif / dx^2 - adv / (2 * dx)
  up  <- dif / dx^2 + adv / (2 * dx)
  ce  <- -2 * dif / dx^2

  A <- matrix(0, nrow = nx, ncol = nx)
  for (i in 2:(nx - 1)) {
    A[i, i - 1] <- -dt * lo
    A[i, i]     <- 1 - dt * ce
    A[i, i + 1] <- -dt * up
  }
  A[1, 1]   <- 1
  A[nx, nx] <- 1

  ## The reference pair at every node, column n holding t_n = n dt.  It
  ## depends only on the grid and on (S0, K, T, sigma, r_D), not on alpha, the
  ## funding rates, the regime or the side, so a sweep that reuses the operator
  ## reuses it too.  The key records what it was computed for.
  S  <- exp(x)
  tn <- seq_len(nt) * dt
  Vh <- vapply(tn, function(t) bs_value(t, S, p), numeric(nx))
  Zh <- vapply(tn, function(t) bs_control(t, S, p), numeric(nx))

  list(x = x, S = S, dx = dx, dt = dt, nx = nx, nt = nt, Ainv = solve(A),
       Vh = Vh, Zh = Zh, ref_key = pde_ref_key(p))
}

pde_ref_key <- function(p) paste(p$S0, p$K, p$Tmat, p$sigma, p$rD, sep = "|")

# -----------------------------------------------------------------------------

## pde_solve(): returns the XVA surface and its value at (0, S0) for each
## regime carried by p$beta_mode.  With keep = TRUE the solution at every time
## step is returned as well: hist[[regime]] is an nx x (nt + 1) matrix whose
## column n + 1 holds u(t_n, .), used by the replication simulator of 18.

# -----------------------------------------------------------------------------
pde_solve <- function(p, op = NULL, nx = 801L, nt = 1000L, xwidth = 6,
                      keep = FALSE) {

  if (is.null(op)) op <- pde_operator(p, nx, nt, xwidth)
  gfun <- get_driver(p)

  dynamic <- identical(p$beta_mode, "dynamic")
  regimes <- if (dynamic) c(0L, 1L) else
             if (identical(p$beta_mode, "frozen1")) 1L else 0L

  U <- lapply(regimes, function(r) numeric(op$nx))    # terminal condition is 0
  names(U) <- as.character(regimes)
  cached <- identical(op$ref_key, pde_ref_key(p))
  RHS <- matrix(0, op$nx, length(regimes))
  hist <- if (keep) lapply(U, function(u) matrix(0, op$nx, op$nt + 1L)) else NULL

  for (n in op$nt:1) {
    tn <- n * op$dt
    Vh <- if (cached) op$Vh[, n] else bs_value(tn, op$S, p)
    Zh <- if (cached) op$Zh[, n] else bs_control(tn, op$S, p)

    for (j in seq_along(regimes)) {
      r   <- regimes[j]
      u   <- U[[as.character(r)]]
      z   <- p$sigma * d_dx(u, op$dx)                 # Zbreve = sigma du/dx
      src <- gfun(u, z, rep(r, op$nx), Vh, Zh, p) * op$dt
      if (dynamic) {
        rate  <- if (r == 0L) p$aU else p$aV
        other <- U[[as.character(1L - r)]]
        src   <- src + op$dt * rate * (other - u)
      }
      rhs <- u + src
      rhs[1] <- rhs[2]; rhs[op$nx] <- rhs[op$nx - 1]  # extrapolation rows
      RHS[, j] <- rhs
    }
    V <- op$Ainv %*% RHS                              # both regimes in one product
    for (j in seq_along(regimes)) {
      v <- V[, j]
      v[1]      <- 2 * v[2] - v[3]
      v[op$nx]  <- 2 * v[op$nx - 1] - v[op$nx - 2]
      U[[as.character(regimes[j])]] <- v
      if (keep) hist[[as.character(regimes[j])]][, n] <- v
    }
  }

  x0  <- log(p$S0)
  val <- sapply(regimes, function(r) approx(op$x, U[[as.character(r)]], x0)$y)
  names(val) <- as.character(regimes)
  list(value = val, U = U, x = op$x, S = op$S, regimes = regimes, hist = hist,
       dt = op$dt, nt = op$nt)
}

## Convenience wrapper returning a single number for the regime in force.
pde_value <- function(p, op = NULL, ...) {
  res <- pde_solve(p, op, ...)
  key <- if (identical(p$beta_mode, "frozen1")) "1" else "0"
  as.numeric(res$value[[key]])
}
