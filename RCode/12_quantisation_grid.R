#-------------------------------------------------------------------------------

##  12_quantisation_grid.R
##  Recursive marginal quantisation of the forward process.
##
##  Sources: Pages and Sagna (2015) for the recursion; McWalter, Rudd, Kienitz
##  and Platen (2018) for the Newton-Raphson form, the gradient and the Hessian
##  (their Equations (5)-(14) and Appendix A); Rudd et al., "Robust product
##  Markovian quantization" (2020) for the Lloyd fallback.
##
##  WHAT IS QUANTISED.  Only the stock.  The status process already lives on the
##  two points {0, 1}, and because W is independent of beta the transition
##  kernel of the pair (S, beta) is the product of the two kernels -- exactly,
##  with no approximation in beta at all (13_quantisation_solver.R).

#------------------------------------------------------------------------------
std_normal_guess <- function(N) {                  # McWalter et al., A.1.1
  n <- seq_len(N); 5.5 * n / (N + 1) - 2.75
}

## Thomas algorithm for a tridiagonal system; returns NULL on a vanishing pivot.
tridiag_solve <- function(lower, diag, upper, rhs) {
  n  <- length(diag)
  cp <- numeric(n); dp <- numeric(n)
  if (!is.finite(diag[1]) || abs(diag[1]) < 1e-300) return(NULL)
  cp[1] <- if (n > 1) upper[1] / diag[1] else 0
  dp[1] <- rhs[1] / diag[1]
  if (n > 1) for (i in 2:n) {
    den <- diag[i] - lower[i - 1] * cp[i - 1]
    if (!is.finite(den) || abs(den) < 1e-300) return(NULL)
    cp[i] <- if (i < n) upper[i] / den else 0
    dp[i] <- (rhs[i] - lower[i - 1] * dp[i - 1]) / den
  }
  x <- numeric(n); x[n] <- dp[n]
  if (n > 1) for (i in (n - 1):1) x[i] <- dp[i] - cp[i] * x[i + 1]
  x
}

## Standardised region boundaries of every mixture component against grid g.
## Cell j runs from lo_j = r_{j-1} to hi_j = r_j, so the lower boundary of cell
## j + 1 is the upper boundary of cell j: Phi and phi are evaluated once on the
## N - 1 interior boundaries and shared, with the limits 0 and 1 (Phi) and 0
## (phi) at -Inf and +Inf.  The values are identical to evaluating both ends.
rmq_pieces <- function(g, c, m) {
  N  <- length(g)
  r  <- if (N > 1) 0.5 * (g[-1] + g[-N]) else numeric(0)
  z  <- outer(-c, r, "+") / m                      # [i, j] = (r_j - c_i) / m
  Fz <- pnorm(z); fz <- dnorm(z)
  list(Fa = cbind(0, Fz), Fb = cbind(Fz, 1), fa = cbind(0, fz), fb = cbind(fz, 0))
}

## Gradient of the distortion and its tridiagonal.Used by the Newton step below 
## and, at the converged grid, by the implicit-function-theorem adjoint of
## 14_quantisation_adjoint.R, so that both use exactly the same Hessian.
rmq_grad_hess <- function(g, c, m, pk, pc = rmq_pieces(g, c, m)) {
  N   <- length(g)
  dF  <- pc$Fb - pc$Fa
  grad <- 2 * colSums(pk * (outer(-c, g, "+") * dF + m * (pc$fb - pc$fa)))
  ub  <- colSums(pk * pc$fb) / m                   # mixture density at r^{j+}
  ua  <- colSums(pk * pc$fa) / m                   # mixture density at r^{j-}
  dia <- 2 * colSums(pk * dF)
  if (N > 1) {
    dg  <- diff(g)
    dia[1:(N - 1)] <- dia[1:(N - 1)] + 0.5 * ub[1:(N - 1)] * (-dg)
    dia[2:N]       <- dia[2:N]       + 0.5 * ua[2:N]       * (-dg)
    off <- 0.5 * ub[1:(N - 1)] * (-dg)
  } else off <- numeric(0)
  list(grad = grad, dia = dia, off = off)
}

## One RMQ step: quantise sum_i pk_i N(c_i, m^2) on N points.
rmq_step <- function(pk, c, m, guess, newton_max = 30L, lloyd_max = 500L,
                     tol = 1e-12) {
  g <- sort(guess)
  N <- length(g)
  method <- "newton"
  for (it in seq_len(newton_max)) {
    gh   <- rmq_grad_hess(g, c, m, pk)
    step <- tridiag_solve(gh$off, gh$dia, gh$off, gh$grad)
    if (is.null(step) || any(!is.finite(step))) { method <- "lloyd"; break }
    gn <- g - step
    if (N > 1 && any(diff(gn) <= 0)) { method <- "lloyd"; break }
    g <- gn
    if (max(abs(step)) < tol) break
  }
  if (method == "lloyd") {
    for (it in seq_len(lloyd_max)) {
      pc  <- rmq_pieces(g, c, m)
      dF  <- pc$Fb - pc$Fa
      num <- colSums(pk * (c * dF + m * (pc$fa - pc$fb)))
      den <- colSums(pk * dF)
      gn  <- ifelse(den > 1e-300, num / pmax(den, 1e-300), g)
      if (max(abs(gn - g)) < tol) { g <- gn; break }
      g <- gn
    }
  }
  pc <- rmq_pieces(g, c, m)
  list(grid = g, trans = pc$Fb - pc$Fa, method = method)
}

# -----------------------------------------------------------------------------

## build_rmq(): the whole quantisation tree for n steps of N points.
## The tree depends only on (S0, sigma, r_D, T, n, N) -- not on alpha, the
## funding rates, the side or the regime -- so a study builds it once and reuses
## it through the cache below.  Its build time is nevertheless charged to every
## solve, so that the cost measure is not flattered by the cache.

# -----------------------------------------------------------------------------
build_rmq <- function(p, n, N) {
  t0 <- proc.time()[["elapsed"]]
  dt <- p$Tmat / n
  b  <- p$rD - 0.5 * p$sigma^2
  m  <- p$sigma * sqrt(dt)
  x0 <- log(p$S0)
  z  <- std_normal_guess(N)
  grids <- vector("list", n + 1L); probs <- vector("list", n + 1L)
  trans <- vector("list", n);      method <- character(n)
  grids[[1]] <- x0; probs[[1]] <- 1
  for (k in seq_len(n)) {
    tk1   <- k * dt
    guess <- x0 + b * tk1 + p$sigma * sqrt(tk1) * z
    st    <- rmq_step(probs[[k]], grids[[k]] + b * dt, m, guess)
    grids[[k + 1L]] <- st$grid
    trans[[k]]      <- st$trans
    probs[[k + 1L]] <- as.vector(probs[[k]] %*% st$trans)
    method[k]       <- st$method
  }
  ## moment diagnostic at T: quantisation preserves the mean of a stationary
  ## grid exactly and understates the variance by the distortion
  gT <- grids[[n + 1L]]; pT <- probs[[n + 1L]]
  mq <- sum(pT * gT); vq <- sum(pT * gT^2) - mq^2
  list(grids = grids, probs = probs, trans = trans, method = method,
       n = n, N = N, dt = dt, b = b,
       mean_err = mq - (x0 + b * p$Tmat),
       var_rel_err = vq / (p$sigma^2 * p$Tmat) - 1,
       seconds = proc.time()[["elapsed"]] - t0)
}

.rmq_cache <- new.env()
get_rmq <- function(p, n, N) {
  key <- paste(n, N, p$S0, p$sigma, p$rD, p$Tmat, sep = "|")
  if (!exists(key, envir = .rmq_cache, inherits = FALSE)) {
    assign(key, build_rmq(p, n, N), envir = .rmq_cache)
  }
  get(key, envir = .rmq_cache, inherits = FALSE)
}
