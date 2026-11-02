#------------------------------------------------------------------------------
##  08_linear_benchmark.R

## Standard normal quadrature: nodes and weights of a truncated trapezoid rule
## on [-zmax, zmax], accurate to better than 1e-12 for zmax = 8.
gauss_normal_nodes <- function(nq = 1201L, zmax = 8) {
  z <- seq(-zmax, zmax, length.out = nq)
  w <- dnorm(z)
  h <- z[2] - z[1]
  w <- w * h
  w[1] <- w[1] / 2; w[nq] <- w[nq] / 2
  w <- w / sum(w)                                   # renormalise to unit mass
  list(z = z, w = w)
}

# -----------------------------------------------------------------------------

## linear_benchmark(): evaluation of (L2).
##   beta    the frozen regime, 0 or 1
##   branch  which funding rate the sign of A selects

# -----------------------------------------------------------------------------
linear_benchmark <- function(p, beta = NULL, branch = "auto",
                             ns = 400L, nq = 1201L, check = FALSE) {

  if (is.null(beta)) beta <- if (identical(p$beta_mode, "frozen1")) 1L else 0L

  rho_f <- switch(branch,
    rf_plus  = p$rf_p,
    rf_minus = p$rf_m,
    auto = {
      ## In the crisis regime A carries -Z/sigma, of order the stock inventory,
      ## which dominates every other term, so A < 0 and rho_f = r_f^-.  In the
      ## calm regime A = -u + (1 - L_I)(1 - alpha) Vhat, which is positive for
      ## alpha < 1 and negative at alpha = 1.
      if (beta == 1L) p$rf_m else if (p$alpha >= 1 - 1e-12) p$rf_m else p$rf_p
    },
    stop("linear_benchmark(): unknown branch '", branch, "'")
  )

  lc <- linear_coefficients(p, beta, rho_f)
  mu <- p$rD + p$sigma * lc$c2                      # Qtilde drift of the stock

  q  <- gauss_normal_nodes(nq)
  sg <- seq(0, p$Tmat, length.out = ns + 1L)
  val <- numeric(ns + 1L)

  for (k in seq_along(sg)) {
    s  <- sg[k]
    sq <- if (s > 0) sqrt(s) else 0
    S  <- p$S0 * exp((mu - 0.5 * p$sigma^2) * s + p$sigma * sq * q$z)
    Zh <- bs_control(s, S, p)
    c0 <- linear_c0(s, S, p, rho_f)
    val[k] <- exp(lc$c1 * s) * sum(q$w * (lc$c2 * Zh + c0))
  }

  ## trapezoidal outer quadrature
  h  <- sg[2] - sg[1]
  out <- h * (sum(val) - 0.5 * (val[1] + val[ns + 1L]))

  res <- list(value = out, rho_f = rho_f, c1 = lc$c1, c2 = lc$c2, drift = mu)
  if (check) res <- c(res, linear_branch_check(p, rho_f))
  res
}

# -----------------------------------------------------------------------------

## linear_branch_check(): is the closed form valid at p?
##
## The closed form fixes the funding rate rho_f (the sign of A) and the repo rate
## of the Z > 0 branch everywhere.  Section 4.7.2 requires the assumed branch to
## be checked after the fact.  Where A has the other sign, the true driver
## exceeds the assumed one by (r_f^- - r_f^+)|A|, and where Z < 0 it differs by
## (r_r^- - r_r^+)(1 - iota) Z^- / sigma.  Summing these over the nodes of the
## quantisation tree, weighted by the node probabilities and dt, estimates the
## error of the closed form.  The tree solution supplies u and z at the nodes.
## branch_ok holds when the estimate is below rtol times the value.

# -----------------------------------------------------------------------------
linear_branch_check <- function(p, rho_f, n = QUANT_SETTINGS$n,
                                N = QUANT_SETTINGS$N, rtol = 1e-5) {
  beta <- if (identical(p$beta_mode, "frozen1")) 1L else 0L
  q    <- quantisation_bsde(p, n, N, keep = TRUE)
  nds  <- get_nodes(p, n, N)
  plus <- isTRUE(all.equal(rho_f, p$rf_p)) && !isTRUE(all.equal(rho_f, p$rf_m))
  err  <- 0; wrong <- 0
  for (k in seq_len(n)) {
    nd <- nds[[k]]; u <- q$U[[k]][, 1]; Z <- q$Z[[k]][, 1] + nd$Zh
    iota <- as.numeric(Z > 0 & beta == 1)
    A  <- u - (iota / p$sigma) * Z + (theta_tilde_I(nd$Vh, p) - u) +
          (theta_tilde_C(nd$Vh, p) - u) + (1 - p$alpha) * nd$Vh
    bad <- if (plus) A < 0 else A > 0
    pk  <- q$probs[[k]]
    err <- err + q$n^-1 * p$Tmat * sum(pk * (
             (p$rf_m - p$rf_p) * abs(A) * bad +
             (p$rr_m - p$rr_p) * (1 - iota) * pmax(-Z, 0) / p$sigma))
    wrong <- wrong + sum(pk * bad) / n
  }
  list(branch_err = err, wrong_mass = wrong,
       branch_ok = err <= rtol * abs(q$U0))
}
