
##  03_paths.R
##  Forward simulation of the state process under the valuation measure Q.
##
##  The routine returns everything the backward sweep needs, so that the
##  training loop never has to recompute a Brownian increment.
## -----------------------------------------------------------------------------

simulate_paths <- function(p, M, N) {

  dt <- p$Tmat / N
  tg <- seq(0, p$Tmat, length.out = N + 1)

  ## Brownian increments under Q
  dW <- matrix(rnorm(M * N, mean = 0, sd = sqrt(dt)), nrow = M, ncol = N)

  ## exact GBM.  An explicit accumulation loop is used rather than apply() so
  ## that the shape of the result is unambiguous for every N, including N = 1.
  incr <- (p$rD - 0.5 * p$sigma^2) * dt + p$sigma * dW
  logS <- matrix(0, nrow = M, ncol = N)
  acc  <- rep(0, M)
  for (n in seq_len(N)) {
    acc <- acc + incr[, n]
    logS[, n] <- acc
  }
  S <- cbind(rep(p$S0, M), p$S0 * exp(logS))

  ## status process
  if (identical(p$beta_mode, "frozen0")) {
    B  <- matrix(0, nrow = M, ncol = N + 1)
    dJt <- matrix(0, nrow = M, ncol = N)
  } else if (identical(p$beta_mode, "frozen1")) {
    B  <- matrix(1, nrow = M, ncol = N + 1)
    dJt <- matrix(0, nrow = M, ncol = N)
  } else if (identical(p$beta_mode, "dynamic")) {
    sim <- sim_beta_chain(p, M, N, start = 0)
    B   <- sim$B
    dJt <- compensated_jumps(p, B, sim$dJ, N)
  } else {
    stop("simulate_paths(): unknown beta_mode '", p$beta_mode, "'")
  }

  ## Vhat and Zhat are needed at every node; precompute once per batch.
  Vh <- matrix(0, nrow = M, ncol = N + 1)
  Zh <- matrix(0, nrow = M, ncol = N + 1)
  for (n in seq_len(N + 1)) {
    Vh[, n] <- bs_value(tg[n], S[, n], p)
    Zh[, n] <- bs_control(tg[n], S[, n], p)
  }

  list(t = tg, dt = dt, S = S, B = B, dW = dW, dJt = dJt, Vh = Vh, Zh = Zh)
}
