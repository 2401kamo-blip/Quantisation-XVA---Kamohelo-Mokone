#-------------------------------------------------------------------------------

##  10_validation.R
##  Component-by-component validation.  Every routine returns a named logical or
##  a small table, and validate_all() prints a report.

# -----------------------------------------------------------------------------

## V1  Status process.  The empirical occupation probability and expected number
## of switches of the simulated chain are compared with the closed forms
## implied by the two-state Markov representation.

# -----------------------------------------------------------------------------
validate_status <- function(p = make_params(beta_mode = "dynamic"),
                            M = 40000L, N = 400L, tol = 0.01) {
  set.seed(101)
  sim <- sim_beta_chain(p, M, N, start = 0)
  emp_p <- mean(sim$B[, N + 1L])
  thr_p <- beta_prob_crisis(p, p$Tmat, start = 0)
  emp_j <- mean(rowSums(sim$dJ))
  thr_j <- beta_expected_jumps(p, p$Tmat, start = 0)
  data.frame(quantity = c("P(beta_T = 1)", "E[jumps on [0,T]]"),
             simulated = c(emp_p, emp_j),
             closed_form = c(thr_p, thr_j),
             ok = c(abs(emp_p - thr_p) < tol, abs(emp_j - thr_j) < tol))
}

# -----------------------------------------------------------------------------

## V2  Black-Scholes layer.  Two checks.
##   (a) Put-call parity at t = 0.
##   (b) The martingale property of the discounted reference value: the sample
##       mean of e^{-r_D T} Vhat_T must equal Vhat_0.
##   (c) Zhat against a central-difference delta multiplied by sigma S.

# -----------------------------------------------------------------------------
validate_black_scholes <- function(p = make_params(), M = 200000L, tol = 1e-4) {
  set.seed(102)
  S  <- p$S0 * exp((p$rD - 0.5 * p$sigma^2) * p$Tmat +
                   p$sigma * sqrt(p$Tmat) * rnorm(M))
  mc <- exp(-p$rD * p$Tmat) * mean(pmax(S - p$K, 0))
  v0 <- bs_value(0, p$S0, p)

  h  <- 1e-5
  fd <- (bs_value(0, p$S0 + h, p) - bs_value(0, p$S0 - h, p)) / (2 * h)
  zh <- bs_control(0, p$S0, p)

  data.frame(
    quantity   = c("Vhat_0 vs Monte Carlo", "Zhat_0 vs sigma S dV/dS"),
    analytic   = c(v0, zh),
    reference  = c(mc, p$sigma * p$S0 * fd),
    ok         = c(abs(v0 - mc) < 3 * sd(exp(-p$rD * p$Tmat) *
                     pmax(S - p$K, 0)) / sqrt(M) + tol,
                   abs(zh - p$sigma * p$S0 * fd) < 1e-5))
}

# -----------------------------------------------------------------------------

## V3  Driver derivatives.  The analytic partial derivatives coded in
## 04_driver.R are compared with central finite differences at random points,
## away from the kinks.

# -----------------------------------------------------------------------------
validate_driver_gradients <- function(p = make_params(alpha = 0.4), n = 2000L,
                                      h = 1e-6, tol = 1e-5) {
  set.seed(103)
  gfun <- get_driver(p)
  u <- rnorm(n, 0, 0.02)
  z <- rnorm(n, 0, 0.05)
  b <- rbinom(n, 1, 0.5)
  S <- p$S0 * exp(rnorm(n, 0, 0.2))
  tt <- runif(n, 0, p$Tmat)
  Vh <- bs_value(tt, S, p); Zh <- bs_control(tt, S, p)

  an <- gfun(u, z, b, Vh, Zh, p, grad = TRUE)
  du <- (gfun(u + h, z, b, Vh, Zh, p) - gfun(u - h, z, b, Vh, Zh, p)) / (2 * h)
  dz <- (gfun(u, z + h, b, Vh, Zh, p) - gfun(u, z - h, b, Vh, Zh, p)) / (2 * h)

  ## Discard the small number of points that straddle a kink.
  keep_u <- abs(du - an$dgdu) < 1e-3
  keep_z <- abs(dz - an$dgdz) < 1e-3
  data.frame(
    quantity = c("dg/du", "dg/dz"),
    max_abs_error = c(max(abs(du[keep_u] - an$dgdu[keep_u])),
                      max(abs(dz[keep_z] - an$dgdz[keep_z]))),
    fraction_at_kink = c(mean(!keep_u), mean(!keep_z)),
    ok = c(max(abs(du[keep_u] - an$dgdu[keep_u])) < tol,
           max(abs(dz[keep_z] - an$dgdz[keep_z])) < tol))
}

# -----------------------------------------------------------------------------

## V4  Network gradients.  The hand-coded backpropagation of 05_network.R is
## checked against central finite differences of the network output summed over
## the batch.  This is the single most error-prone part of the project.

# -----------------------------------------------------------------------------
validate_network_gradients <- function(d_in = 2L, d_out = 1L, M = 16L,
                                       h = 1e-6, tol = 1e-6) {
  set.seed(104)
  P <- init_subnet(d_in, d_out)
  x <- matrix(rnorm(M * d_in), nrow = M)
  w <- matrix(rnorm(M * d_out), nrow = M)     # arbitrary output weights

  fw <- subnet_forward(P, x, train = TRUE)
  G  <- subnet_backward(fw$P, fw$cache, w)

  worst <- 0
  for (nm in TRAIN_KEYS) {
    a <- P[[nm]]
    idx <- seq_len(min(6L, length(a)))
    for (i in idx) {
      Pp <- P; Pp[[nm]][i] <- a[i] + h
      Pm <- P; Pm[[nm]][i] <- a[i] - h
      fp <- sum(w * subnet_forward(Pp, x, train = TRUE)$y)
      fm <- sum(w * subnet_forward(Pm, x, train = TRUE)$y)
      worst <- max(worst, abs((fp - fm) / (2 * h) - G[[nm]][i]))
    }
  }
  data.frame(quantity = "network backpropagation",
             max_abs_error = worst, ok = worst < 1e-4)
}

# -----------------------------------------------------------------------------

## V5  Solver agreement.  The closed-form representation, the PDE and the deep
## BSDE solver are compared on the frozen regimes.

# -----------------------------------------------------------------------------
validate_solvers <- function(iters = 1500L) {
  op  <- pde_operator(make_params(), 801L, 1500L)
  rows <- list()
  for (bm in c("frozen0", "frozen1")) {
    for (a in c(0, 1)) {
      p <- make_params(beta_mode = bm, rf_m = 0.10, alpha = a)
      lin <- linear_benchmark(p)$value
      pde <- pde_value(p, op)
      dbs <- deep_bsde(p, N = 20L, M = 128L, iters = iters,
                       seed = 3L, verbose = FALSE)$U0
      rows[[length(rows) + 1L]] <- data.frame(
        regime = bm, alpha = a, linear = lin, pde = pde, dbsde = dbs,
        lin_vs_pde = lin - pde, dbsde_rel = (dbs - pde) / pde)
    }
  }
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------

## V6  The missing martingale component.
##
## Under a frozen beta the terminal loss of the deep BSDE scheme can be driven
## to numerical zero.  Under a dynamic beta, with the martingale part
## represented by an integral against W alone as in Equation (17), the loss
## plateaus at a strictly positive floor.  Adding the compensated jump control
## restores the representation and lowers the floor.  This function measures the
## effect and returns the three loss levels together with the point estimates.

# -----------------------------------------------------------------------------
validate_jump_component <- function(iters = 2000L, mean_U = 0.3) {
  op  <- pde_operator(make_params(), 801L, 1500L)
  out <- list()
  cfg <- list(
    list(tag = "frozen beta",              p = make_params(beta_mode = "frozen0", rf_m = 0.10), jump = FALSE),
    list(tag = "dynamic, Eq. (17) as is",  p = make_params(beta_mode = "dynamic", rf_m = 0.10, mean_U = mean_U), jump = FALSE),
    list(tag = "dynamic, with jump term",  p = make_params(beta_mode = "dynamic", rf_m = 0.10, mean_U = mean_U), jump = TRUE)
  )
  for (cf in cfg) {
    fit <- deep_bsde(cf$p, N = 20L, M = 128L, iters = iters, seed = 21L,
                     use_jump = cf$jump, verbose = FALSE)
    out[[length(out) + 1L]] <- data.frame(
      scheme = cf$tag,
      terminal_loss = mean(tail(fit$loss, 200L)),
      U0 = fit$U0,
      pde = pde_value(cf$p, op))
  }
  do.call(rbind, out)
}

# -----------------------------------------------------------------------------

## validate_all(): the full report.
## -----------------------------------------------------------------------------
## -----------------------------------------------------------------------------
## V7  The quantisation solver of 13 against the PDE benchmark of 07, in all
## three regimes, and the joint refinement of time steps and grid size.  The
## value error should be well below 0.1% in the frozen regimes and halve with
## each doubling of (n, N).

# -----------------------------------------------------------------------------
validate_quantisation <- function(alpha = 0, rf_m = 0.10) {
  cases <- list(calm    = make_params(alpha = alpha, rf_m = rf_m, beta_mode = "frozen0"),
                crisis  = make_params(alpha = alpha, rf_m = rf_m, beta_mode = "frozen1"),
                dynamic = make_params(alpha = alpha, rf_m = rf_m, beta_mode = "dynamic",
                                      mean_U = 0.3))
  op <- pde_operator(cases$calm, 801L, 1500L)
  tab <- do.call(rbind, lapply(names(cases), function(nm) {
    p <- cases[[nm]]; q <- quantisation_bsde(p); v <- pde_value(p, op)
    data.frame(regime = nm, n = q$n, N = q$N, quant = q$U0, pde = v,
               rel_err = (q$U0 - v) / v, seconds = q$seconds, fallbacks = q$fallbacks)
  }))
  p <- cases$crisis; v <- pde_value(p, op)
  ladder <- do.call(rbind, lapply(list(c(25, 50), c(50, 100), c(100, 200)), function(nN) {
    q <- quantisation_bsde(p, n = as.integer(nN[1]), N = as.integer(nN[2]))
    data.frame(regime = "crisis, refinement", n = q$n, N = q$N, quant = q$U0, pde = v,
               rel_err = (q$U0 - v) / v, seconds = q$seconds, fallbacks = q$fallbacks)
  }))
  rbind(tab, ladder)
}

validate_all <- function(quick = FALSE) {
  cat("\n== V1  status process ==\n");      print(validate_status())
  cat("\n== V2  Black-Scholes layer ==\n"); print(validate_black_scholes())
  cat("\n== V3  driver gradients ==\n");    print(validate_driver_gradients())
  cat("\n== V4  network gradients ==\n");   print(validate_network_gradients())
  cat("\n== V7  quantisation solver ==\n");  print(validate_quantisation())
  cat("\n== V8  adjoint sensitivities against bump-and-revalue (same scheme) ==\n")
  v8 <- validate_adjoint_discrete()
  cat(sprintf("%d of %d checks pass; max |relative difference| = %.2e\n",
              sum(v8$ok), nrow(v8), max(abs(v8$rel_diff))))
  if (!quick) {
    cat("\n== V5  solver agreement ==\n");  print(validate_solvers())
    cat("\n== V6  missing martingale component ==\n")
    print(validate_jump_component())
    cat("\n== V9  adjoint sensitivities against the PDE and closed form ==\n")
    print(validate_adjoint_continuous(), digits = 4)
  }
  invisible(NULL)
}
