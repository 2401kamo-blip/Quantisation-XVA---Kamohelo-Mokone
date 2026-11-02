#-------------------------------------------------------------------------------

##  11_deep_bsde_jumps.R
##  A second deep solver: the jump-augmented scheme of
##      Andersson, K., Gnoatto, A., Patacca, M. and Picarelli, A. (2025)
##      "A deep solver for BSDEs with jumps", SIAM J. Financial Math. 16(3),
##      875--911,
##  specialised to the reduced XVA BSDE of Pang and Sturm.
##
##------------------------------------------------------------------------------
deep_bsde_jumps <- function(p,
                            N        = 20L,
                            M        = 256L,
                            iters    = 1500L,
                            lr_start = 1e-2,
                            lr_end   = 1e-3,
                            K        = 1e3,
                            tied     = FALSE,
                            compensator = c("intensity", "exact"),
                            z_scale  = 0.1,
                            j_scale  = 0.01,
                            seed     = NULL,
                            verbose  = TRUE,
                            report_every = 500L) {

  if (!identical(p$beta_mode, "dynamic")) {
    stop("deep_bsde_jumps(): the status process must be dynamic; with a frozen ",
         "beta there is no jump component and 06_deep_bsde.R is the right solver.")
  }
  if (!is.null(seed)) set.seed(seed)

  compensator <- match.arg(compensator)
  d_in <- 2L                       # the state is the pair (S, beta)
  dt   <- p$Tmat / N
  trx  <- beta_transition(p, dt)   # exact one-step switch probabilities
  gfun <- get_driver(p)

  ## ---- parameters ----------------------------------------------------------
  Znets <- vector("list", N - 1L)  # diffusive control
  Vnets <- vector("list", N - 1L)  # jump size nu
  Cnets <- vector("list", N - 1L)  # compensator
  for (k in seq_len(N - 1L)) {
    Znets[[k]] <- init_subnet(d_in, 1L, out_scale = z_scale)
    Vnets[[k]] <- init_subnet(d_in, 1L, out_scale = j_scale)
    Cnets[[k]] <- init_subnet(d_in, 1L, out_scale = j_scale)
  }

  ## At time zero the state is deterministic, so the three controls are scalars.
  theta <- list(U0 = 0, Z0 = 0, V0 = 0, C0 = 0)
  st_theta <- adam_init(theta)
  st_Z <- st_V <- st_C <- vector("list", N - 1L)
  for (k in seq_len(N - 1L)) {
    st_Z[[k]] <- adam_init(Znets[[k]][TRAIN_KEYS])
    st_V[[k]] <- adam_init(Vnets[[k]][TRAIN_KEYS])
    st_C[[k]] <- adam_init(Cnets[[k]][TRAIN_KEYS])
  }

  mse_hist <- numeric(iters)
  pen_hist <- numeric(iters)

  for (it in seq_len(iters)) {

    lr   <- lr_start * (lr_end / lr_start)^((it - 1) / max(iters - 1, 1))
    path <- simulate_paths(p, M, N)

    ## forward jump of the status process, and its analytic intensity
    dB  <- path$B[, 2:(N + 1L), drop = FALSE] - path$B[, 1:N, drop = FALSE]
    lam <- ifelse(path$B[, 1:N, drop = FALSE] == 0, p$aU, p$aV)
    ## Tied compensator per unit of nu: the intensity a(beta)(1 - 2 beta) of
    ## (4.38)-(4.39), or the exact grid value E[Delta B | beta] / dt, which is
    ## p01 / dt from the calm regime and -p10 / dt from the crisis.  The grid
    ## chain only observes the net change of beta over a step, so the intensity
    ## form leaves a drift of order s dt / 2 relative in the jump term.
    cfac <- if (compensator == "exact") {
      ifelse(path$B[, 1:N, drop = FALSE] == 0, trx$p01 / dt, -trx$p10 / dt)
    } else lam * (1 - 2 * path$B[, 1:N, drop = FALSE])

    U   <- rep(theta$U0, M)
    Us  <- vector("list", N + 1L); Us[[1]] <- U
    store <- vector("list", N)
    pen   <- numeric(N)

    ## ---- forward sweep, Equation (J1) --------------------------------------
    for (n in seq_len(N)) {
      x <- cbind(path$S[, n], path$B[, n])
      if (n == 1L) {
        Z  <- rep(theta$Z0, M); nu <- rep(theta$V0, M); Cn <- rep(theta$C0, M)
        cz <- cv <- cc <- NULL
      } else {
        fz <- subnet_forward(Znets[[n - 1L]], x, train = TRUE)
        Znets[[n - 1L]] <- fz$P; Z <- fz$y[, 1]; cz <- fz$cache
        fv <- subnet_forward(Vnets[[n - 1L]], x, train = TRUE)
        Vnets[[n - 1L]] <- fv$P; nu <- fv$y[, 1]; cv <- fv$cache
        fc <- subnet_forward(Cnets[[n - 1L]], x, train = TRUE)
        Cnets[[n - 1L]] <- fc$P; Cn <- fc$y[, 1]; cc <- fc$cache
      }
      if (tied) Cn <- nu * cfac[, n]

      dU <- nu * dB[, n]                                   # Remark 3.1 structure
      gv <- gfun(U, Z, path$B[, n], path$Vh[, n], path$Zh[, n], p)
      U  <- U - gv * dt + Z * path$dW[, n] + dU - Cn * dt

      pen[n] <- mean(dU - Cn * dt)                         # batch martingale residual
      store[[n]] <- list(Z = Z, nu = nu, Cn = Cn, cz = cz, cv = cv, cc = cc)
      Us[[n + 1L]] <- U
    }

    mse_hist[it] <- mean(U^2)
    pen_hist[it] <- sum(pen^2)

    ## ---- backward sweep ----------------------------------------------------
    ## Differentiating (J2) through (J1),
    ##   dL/d(DeltaU_n) = dL/dU_{n+1} + 2 K pen_n / M
    ##   dL/dC_n        = -dt * ( dL/dU_{n+1} + 2 K pen_n / M )
    ##   dL/dnu_n       = dL/d(DeltaU_n) * Dbeta_n
    dU_q   <- 2 * U / M
    gtheta <- list(U0 = 0, Z0 = 0, V0 = 0, C0 = 0)
    gZ <- gV <- gC <- vector("list", N - 1L)

    for (n in rev(seq_len(N))) {
      s  <- store[[n]]
      gr <- gfun(Us[[n]], s$Z, path$B[, n], path$Vh[, n], path$Zh[, n],
                 p, grad = TRUE)
      dpen <- 2 * K * pen[n] / M
      d_dU <- dU_q + dpen
      d_C  <- -(dU_q + dpen) * dt
      d_nu <- d_dU * dB[, n]
      if (tied) {
        d_nu <- d_nu + d_C * cfac[, n]
        d_C  <- rep(0, M)
      }
      d_Z  <- dU_q * (-gr$dgdz * dt + path$dW[, n])
      dU_q <- dU_q * (1 - gr$dgdu * dt)

      if (n == 1L) {
        gtheta$Z0 <- sum(d_Z); gtheta$V0 <- sum(d_nu)
        gtheta$C0 <- sum(d_C); gtheta$U0 <- sum(dU_q)
      } else {
        gZ[[n - 1L]] <- subnet_backward(Znets[[n - 1L]], s$cz, matrix(d_Z,  ncol = 1L))
        gV[[n - 1L]] <- subnet_backward(Vnets[[n - 1L]], s$cv, matrix(d_nu, ncol = 1L))
        if (!tied) {
          gC[[n - 1L]] <- subnet_backward(Cnets[[n - 1L]], s$cc, matrix(d_C, ncol = 1L))
        }
      }
    }

    ## ---- Adam updates ------------------------------------------------------
    up <- adam_step(theta, gtheta, st_theta, lr)
    theta <- up$par; st_theta <- up$state

    for (k in seq_len(N - 1L)) {
      for (which in 1:3) {
        if (which == 3L && tied) next
        nets  <- switch(which, Znets, Vnets, Cnets)
        state <- switch(which, st_Z[[k]], st_V[[k]], st_C[[k]])
        grad  <- switch(which, gZ[[k]], gV[[k]], gC[[k]])
        sub <- nets[[k]][TRAIN_KEYS]
        up  <- adam_step(sub, grad, state, lr)
        for (nm in TRAIN_KEYS) nets[[k]][[nm]] <- up$par[[nm]]
        if (which == 1L) { Znets[[k]] <- nets[[k]]; st_Z[[k]] <- up$state }
        if (which == 2L) { Vnets[[k]] <- nets[[k]]; st_V[[k]] <- up$state }
        if (which == 3L) { Cnets[[k]] <- nets[[k]]; st_C[[k]] <- up$state }
      }
    }

    if (verbose && (it %% report_every == 0L || it == 1L)) {
      cat(sprintf("  iter %6d   lr %.2e   mse %.4e   pen %.4e   U0 %.7f\n",
                  it, lr, mse_hist[it], pen_hist[it], theta$U0))
    }
  }

  list(U0 = theta$U0, Z0 = theta$Z0, mse = mse_hist, penalty = pen_hist,
       N = N, M = M, K = K, tied = tied, compensator = compensator)
}

# -----------------------------------------------------------------------------

## compare_deep_solvers()
##
## Runs both deep solvers and the partial differential equation benchmark over a
## grid of mean calm-regime lengths, and returns the signed relative errors.
## The expected pattern is that the solver of 06 is biased downwards by an
## amount that grows with the switching intensity, while the solver of this file
## scatters around zero.

# -----------------------------------------------------------------------------
compare_deep_solvers <- function(mean_U = c(0.3, 1.0, 1.7),
                                 alphas = c(0, 0.5, 1),
                                 rf_m   = 0.10,
                                 iters  = 1500L, seed = 7L) {
  op   <- pde_operator(make_params(), 801L, 1500L)
  grid <- expand.grid(alpha = alphas, mean_U = mean_U)
  ## the cases are independent, and each training run sets its own seed
  rows <- par_lapply(seq_len(nrow(grid)), function(i) {
    m   <- grid$mean_U[i]; a <- grid$alpha[i]
    p   <- make_params(beta_mode = "dynamic", rf_m = rf_m, mean_U = m, alpha = a)
    ref <- pde_value(p, op)
    ehj <- deep_bsde(p, N = 20L, M = 128L, iters = iters, seed = seed,
                     use_jump = FALSE, verbose = FALSE)$U0
    agp <- deep_bsde_jumps(p, N = 20L, M = 256L, iters = iters, K = 1e3,
                           seed = seed, verbose = FALSE)$U0
    data.frame(
      lambda_U = m, alpha = a, pde = ref,
      han_jentzen_e = ehj, andersson = agp,
      err_hje = (ehj - ref) / ref, err_agpp = (agp - ref) / ref)
  })
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------

## compare_jump_schemes(): the four treatments of the jump term, several seeds
## each, because a single training run is not informative.
##   none    Equation (17) exactly as published (06, use_jump = FALSE)
##   single  the ad hoc repair of 06: one extra output multiplying the
##           compensated increment (06, use_jump = TRUE)
##   tied    this solver with the known compensator imposed
##   untied  this solver with the compensator learnt
## Errors are relative to the PDE benchmark of 07.  With the defaults this is
## 16 trainings; expect roughly half an hour on one core.

# -----------------------------------------------------------------------------
compare_jump_schemes <- function(mean_U = 0.3, alpha = 0, rf_m = 0.10,
                                 seeds = 1:4, iters = 1500L) {
  p   <- make_params(beta_mode = "dynamic", rf_m = rf_m, mean_U = mean_U, alpha = alpha)
  ref <- pde_value(p)
  runs <- list(
    none   = function(s) deep_bsde(p, iters = iters, use_jump = FALSE, seed = s,
                                   verbose = FALSE)$U0,
    single = function(s) deep_bsde(p, iters = iters, use_jump = TRUE, seed = s,
                                   verbose = FALSE)$U0,
    tied   = function(s) deep_bsde_jumps(p, iters = iters, tied = TRUE, seed = s,
                                         verbose = FALSE)$U0,
    untied = function(s) deep_bsde_jumps(p, iters = iters, tied = FALSE, seed = s,
                                         verbose = FALSE)$U0)
  ## every (treatment, seed) training run is independent
  jobs <- expand.grid(seed = seeds, treatment = names(runs), stringsAsFactors = FALSE)
  U0s  <- unlist(par_lapply(seq_len(nrow(jobs)), function(i)
            runs[[jobs$treatment[i]]](jobs$seed[i])))
  tab <- do.call(rbind, lapply(names(runs), function(nm) {
    e <- (U0s[jobs$treatment == nm] - ref) / ref
    data.frame(treatment = nm, mean_err = mean(e), sd_err = sd(e),
               min_err = min(e), max_err = max(e))
  }))
  attr(tab, "benchmark") <- ref
  tab
}
