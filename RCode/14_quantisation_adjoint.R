#-------------------------------------------------------------------------------

##  14_quantisation_adjoint.R
##  Adjoint algorithmic differentiation (AAD) of the quantisation solver.
##
##  Sources: Capriotti and Giles, "Algorithmic differentiation: adjoint Greeks
##  made easy" (the adjoint design paradigm); Capriotti and Lee, "Adjoint credit
##  risk management", and Henrard (2014) (the implicit function theorem for an
##  iterative calibration step); Capriotti, Jiang and Macrina (2017) (why the
##  conditional-expectation machinery must not be frozen for XVA).
##
##  WHAT IS DIFFERENTIATED.  The output is the seller's or buyer's XVA_0 as
##  computed by 13_quantisation_solver.R, with the published scheme (zmode =
##  "cgg").  The adjoint is the exact derivative of that DISCRETE quantity, so it
##  must agree with bump-and-revalue on the same scheme to the accuracy of the
##  finite differences (15_sensitivity_validation.R checks this).
##
##  THE FIVE SENSITIVITIES.
##      S0      delta   enters only through Vhat and Zhat at the nodes
##      sigma   vega    enters the grid, the transition probabilities, the
##                      control estimator (Q1), Vhat, Zhat and the driver
##      rf_m    r_f^-   the driver only
##      hI      h_I^Q   the driver only
##      hC      h_C^Q   the driver only
##
##
## ----------------------------------------------------------------------------

SENS_PARAMS <- c("S0", "sigma", "rf_m", "hI", "hC")
## The extension: the collateralisation level and the two loss rates.  Like
## r_f^-, h_I and h_C they enter the driver only, so their adjoints are
## accumulated in the same sweep; the default output stays the five above.
SENS_PARAMS_EXT <- c(SENS_PARAMS, "alpha", "LI", "LC")

quantisation_adjoint <- function(p, n = QUANT_SETTINGS$n, N = QUANT_SETTINGS$N,
                                 grid = c("full", "frozen"), params = SENS_PARAMS) {
  grid <- match.arg(grid)

# ----------------forward pass (stores U_k and Z_k) -------------------------
  rq  <- get_rmq(p, n, N)
  t0  <- proc.time()[["elapsed"]]
  fwd <- quantisation_bsde(p, n, N, zmode = "cgg", keep = TRUE)
  nds <- get_nodes(p, n, N)
  t_forward <- proc.time()[["elapsed"]] - t0

  t1   <- proc.time()[["elapsed"]]
  gadj <- get_driver_adjoint(p)
  dt <- rq$dt; b <- rq$b; sig <- p$sigma
  if (identical(p$beta_mode, "dynamic")) {
    tr   <- beta_transition(p, dt)
    Pb   <- matrix(c(1 - tr$p01, tr$p10, tr$p01, 1 - tr$p10), 2, 2)
    regs <- c(0L, 1L)
  } else {
    Pb   <- matrix(1, 1, 1)
    regs <- if (identical(p$beta_mode, "frozen1")) 1L else 0L
  }
  R <- length(regs)

  sig_bar  <- 0                                   # explicit sigma adjoint
  par_bar  <- c(rf_m = 0, hI = 0, hC = 0, alpha = 0, LI = 0, LC = 0)
  dlogS0   <- 0                                   # delta through Vhat, Zhat
  x_bar    <- lapply(rq$grids, function(g) numeric(length(g)))
  P_bar    <- vector("list", n)

  # ---- sweep 1: the backward scheme in reverse ---------------------------------
  U_bar <- matrix(0, nrow = 1L, ncol = R); U_bar[1, 1] <- 1       # d XVA_0 / d U_0
  for (k in 0:(n - 1L)) {
    P  <- rq$trans[[k + 1L]]
    nd <- nds[[k + 1L]]
    x  <- rq$grids[[k + 1L]]
    D  <- nd$D
    UPb  <- fwd$U[[k + 2L]] %*% t(Pb)
    EU   <- P %*% UPb
    EUdx <- nd$PD %*% UPb
    Z    <- fwd$Z[[k + 1L]]
    S  <- nd$S
    Vh <- nd$Vh
    Zh <- nd$Zh

    EU_bar <- U_bar
    Z_bar  <- matrix(0, nrow = nrow(U_bar), ncol = R)
    Vh_bar <- numeric(length(x)); Zh_bar <- numeric(length(x))
    for (j in seq_len(R)) {
      a <- gadj(EU[, j], Z[, j], rep(regs[j], length(x)), Vh, Zh, p,
                gbar = dt * U_bar[, j])
      EU_bar[, j] <- EU_bar[, j] + a$u
      Z_bar[, j]  <- a$z
      Vh_bar  <- Vh_bar + a$Vh
      Zh_bar  <- Zh_bar + a$Zh
      sig_bar <- sig_bar + sum(a$sigma)
      par_bar <- par_bar + c(sum(a$rf_m), sum(a$hI), sum(a$hC),
                             sum(a$alpha), sum(a$LI), sum(a$LC))
    }

    ## (Q1): Z = EUdx / (sigma dt) - (b / sigma) EU,  b = r_D - sigma^2 / 2
    EU_bar   <- EU_bar - (b / sig) * Z_bar
    EUdx_bar <- Z_bar / (sig * dt)
    sig_bar  <- sig_bar + sum(Z_bar * (-EUdx / (sig^2 * dt) + (b / sig^2 + 1) * EU))

    ## Vhat and Zhat at the nodes
    bp <- nd$bp
    sig_bar <- sig_bar + sum(Vh_bar * bp$dV_dsigma + Zh_bar * bp$dZ_dsigma)
    xS <- (Vh_bar * bp$dV_dS + Zh_bar * bp$dZ_dS) * S
    dlogS0 <- dlogS0 + sum(xS)
    x_bar[[k + 1L]] <- x_bar[[k + 1L]] + xS

    ## the two matrix products
    M1 <- EU_bar   %*% t(UPb)
    M2 <- EUdx_bar %*% t(UPb)
    P_bar[[k + 1L]] <- M1 + M2 * D
    D_bar <- M2 * P
    x_bar[[k + 1L]] <- x_bar[[k + 1L]] - rowSums(D_bar)
    x_bar[[k + 2L]] <- x_bar[[k + 2L]] + colSums(D_bar)
    U_bar <- t(P) %*% EU_bar %*% Pb + t(nd$PD) %*% EUdx_bar %*% Pb
  }
  t_sweep1 <- proc.time()[["elapsed"]] - t1

  # ---- sweep 2: the grid, by the implicit function theorem ---------------------
  t2 <- proc.time()[["elapsed"]]
  sig_grid <- 0
  dlogS0_root <- NA_real_
  if (grid == "full") {
    m <- sig * sqrt(dt)
    m_bar <- 0; b_bar <- 0
    p_bar <- lapply(rq$probs, function(v) numeric(length(v)))
    fin <- function(v) ifelse(is.finite(v), v, 0)
    for (k in (n - 1L):0L) {
      g1 <- rq$grids[[k + 2L]]; g0 <- rq$grids[[k + 1L]]
      pk <- rq$probs[[k + 1L]]; Pk <- rq$trans[[k + 1L]]
      N1 <- length(g1)
      cc <- g0 + b * dt

      ## weights: p_{k+1} = p_k P_k
      p_bar[[k + 1L]] <- p_bar[[k + 1L]] + as.vector(Pk %*% p_bar[[k + 2L]])
      Pbar_k <- P_bar[[k + 1L]] + outer(pk, p_bar[[k + 2L]])

      ## transitions: P_k = Phi(B) - Phi(A)
      r  <- 0.5 * (g1[-1] + g1[-N1])
      lo <- c(-Inf, r); hi <- c(r, Inf)
      ## The pieces of rmq_pieces() at the converged grid, computed once and
      ## shared with the Hessian below.  fA0 and fB0 are the standardised
      ## boundaries A and B with the infinite ends replaced by 0.
      pc  <- rmq_pieces(g1, cc, m)
      z   <- outer(-cc, r, "+") / m
      fA  <- pc$fa; fB <- pc$fb
      fA0 <- cbind(0, z);  fB0 <- cbind(z, 0)
      Bb <- Pbar_k * fB
      Ab <- -Pbar_k * fA
      hi_bar <- colSums(Bb) / m
      lo_bar <- colSums(Ab) / m
      c_bar  <- -(rowSums(Bb) + rowSums(Ab)) / m
      m_bar  <- m_bar - sum(Bb * fB0 + Ab * fA0) / m
      g1_bar <- x_bar[[k + 2L]]
      if (N1 > 1) {
        j <- 1:(N1 - 1)
        g1_bar[j]     <- g1_bar[j]     + 0.5 * hi_bar[j]     + 0.5 * lo_bar[j + 1]
        g1_bar[j + 1] <- g1_bar[j + 1] + 0.5 * hi_bar[j]     + 0.5 * lo_bar[j + 1]
      }

      ## the grid: grad D(g1; c, p_k, m) = 0, implicit function theorem
      gh  <- rmq_grad_hess(g1, cc, m, pk, pc)
      lam <- tridiag_solve(gh$off, gh$dia, gh$off, g1_bar)
      if (is.null(lam)) stop("quantisation_adjoint(): singular Hessian at step ", k)
      dPhi <- pc$Fb - pc$Fa
      G    <- matrix(g1, nrow = length(cc), ncol = N1, byrow = TRUE)
      Cm   <- matrix(cc, nrow = length(cc), ncol = N1)
      Hm   <- matrix(fin(hi), nrow = length(cc), ncol = N1, byrow = TRUE)
      Lm   <- matrix(fin(lo), nrow = length(cc), ncol = N1, byrow = TRUE)
      f    <- (G - Cm) * dPhi + m * (fB - fA)
      fc   <- -dPhi + ((Hm - G) * fB - (Lm - G) * fA) / m
      fm   <- -(G - Cm) * (fB0 * fB - fA0 * fA) / m + (fB - fA) +
              (fB0^2 * fB - fA0^2 * fA)
      c_bar <- c_bar - 2 * pk * as.vector(fc %*% lam)
      p_bar[[k + 1L]] <- p_bar[[k + 1L]] - 2 * as.vector(f %*% lam)
      m_bar <- m_bar - 2 * sum(pk * as.vector(fm %*% lam))

      ## centres: c = Gamma_k + b dt
      x_bar[[k + 1L]] <- x_bar[[k + 1L]] + c_bar
      b_bar <- b_bar + dt * sum(c_bar)
    }
    sig_grid    <- sqrt(dt) * m_bar - sig * b_bar
    dlogS0_root <- sum(x_bar[[1L]])              # the root is log S0
  }
  t_sweep2 <- proc.time()[["elapsed"]] - t2

  grads <- c(S0 = dlogS0 / p$S0, sigma = sig_bar + sig_grid,
             rf_m = unname(par_bar["rf_m"]), hI = unname(par_bar["hI"]),
             hC = unname(par_bar["hC"]), alpha = unname(par_bar["alpha"]),
             LI = unname(par_bar["LI"]), LC = unname(par_bar["LC"]))[params]
  list(U0 = fwd$U0, grads = grads, grid = grid,
       vega_explicit = sig_bar, vega_grid = sig_grid,
       delta_root = if (is.na(dlogS0_root)) NA_real_ else dlogS0_root / p$S0,
       seconds = c(grid_build = rq$seconds, forward = t_forward,
                   adjoint_scheme = t_sweep1, adjoint_grid = t_sweep2))
}
