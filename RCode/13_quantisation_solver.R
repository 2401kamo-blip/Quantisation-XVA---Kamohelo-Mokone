#-------------------------------------------------------------------------------

##  13_quantisation_solver.R
##  The quantisation solver for the reduced XVA BSDE, Equation (17).
##
##  Source: Callegaro, Gnoatto and Grasselli (2023), "A fully quantization-based
##  scheme for FBSDEs", scheme (3.1), on the grids of 12_quantisation_grid.R.
##
## -----------------------------------------------------------------------------

QUANT_SETTINGS <- list(n = 50L, N = 100L, zmode = "cgg")

# -----------------------------------------------------------------------------

## Node quantities that depend only on the tree and the reference pair -- the
## prices, Vhat, Zhat and their partial derivatives at the nodes, the increments
## Delta x_k and the product Pi_k o Delta x_k -- and not on alpha, the funding
## rates, the side or the regime.  They are cached alongside the tree, so a sweep
## over alpha, r_f^-, regime or side computes them once.  The arithmetic is the
## same as computing them inline, so results are unchanged.

# -----------------------------------------------------------------------------
.node_cache <- new.env()
get_nodes <- function(p, n, N) {
  key <- paste(n, N, p$S0, p$sigma, p$rD, p$Tmat, p$K, sep = "|")
  if (!exists(key, envir = .node_cache, inherits = FALSE)) {
    rq <- get_rmq(p, n, N)
    nodes <- lapply(0:(n - 1L), function(k) {
      x  <- rq$grids[[k + 1L]]
      S  <- exp(x)
      tk <- k * rq$dt
      D  <- outer(-x, rq$grids[[k + 2L]], "+")
      list(S = S, Vh = bs_value(tk, S, p), Zh = bs_control(tk, S, p),
           bp = bs_partials(tk, S, p), D = D, PD = rq$trans[[k + 1L]] * D)
    })
    assign(key, nodes, envir = .node_cache)
  }
  get(key, envir = .node_cache, inherits = FALSE)
}

## Empties both caches, so that a timing charges the full construction.
clear_quant_caches <- function() {
  rm(list = ls(envir = .rmq_cache), envir = .rmq_cache)
  rm(list = ls(envir = .node_cache), envir = .node_cache)
}

quantisation_bsde <- function(p, n = QUANT_SETTINGS$n, N = QUANT_SETTINGS$N,
                              zmode = QUANT_SETTINGS$zmode, keep = FALSE) {
  gfun <- get_driver(p)
  rq   <- get_rmq(p, n, N)                      # cached across alpha, rates, side
  t0   <- proc.time()[["elapsed"]]
  nds  <- get_nodes(p, n, N)                    # likewise
  dt <- rq$dt; b <- rq$b

  dynamic <- identical(p$beta_mode, "dynamic")
  if (dynamic) {
    tr   <- beta_transition(p, dt)
    Pb   <- matrix(c(1 - tr$p01, tr$p10, tr$p01, 1 - tr$p10), 2, 2)   # rows = from
    regs <- c(0L, 1L)
  } else {
    Pb   <- matrix(1, 1, 1)
    regs <- if (identical(p$beta_mode, "frozen1")) 1L else 0L
  }
  R <- length(regs)

  U <- matrix(0, nrow = length(rq$grids[[n + 1L]]), ncol = R)       # Ubreve_n = 0
  U_all <- Z_all <- if (keep) vector("list", n + 1L) else NULL
  if (keep) U_all[[n + 1L]] <- U

  for (k in (n - 1L):0L) {
    P  <- rq$trans[[k + 1L]]
    nd <- nds[[k + 1L]]
    x  <- rq$grids[[k + 1L]]
    EU   <- P %*% U %*% t(Pb)
    EUdx <- nd$PD %*% U %*% t(Pb)
    if (identical(zmode, "cgg")) {
      Z <- EUdx / (p$sigma * dt) - (b / p$sigma) * EU                # (Q1)
    } else {
      Z <- (EUdx - rowSums(nd$PD) * EU) / (p$sigma * dt)
    }
    Vh <- nd$Vh
    Zh <- nd$Zh
    Un <- EU
    for (j in seq_len(R)) {                                           # (Q2)
      Un[, j] <- EU[, j] + dt * gfun(EU[, j], Z[, j], rep(regs[j], length(x)),
                                     Vh, Zh, p)
    }
    U <- Un
    if (keep) { U_all[[k + 1L]] <- U; Z_all[[k + 1L]] <- Z }
  }

  list(U0 = U[1, 1], Z0 = Z[1, 1],           # beta_0 = 0, or the frozen regime
       regimes = regs, U = U_all, Z = Z_all,
       grids = rq$grids, probs = rq$probs,
       seconds = rq$seconds + (proc.time()[["elapsed"]] - t0),
       grid_seconds = rq$seconds, fallbacks = sum(rq$method == "lloyd"),
       var_rel_err = rq$var_rel_err, n = n, N = N, zmode = zmode)
}

## Convenience wrapper, the analogue of pde_value() in 07.
quant_value <- function(p, ...) quantisation_bsde(p, ...)$U0
