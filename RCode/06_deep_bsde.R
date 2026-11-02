#-------------------------------------------------------------------------------

##  06_deep_bsde.R
##  The deep BSDE method of E, Han and Jentzen (2017), reference [40] of Pang
##  and Sturm, applied to the reduced XVA BSDE (17).

##------------------------------------------------------------------------------
deep_bsde <- function(p,
                      N        = 25L,
                      M        = 256L,
                      iters    = 4000L,
                      lr_start = 1e-2,
                      lr_end   = 1e-3,
                      use_jump = NULL,
                      seed     = NULL,
                      verbose  = TRUE,
                      report_every = 500L) {

  if (!is.null(seed)) set.seed(seed)

  dynamic <- identical(p$beta_mode, "dynamic")
  if (is.null(use_jump)) use_jump <- dynamic
  if (use_jump && !dynamic) use_jump <- FALSE          # nothing to compensate

  d_in  <- if (dynamic) 2L else 1L
  d_out <- if (use_jump) 2L else 1L
  dt    <- p$Tmat / N
  gfun  <- get_driver(p)

  ## ---- parameters ----------------------------------------------------------
  nets <- vector("list", N - 1L)
  for (k in seq_len(N - 1L)) nets[[k]] <- init_subnet(d_in, d_out)

  theta <- list(U0 = 0, Z0 = rep(0, d_out))
  st_theta <- adam_init(theta)
  st_nets  <- vector("list", N - 1L)
  for (k in seq_len(N - 1L)) st_nets[[k]] <- adam_init(nets[[k]][TRAIN_KEYS])

  loss_hist <- numeric(iters)

  for (it in seq_len(iters)) {

    lr <- lr_start * (lr_end / lr_start)^((it - 1) / max(iters - 1, 1))
    path <- simulate_paths(p, M, N)

    ## ---- forward sweep, Equation (DB1)/(DB3) -------------------------------
    U      <- rep(theta$U0, M)
    Us     <- vector("list", N + 1L); Us[[1]] <- U
    Zs     <- vector("list", N)
    caches <- vector("list", N)

    for (n in seq_len(N)) {
      if (n == 1L) {
        Zb <- matrix(rep(theta$Z0, each = M), nrow = M, ncol = d_out)
      } else {
        x  <- if (d_in == 1L) matrix(path$S[, n], ncol = 1L)
              else cbind(path$S[, n], path$B[, n])
        fw <- subnet_forward(nets[[n - 1L]], x, train = TRUE)
        nets[[n - 1L]] <- fw$P
        Zb <- fw$y
        caches[[n]] <- fw$cache
      }
      Zs[[n]] <- Zb

      gv <- gfun(U, Zb[, 1], path$B[, n], path$Vh[, n], path$Zh[, n], p)
      U  <- U - gv * dt + Zb[, 1] * path$dW[, n]
      if (d_out == 2L) U <- U + Zb[, 2] * path$dJt[, n]
      Us[[n + 1L]] <- U
    }

    loss_hist[it] <- mean(U^2)

    ## ---- backward sweep ----------------------------------------------------
    ## dL/dU_N = 2 U_N / M.  Differentiating (DB1),
    ##   dL/dU_n  = dL/dU_{n+1} ( 1 - dg/du dt )
    ##   dL/dZ_n  = dL/dU_{n+1} ( -dg/dz dt + dW_n )
    ##   dL/dZb_n = dL/dU_{n+1} dJtilde_n
    dU     <- 2 * U / M
    gtheta <- list(U0 = 0, Z0 = rep(0, d_out))
    gnets  <- vector("list", N - 1L)

    for (n in rev(seq_len(N))) {
      Zb <- Zs[[n]]
      gr <- gfun(Us[[n]], Zb[, 1], path$B[, n], path$Vh[, n], path$Zh[, n],
                 p, grad = TRUE)
      dZ <- matrix(0, nrow = M, ncol = d_out)
      dZ[, 1] <- dU * (-gr$dgdz * dt + path$dW[, n])
      if (d_out == 2L) dZ[, 2] <- dU * path$dJt[, n]
      dU <- dU * (1 - gr$dgdu * dt)

      if (n == 1L) {
        gtheta$Z0 <- colSums(dZ)
        gtheta$U0 <- sum(dU)
      } else {
        gnets[[n - 1L]] <- subnet_backward(nets[[n - 1L]], caches[[n]], dZ)
      }
    }

    ## ---- Adam updates ------------------------------------------------------
    up    <- adam_step(theta, gtheta, st_theta, lr)
    theta <- up$par; st_theta <- up$state

    for (k in seq_len(N - 1L)) {
      sub <- nets[[k]][TRAIN_KEYS]
      up  <- adam_step(sub, gnets[[k]], st_nets[[k]], lr)
      st_nets[[k]] <- up$state
      for (nm in TRAIN_KEYS) nets[[k]][[nm]] <- up$par[[nm]]
    }

    if (verbose && (it %% report_every == 0L || it == 1L)) {
      cat(sprintf("  iter %6d   lr %.2e   loss %.4e   U0 %.7f\n",
                  it, lr, loss_hist[it], theta$U0))
    }
  }

  list(U0 = theta$U0, Z0 = theta$Z0, loss = loss_hist,
       nets = nets, N = N, M = M, use_jump = use_jump)
}
