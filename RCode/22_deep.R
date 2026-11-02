# ------------------------------------------------------------------------------

##  22_deep.R
##  Suite (f): the deep solvers against the quantisation solver.
##
##  SOLVERS.
##    plain   the deep BSDE solver of 06 with Equation (4.36) as published, no
##            jump term (use_jump = FALSE);
##    tied    the jump-augmented solver of 11 with the compensator imposed,
##            C = a(beta)(1 - 2 beta) nu (Section 4.6.3);
##    untied  the jump-augmented solver with the compensator learnt;
##    tied_exact  the tied solver with the exact grid compensator
##            E[Delta B | beta] / dt (p01 / dt, -p10 / dt) instead of the
##            intensity: the test of finding E5 of the audit.
##  DESIGN.
##    main      plain, tied, untied at lambda_U in {0.3, 1.0, 1.7}, four seeds;
##    K sweep   tied and untied at lambda_U = 0.3 for K in {1, 10, 1e2, 1e3, 1e4},
##              four seeds (the tied batch penalty is not identically zero);
##    E5        tied_exact at every lambda_U, four seeds.
##  A single seed is never reported alone: every configuration has four, and
##  the summary gives mean, sd and range.  Training settings are those of the
##  replication record (N = 20, 1500 iterations; M = 128 for 06, 256 for 11).
##
##  REFERENCES.  The quantisation solver extrapolated by Richardson from
##  (50, 100) and (100, 200), and the PDE on 1601 x 3000, for XVA_0 and Z_0.
##  The deep solvers discretise time with only N = 20 steps, and the explicit
##  scheme is biased at that step whatever solves it.  The same scheme at
##  n = N, solved by quantisation on a fine grid (N_grid = 400, converged to
##  1e-7), is therefore recorded too: the error against it is the LEARNING
##  error of a deep solver, and the remainder is the time-discretisation error.
##  Errors are reported in XVA_0, in Z_0 and in the stock position at t = 0,
##  xi^S_0 = (Zbreve_0 + Zhat_0) / (sigma S0), so the Z_0 error in shares is
##  (Zbreve_0 - Zbreve_0^ref) / (sigma S0).
##
##  CACHING.  Every training run is saved to output/cache/deep/ under a key of
##  its job and the source of the solver modules, so an interrupted suite
##  resumes where it stopped; run_suite() caches the whole as usual.
##
##  OUTPUT (output/data/): deep_runs.csv (one row per training run),
##  deep_summary.csv (per configuration), deep_reference.csv.
##
##  USAGE.  Rscript main.R deep [force]
# ------------------------------------------------------------------------------

DEEP_SETTINGS <- list(
  alpha   = 0, rf_m = 0.10,
  lambda_U = c(0.3, 1.0, 1.7),
  seeds   = 1:4,
  K_main  = 1e3,
  K_sweep = c(1, 10, 1e2, 1e3, 1e4), K_sweep_lambda = 0.3,
  N = 20L, iters = 1500L, M_plain = 128L, M_jump = 256L,
  ref_levels = list(c(50L, 100L), c(100L, 200L)),
  pde_ref = c(1601L, 3000L),
  disc_grid = 400L,                  # grid size of the same-step reference
  seed = 22L)

DEEP_USED_BY <- paste("Ch5 deep-versus-quantisation table (mean, sd, time; errors in XVA_0,",
  "Z_0 and shares); penalty-sweep figure; exact-compensator (E5) table")

DEEP_CACHE <- file.path("output", "cache", "deep")

deep_jobs <- function(s = DEEP_SETTINGS) {
  main <- expand.grid(solver = c("plain", "tied", "untied"), lambda_U = s$lambda_U,
                      K = s$K_main, seed = s$seeds, stringsAsFactors = FALSE)
  main$K[main$solver == "plain"] <- NA
  sweep <- expand.grid(solver = c("tied", "untied"), lambda_U = s$K_sweep_lambda,
                       K = s$K_sweep, seed = s$seeds, stringsAsFactors = FALSE)
  e5 <- expand.grid(solver = "tied_exact", lambda_U = s$lambda_U, K = s$K_main,
                    seed = s$seeds, stringsAsFactors = FALSE)
  j <- unique(rbind(main, sweep, e5))
  j$design <- ifelse(j$solver == "tied_exact", "E5",
               ifelse(!is.na(j$K) & j$K != s$K_main, "K sweep", "main"))
  j[order(j$solver == "plain"), ]          # long jobs first, for load balance
}

deep_job_key <- function(job, p0, s) {
  src <- file.path("R", c("01_status_process.R", "03_paths.R", "04_driver.R",
                          "05_network.R", "06_deep_bsde.R", "11_deep_bsde_jumps.R"))
  f <- tempfile(); on.exit(unlink(f))
  writeBin(serialize(list(job[c("solver", "lambda_U", "K", "seed")], p0,
                          s[c("alpha", "rf_m", "N", "iters", "M_plain", "M_jump")],
                          unname(tools::md5sum(src))), NULL), f)
  unname(tools::md5sum(f))
}

deep_run <- function(job, p0, s) {
  dir.create(DEEP_CACHE, showWarnings = FALSE, recursive = TRUE)
  fc <- file.path(DEEP_CACHE, paste0(deep_job_key(job, p0, s), ".rds"))
  if (file.exists(fc)) return(readRDS(fc))
  p  <- with_params(p0, beta_mode = "dynamic", alpha = s$alpha, rf_m = s$rf_m,
                    mean_U = job$lambda_U)
  t0 <- proc.time()[["elapsed"]]
  fit <- switch(job$solver,
    plain  = deep_bsde(p, N = s$N, M = s$M_plain, iters = s$iters, use_jump = FALSE,
                       seed = job$seed, verbose = FALSE),
    tied   = deep_bsde_jumps(p, N = s$N, M = s$M_jump, iters = s$iters, K = job$K,
                             tied = TRUE, seed = job$seed, verbose = FALSE),
    untied = deep_bsde_jumps(p, N = s$N, M = s$M_jump, iters = s$iters, K = job$K,
                             tied = FALSE, seed = job$seed, verbose = FALSE),
    tied_exact = deep_bsde_jumps(p, N = s$N, M = s$M_jump, iters = s$iters, K = job$K,
                                 tied = TRUE, compensator = "exact", seed = job$seed,
                                 verbose = FALSE))
  loss <- if (job$solver == "plain") fit$loss else fit$mse
  out <- data.frame(job, U0 = fit$U0, Z0 = fit$Z0[1],
                    final_mse = mean(tail(loss, 100L)),
                    final_penalty = if (job$solver == "plain") NA_real_ else mean(tail(fit$penalty, 100L)),
                    seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE)
  saveRDS(out, fc)
  out
}

deep_reference <- function(p0, s) {
  op <- pde_operator(p0, s$pde_ref[1], s$pde_ref[2])
  do.call(rbind, lapply(s$lambda_U, function(l) {
    p <- with_params(p0, beta_mode = "dynamic", alpha = s$alpha, rf_m = s$rf_m, mean_U = l)
    q <- lapply(s$ref_levels, function(lv) quantisation_bsde(p, lv[1], lv[2]))
    pv <- pde_value_control(p, op)
    qd <- quantisation_bsde(p, s$N, s$disc_grid)
    data.frame(lambda_U = l,
               disc_U0 = qd$U0, disc_Z0 = qd$Z0,
               quant_U0 = 2 * q[[2]]$U0 - q[[1]]$U0, quant_Z0 = 2 * q[[2]]$Z0 - q[[1]]$Z0,
               pde_U0 = pv[["value"]], pde_Z0 = pv[["control"]],
               Zhat0 = bs_control(0, p$S0, p), sigmaS0 = p$sigma * p$S0)
  }))
}

deep_summary <- function(runs) {
  key <- c("design", "solver", "lambda_U", "K")
  runs$K_lab <- ifelse(is.na(runs$K), "-", format(runs$K, scientific = TRUE))
  g <- split(runs, list(runs$design, runs$solver, runs$lambda_U, runs$K_lab), drop = TRUE)
  do.call(rbind, lapply(g, function(x) {
    st <- function(v) c(mean = mean(v), sd = sd(v), min = min(v), max = max(v))
    eU <- st(100 * x$err_U0_rel); eZ <- st(100 * x$err_Z0_rel); eS <- st(x$err_shares)
    lU <- st(100 * x$err_U0_rel_disc); lZ <- st(100 * x$err_Z0_rel_disc)
    data.frame(x[1, key], seeds = nrow(x),
               U0_mean = mean(x$U0), U0_sd = sd(x$U0),
               errU_pct_mean = eU[["mean"]], errU_pct_sd = eU[["sd"]],
               errU_pct_min = eU[["min"]], errU_pct_max = eU[["max"]],
               errZ_pct_mean = eZ[["mean"]], errZ_pct_sd = eZ[["sd"]],
               err_shares_mean = eS[["mean"]], err_shares_sd = eS[["sd"]],
               learnU_pct_mean = lU[["mean"]], learnU_pct_sd = lU[["sd"]],
               learnZ_pct_mean = lZ[["mean"]], learnZ_pct_sd = lZ[["sd"]],
               final_mse_mean = mean(x$final_mse), seconds_mean = mean(x$seconds),
               row.names = NULL, stringsAsFactors = FALSE)
  }))
}

# -----------------------------------------------------------------------------

## run_deep(): the suite.

# -----------------------------------------------------------------------------
run_deep <- function(p0 = make_params(), s = DEEP_SETTINGS, force = FALSE) {
  run_suite("deep", "22_deep.R", p0 = p0, settings = s, seed = s$seed,
            used_by = DEEP_USED_BY, force = force, run = function() {
    for (lv in s$ref_levels) invisible(get_nodes(p0, lv[1], lv[2]))
    ref  <- deep_reference(p0, s)
    jobs <- deep_jobs(s)
    runs <- do.call(rbind, par_lapply(seq_len(nrow(jobs)), function(i) deep_run(jobs[i, ], p0, s)))
    r <- ref[match(runs$lambda_U, ref$lambda_U), ]
    runs$ref_U0 <- r$quant_U0; runs$ref_Z0 <- r$quant_Z0
    runs$err_U0_rel <- runs$U0 / r$quant_U0 - 1
    runs$err_U0_rel_pde <- runs$U0 / r$pde_U0 - 1
    runs$err_U0_rel_disc <- runs$U0 / r$disc_U0 - 1     # learning error
    runs$err_Z0_rel_disc <- runs$Z0 / r$disc_Z0 - 1
    runs$err_Z0_rel <- runs$Z0 / r$quant_Z0 - 1
    runs$err_shares <- (runs$Z0 - r$quant_Z0) / r$sigmaS0
    runs$shares_ref <- (r$quant_Z0 + r$Zhat0) / r$sigmaS0
    list(runs = runs, summary = deep_summary(runs), reference = ref)
  })
}

## Console report used by main.R.
report_deep <- function(res) {
  cat("\nreferences (quantisation, Richardson; PDE 1601 x 3000)\n")
  print(res$reference, row.names = FALSE, digits = 6)
  sm <- res$summary
  show <- function(d, title) {
    cat("\n", title, "\n", sep = "")
    print(d[, c("solver", "lambda_U", "K", "seeds", "errU_pct_mean", "errU_pct_sd",
                "errU_pct_min", "errU_pct_max", "errZ_pct_mean", "errZ_pct_sd",
                "err_shares_mean", "err_shares_sd", "learnU_pct_mean", "learnU_pct_sd",
                "learnZ_pct_mean", "seconds_mean")], row.names = FALSE, digits = 3)
  }
  show(sm[sm$design == "main", ][order(sm$lambda_U[sm$design == "main"]), ],
       "main: errors against quantisation (XVA_0 and Z_0 in %, shares in units)")
  ks <- sm[sm$solver %in% c("tied", "untied") & sm$lambda_U == 0.3, ]
  show(ks[order(ks$solver, ks$K), ], "penalty sweep at lambda_U = 0.3")
  e5 <- sm[sm$solver %in% c("tied", "tied_exact") & sm$K %in% 1e3, ]
  show(e5[order(e5$lambda_U, e5$solver), ], "E5: intensity against exact grid compensator (tied, K = 1e3)")
  rf <- res$reference
  cat("\ntime-discretisation error of the scheme at N =", res$meta$settings$N, "steps (% against Richardson):\n")
  print(data.frame(lambda_U = rf$lambda_U, value = 100 * (rf$disc_U0 / rf$quant_U0 - 1),
                   Z0 = 100 * (rf$disc_Z0 / rf$quant_Z0 - 1)), row.names = FALSE, digits = 3)
  cat(sprintf("\n%d training runs (%.0f min single-core equivalent); suite time %.0f s on %d core(s)\n",
              nrow(res$runs), sum(res$runs$seconds) / 60, res$meta$seconds, res$meta$cores))
}
