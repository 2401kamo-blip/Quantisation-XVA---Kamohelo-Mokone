# ------------------------------------------------------------------------------

##  24_pang_sturm.R
##  Suite (h): the replication of Pang and Sturm (2020), Figures 4 and 5, with
##  this project's solvers, for the appendix.
##
##  FIGURE 4.  XVA^+ against alpha in the frozen calm and frozen crisis regimes
##  for r_f^- in {0.05, 0.075, 0.10}: the closed form (reported only where its
##  branch holds, linear_branch_check() in 08), the PDE (801 x 1500), the
##  quantisation solver (Richardson from (50,100) and (100,200)) on the fine
##  alpha grid, and the deep BSDE solver of the article (06, Equation (4.36) as
##  published) at alpha in {0, 0.5, 1} with four seeds.
##
##  FIGURE 5.  The frozen calm regime against the dynamic status process for
##  lambda_U in {0.3, 1.0, 1.7} at r_f^- = 0.05 (the benchmark, where the
##  curves coincide because the calm and crisis drivers are identical) and
##  r_f^- = 0.10 (where they differ).  The same deterministic solvers, the deep
##  solver of the article at both rates, and the jump-augmented tied solver
##  with K = 10 (the best treatment of Module 22) at r_f^- = 0.10.  The
##  quantisation solver is also run under the literal reading of Definition 1,
##  lambda as a RATE, to show that the two readings order the curves oppositely
##  (record, item 6).
##
##  PUBLISHED VALUES.  The article prints no numbers; the values below were read
##  off its Figures 4(a) and 4(b) by eye (replication record V4, Table 6) and
##  carry that uncertainty.  Further readings can be appended to
##  PS_PUBLISHED with their source.
##
##  OUTPUT (output/data/): ps_figure4.csv, ps_figure5.csv, ps_deep_runs.csv,
##  ps_deep_points.csv, ps_ratios.csv.
##
##  USAGE.  Rscript main.R pangsturm [force]
# ------------------------------------------------------------------------------

PS_PUBLISHED <- data.frame(
  figure = c("4(a)", "4(a)", "4(b)", "4(b)", "4(b)"),
  regime = c("frozen0", "frozen0", "frozen1", "frozen1", "frozen1"),
  rf_m   = c(0.05, 0.05, 0.05, 0.10, 0.10),
  alpha  = c(0, 1, 0, 0, 1),
  published = c(0.0075, 0.0105, 0.0098, 0.0240, 0.0280),
  source = "read off the published figure; V4 Table 6",
  stringsAsFactors = FALSE)

PS_SETTINGS <- list(
  alphas  = seq(0, 1, by = 0.05),
  deep_alphas = c(0, 0.5, 1),
  fig4_rf_m = c(0.05, 0.075, 0.10),
  fig5_rf_m = c(0.05, 0.10),
  fig5_lambda_U = c(0.3, 1.0, 1.7),
  levels  = list(c(50L, 100L), c(100L, 200L)),
  pde     = c(801L, 1500L),
  seeds   = 1:4,
  deep_N = 20L, deep_iters = 1500L, deep_M_plain = 128L, deep_M_jump = 256L,
  tied_K = 10, tied_rf_m = 0.10,
  seed = 24L)

PS_USED_BY <- paste("Appendix: reproduction of Pang and Sturm Figures 4 and 5 with all solvers;",
  "published-versus-reproduced ratio table; lambda-convention figure")

PS_CACHE <- file.path("output", "cache", "pangsturm")

ps_quant <- function(p, levels) {
  q <- vapply(levels, function(lv) quantisation_bsde(p, lv[1], lv[2])$U0, 0)
  2 * q[2] - q[1]
}

## One curve of a figure: every deterministic solver over the alpha grid.
ps_curve <- function(cs, p0, s, op) {
  p <- with_params(p0, beta_mode = cs$regime, rf_m = cs$rf_m)
  if (!is.na(cs$lambda_U)) p <- with_params(p, mean_U = cs$lambda_U)
  if (identical(cs$convention, "rate")) p <- with_params(p, lambda_convention = "rate")
  Vhat0 <- bs_value(0, p$S0, p)
  do.call(rbind, lapply(s$alphas, function(a) {
    pa <- with_params(p, alpha = a)
    cf <- if (cs$regime != "dynamic") linear_benchmark(pa, check = TRUE) else NULL
    q  <- ps_quant(pa, s$levels)
    data.frame(figure = cs$figure, regime = cs$regime, rf_m = cs$rf_m,
               lambda_U = cs$lambda_U, convention = cs$convention, alpha = a,
               closed_form = if (is.null(cf) || !cf$branch_ok) NA_real_ else cf$value,
               pde = if (cs$convention == "mean") pde_value(pa, op) else NA_real_,
               quantisation = q, quantisation_pct_Vhat0 = 100 * q / Vhat0,
               stringsAsFactors = FALSE)
  }))
}

## Deep training runs, cached one by one as in Module 22.
ps_deep_jobs <- function(s) {
  f4 <- expand.grid(figure = "4", solver = "plain", regime = c("frozen0", "frozen1"),
                    rf_m = s$fig4_rf_m, lambda_U = NA, alpha = s$deep_alphas,
                    seed = s$seeds, stringsAsFactors = FALSE)
  f5 <- expand.grid(figure = "5", solver = "plain", regime = "dynamic", rf_m = s$fig5_rf_m,
                    lambda_U = s$fig5_lambda_U, alpha = s$deep_alphas, seed = s$seeds,
                    stringsAsFactors = FALSE)
  f5t <- expand.grid(figure = "5", solver = "tied", regime = "dynamic", rf_m = s$tied_rf_m,
                     lambda_U = s$fig5_lambda_U, alpha = s$deep_alphas, seed = s$seeds,
                     stringsAsFactors = FALSE)
  j <- rbind(f4, f5, f5t)
  j[order(j$solver != "tied"), ]                     # long jobs first
}

ps_deep_run <- function(job, p0, s) {
  dir.create(PS_CACHE, showWarnings = FALSE, recursive = TRUE)
  src <- file.path("R", c("01_status_process.R", "03_paths.R", "04_driver.R",
                          "05_network.R", "06_deep_bsde.R", "11_deep_bsde_jumps.R"))
  f <- tempfile(); writeBin(serialize(list(job, p0, s[grep("^deep_|^tied_K", names(s))],
                                           unname(tools::md5sum(src))), NULL), f)
  fc <- file.path(PS_CACHE, paste0(unname(tools::md5sum(f)), ".rds")); unlink(f)
  if (file.exists(fc)) return(readRDS(fc))
  p <- with_params(p0, beta_mode = job$regime, rf_m = job$rf_m, alpha = job$alpha)
  if (!is.na(job$lambda_U)) p <- with_params(p, mean_U = job$lambda_U)
  t0 <- proc.time()[["elapsed"]]
  fit <- if (job$solver == "plain") {
    deep_bsde(p, N = s$deep_N, M = s$deep_M_plain, iters = s$deep_iters, use_jump = FALSE,
              seed = job$seed, verbose = FALSE)
  } else {
    deep_bsde_jumps(p, N = s$deep_N, M = s$deep_M_jump, iters = s$deep_iters, K = s$tied_K,
                    tied = TRUE, seed = job$seed, verbose = FALSE)
  }
  out <- data.frame(job, U0 = fit$U0, seconds = proc.time()[["elapsed"]] - t0,
                    stringsAsFactors = FALSE)
  saveRDS(out, fc)
  out
}

# -----------------------------------------------------------------------------

## run_pang_sturm(): the suite.

# -----------------------------------------------------------------------------
run_pang_sturm <- function(p0 = make_params(), s = PS_SETTINGS, force = FALSE) {
  run_suite("ps", "24_pang_sturm.R", p0 = p0, settings = s, seed = s$seed,
            used_by = PS_USED_BY, force = force, run = function() {
    for (lv in s$levels) invisible(get_nodes(p0, lv[1], lv[2]))
    op <- pde_operator(p0, s$pde[1], s$pde[2])
    c4 <- expand.grid(figure = "4", regime = c("frozen0", "frozen1"), rf_m = s$fig4_rf_m,
                      lambda_U = NA, convention = "mean", stringsAsFactors = FALSE)
    c5 <- rbind(
      expand.grid(figure = "5", regime = "frozen0", rf_m = s$fig5_rf_m, lambda_U = NA,
                  convention = "mean", stringsAsFactors = FALSE),
      expand.grid(figure = "5", regime = "dynamic", rf_m = s$fig5_rf_m,
                  lambda_U = s$fig5_lambda_U, convention = c("mean", "rate"),
                  stringsAsFactors = FALSE))
    cs <- rbind(c4, c5)
    curves <- do.call(rbind, par_lapply(seq_len(nrow(cs)), function(i) ps_curve(cs[i, ], p0, s, op)))

    jobs <- ps_deep_jobs(s)
    runs <- do.call(rbind, par_lapply(seq_len(nrow(jobs)), function(i) ps_deep_run(jobs[i, ], p0, s)))
    key  <- c("figure", "solver", "regime", "rf_m", "lambda_U", "alpha")
    runs$lambda_lab <- ifelse(is.na(runs$lambda_U), "-", runs$lambda_U)
    pts <- do.call(rbind, lapply(split(runs, runs[, c("figure", "solver", "regime", "rf_m",
                                                      "lambda_lab", "alpha")], drop = TRUE),
      function(x) data.frame(x[1, key], seeds = nrow(x), U0_mean = mean(x$U0), U0_sd = sd(x$U0),
                             U0_min = min(x$U0), U0_max = max(x$U0),
                             seconds_mean = mean(x$seconds), row.names = NULL)))
    m <- curves[curves$convention == "mean", ]
    ref <- function(r) {
      i <- match(paste(r$regime, r$rf_m, r$alpha, r$lambda_U),
                 paste(m$regime, m$rf_m, m$alpha, m$lambda_U))
      m$quantisation[i]
    }
    pts$quantisation <- ref(pts)
    pts$deep_vs_quant <- pts$U0_mean / pts$quantisation - 1

    pub <- PS_PUBLISHED
    m4 <- m[m$figure == "4", ]
    i  <- match(paste(pub$regime, pub$rf_m, pub$alpha), paste(m4$regime, m4$rf_m, m4$alpha))
    pub$quantisation <- m4$quantisation[i]; pub$pde <- m4$pde[i]; pub$closed_form <- m4$closed_form[i]
    pub$ratio_published_to_reproduced <- pub$published / pub$quantisation

    list(figure4 = curves[curves$figure == "4", ], figure5 = curves[curves$figure == "5", ],
         deep_runs = runs, deep_points = pts, ratios = pub)
  })
}

## Console report used by main.R.
report_pang_sturm <- function(res) {
  f4 <- res$figure4
  cat("\nFigure 4: XVA^+ at alpha = 0, 0.5, 1 (quantisation; PDE and closed form agree)\n")
  h <- f4[f4$alpha %in% c(0, 0.5, 1), ]
  print(reshape(h[, c("regime", "rf_m", "alpha", "quantisation")], idvar = c("regime", "rf_m"),
                timevar = "alpha", direction = "wide"), row.names = FALSE, digits = 5)
  cat(sprintf("max |quantisation / PDE - 1| = %.1e; max |closed form / PDE - 1| = %.1e\n",
              max(abs(f4$quantisation / f4$pde - 1)),
              max(abs(f4$closed_form / f4$pde - 1), na.rm = TRUE)))
  f5 <- res$figure5
  cat("\nFigure 5: XVA^+ at alpha = 0, 0.5, 1 under the two readings of lambda\n")
  h <- f5[f5$alpha %in% c(0, 0.5, 1), ]
  h$curve <- ifelse(h$regime == "frozen0", "beta = 0",
                    sprintf("lambda_U = %.1f (%s)", h$lambda_U, h$convention))
  print(reshape(h[, c("rf_m", "curve", "alpha", "quantisation")], idvar = c("rf_m", "curve"),
                timevar = "alpha", direction = "wide"), row.names = FALSE, digits = 5)
  cat("\ndeep solvers against quantisation (mean of four seeds)\n")
  print(res$deep_points[, c("figure", "solver", "regime", "rf_m", "lambda_U", "alpha",
                            "U0_mean", "U0_sd", "quantisation", "deep_vs_quant")],
        row.names = FALSE, digits = 4)
  cat("\npublished (read off the figures) against reproduced\n")
  print(res$ratios[, c("figure", "rf_m", "alpha", "published", "quantisation", "pde",
                       "ratio_published_to_reproduced")], row.names = FALSE, digits = 4)
  cat(sprintf("\n%d deep runs (%.0f min single-core); suite time %.0f s on %d core(s)\n",
              nrow(res$deep_runs), sum(res$deep_runs$seconds) / 60, res$meta$seconds,
              res$meta$cores))
}
