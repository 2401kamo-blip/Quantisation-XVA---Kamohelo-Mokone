# -----------------------------------------------------------------------------
##  09_experiments.R

## Grid of collateralisation levels.  A fine grid is used for the two
## deterministic benchmarks, a coarse one for the deep BSDE solver, which is by
## far the most expensive of the three.

# -----------------------------------------------------------------------------
ALPHA_FINE   <- seq(0, 1, length.out = 21L)
ALPHA_COARSE <- c(0, 0.5, 1)

## Numerical settings.  Reduce dbsde_iters for a quick pass.
SETTINGS <- list(
  pde_nx     = 801L,
  pde_nt     = 1500L,
  dbsde_N    = 20L,
  dbsde_M    = 128L,
  dbsde_iters = 1500L,
  quant_n    = 50L,
  quant_N    = 100L
)

# -----------------------------------------------------------------------------

## par_lapply(): independent runs (training runs, valuations at different
## parameters) spread over cores.  Every deep training run sets its own seed,
## so the results do not depend on the number of cores or on the schedule.
## Forks on Unix; on Windows, where forking is unavailable, a PSOCK cluster is
## started and the global environment (every function of the project) is copied
## to it.  options(xva.cores = 1) runs serially.

# -----------------------------------------------------------------------------
xva_cores <- function() getOption("xva.cores", max(1L, parallel::detectCores()))

par_lapply <- function(X, FUN, ...) {
  nc <- min(xva_cores(), length(X))
  if (nc <= 1L) return(lapply(X, FUN, ...))
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(nc)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, ls(globalenv(), all.names = TRUE), envir = globalenv())
    return(parallel::parLapplyLB(cl, X, FUN, ...))
  }
  res <- parallel::mclapply(X, FUN, ..., mc.cores = nc, mc.preschedule = FALSE)
  bad <- vapply(res, inherits, NA, what = "try-error")
  if (any(bad)) stop("par_lapply(): ", res[[which(bad)[1]]])
  res
}

# -----------------------------------------------------------------------------

## Suite output.  Every suite of Chapter 5 (modules 17 onwards) runs through
## run_suite(), which
##   * seeds the run,
##   * caches it: the result is reused while the parameter object, the suite's
##     settings and the source of the solver modules (00-16) and of the suite
##     itself are unchanged, and recomputed otherwise or when force = TRUE,
##   * writes every data frame of the result as a tidy CSV and the whole result
##     as an RDS to output/data/, and
##   * records the run in output/data/manifest.csv: seed, run time, cores,
##     R version, cache key, and the figures and tables the data feed.
## No suite plots; the plotting module reads these files.

# -----------------------------------------------------------------------------
OUTPUT_DATA <- file.path("output", "data")

suite_key <- function(suite_file, ...) {
  src <- list.files("R", "^(0[0-9]|1[0-6])_.*\\.R$", full.names = TRUE)
  src <- sort(c(src, file.path("R", suite_file)))
  f <- tempfile(); on.exit(unlink(f))
  writeBin(serialize(list(..., src = unname(tools::md5sum(src))), NULL), f)
  unname(tools::md5sum(f))
}

run_suite <- function(name, suite_file, run, p0, settings, seed,
                      used_by = "", force = FALSE) {
  dir.create(OUTPUT_DATA, showWarnings = FALSE, recursive = TRUE)
  rds <- file.path(OUTPUT_DATA, paste0(name, ".rds"))
  key <- suite_key(suite_file, p0, settings, seed)
  if (!force && file.exists(rds)) {
    old <- readRDS(rds)
    if (identical(old$meta$key, key)) {
      cat(sprintf("%s: cached result (%s, %.0f s when computed)\n",
                  name, format(old$meta$created), old$meta$seconds))
      return(old)
    }
  }
  set.seed(seed)
  t0  <- proc.time()[["elapsed"]]
  res <- run()
  res$meta <- list(suite = name, key = key, seed = seed, created = Sys.time(),
                   seconds = proc.time()[["elapsed"]] - t0, cores = xva_cores(),
                   R = R.version.string, platform = R.version$platform,
                   used_by = used_by, params = p0, settings = settings)
  saveRDS(res, rds)
  for (nm in names(res)) if (is.data.frame(res[[nm]])) {
    write.csv(res[[nm]], file.path(OUTPUT_DATA, sprintf("%s_%s.csv", name, nm)),
              row.names = FALSE)
  }
  record_manifest(res$meta)
  res
}

record_manifest <- function(meta) {
  f   <- file.path(OUTPUT_DATA, "manifest.csv")
  row <- data.frame(suite = meta$suite, created = format(meta$created),
                    seconds = round(meta$seconds, 1), seed = meta$seed,
                    cores = meta$cores, R = meta$R, key = meta$key,
                    used_by = meta$used_by, stringsAsFactors = FALSE)
  if (file.exists(f)) {
    old <- read.csv(f, stringsAsFactors = FALSE)
    row <- rbind(old[old$suite != meta$suite, names(row)], row)
  }
  write.csv(row, f, row.names = FALSE)
}

# -----------------------------------------------------------------------------

## run_experiment(): sweeps alpha for one fixed configuration.

# -----------------------------------------------------------------------------
run_experiment <- function(base, alphas_fine = ALPHA_FINE,
                           alphas_coarse = ALPHA_COARSE,
                           do_dbsde = TRUE, do_quant = TRUE, op = NULL,
                           seed = 3L, settings = SETTINGS) {

  if (is.null(op)) {
    op <- pde_operator(base, settings$pde_nx, settings$pde_nt)
  }

  pde <- numeric(length(alphas_fine))
  lin <- rep(NA_real_, length(alphas_fine))
  lin_err <- rep(NA_real_, length(alphas_fine))     # closed-form branch error estimate
  qnt <- rep(NA_real_, length(alphas_fine))
  for (i in seq_along(alphas_fine)) {
    p <- base; p$alpha <- alphas_fine[i]
    pde[i] <- pde_value(p, op)
    if (!identical(p$beta_mode, "dynamic")) {
      ## the closed form is a reference only where its assumed branch holds
      ## (linear_branch_check() in 08); elsewhere it is reported as NA
      lb <- linear_benchmark(p, check = TRUE)
      lin[i]     <- if (lb$branch_ok) lb$value else NA_real_
      lin_err[i] <- lb$branch_err
    }
    if (do_quant) {                       # cheap enough for the fine grid
      qnt[i] <- quant_value(p, n = settings$quant_n, N = settings$quant_N)
    }
  }

  dbs <- rep(NA_real_, length(alphas_coarse))
  lss <- rep(NA_real_, length(alphas_coarse))
  if (do_dbsde) {
    for (i in seq_along(alphas_coarse)) {
      p <- base; p$alpha <- alphas_coarse[i]
      fit <- deep_bsde(p, N = settings$dbsde_N, M = settings$dbsde_M,
                       iters = settings$dbsde_iters, seed = seed, verbose = FALSE)
      dbs[i] <- fit$U0
      lss[i] <- mean(tail(fit$loss, 100L))
    }
  }

  list(alpha_fine = alphas_fine, pde = pde, linear = lin, linear_branch_err = lin_err,
       quant = qnt,
       alpha_coarse = alphas_coarse, dbsde = dbs, loss = lss)
}

# -----------------------------------------------------------------------------

## experiment_figure4(): the two frozen-regime panels.

# -----------------------------------------------------------------------------
experiment_figure4 <- function(p0 = make_params(), rfm = c(0.05, 0.075, 0.1),
                               do_dbsde = TRUE, settings = SETTINGS) {
  op  <- pde_operator(p0, settings$pde_nx, settings$pde_nt)
  out <- list()
  for (bm in c("frozen0", "frozen1")) {
    for (r in rfm) {
      base <- p0; base$beta_mode <- bm; base$rf_m <- r
      key  <- paste0(bm, "_rfm", r)
      cat("running ", key, "\n", sep = "")
      out[[key]] <- run_experiment(base, do_dbsde = do_dbsde, op = op,
                                   settings = settings)
    }
  }
  out
}

# -----------------------------------------------------------------------------

## experiment_figure5(): frozen calm regime against a dynamic status process.
##
## NOTE.  At the benchmark rates r_f^- = r_r^- = 0.05 and r_f^+ = r_r^+ = 0.05
## the calm and crisis drivers coincide identically, because the only effect of
## the regime is to move the financing of the stock inventory from the repo desk
## to the treasury desk at the same rate.  The dependence on lambda_U therefore
## vanishes.  The experiment is run at r_f^- = 0.10 as well, where the two
## regimes genuinely differ.

# -----------------------------------------------------------------------------
experiment_figure5 <- function(p0 = make_params(), mean_U = c(0.3, 1.0, 1.7),
                               rfm = c(0.05, 0.10), do_dbsde = TRUE,
                               settings = SETTINGS) {
  op  <- pde_operator(p0, settings$pde_nx, settings$pde_nt)
  out <- list()
  for (r in rfm) {
    base <- p0; base$beta_mode <- "frozen0"; base$rf_m <- r
    out[[paste0("frozen0_rfm", r)]] <-
      run_experiment(base, do_dbsde = do_dbsde, op = op, settings = settings)
    for (m in mean_U) {
      base <- make_params(beta_mode = "dynamic", rf_m = r, mean_U = m)
      key  <- paste0("dynamic_rfm", r, "_lamU", m)
      cat("running ", key, "\n", sep = "")
      out[[key]] <- run_experiment(base, do_dbsde = do_dbsde, op = op,
                                   settings = settings)
    }
  }
  out
}

# -----------------------------------------------------------------------------

## Plotting.  Base graphics only, so that the script runs on a bare R install.

# -----------------------------------------------------------------------------
plot_panel <- function(res_list, labels, main = "", ylab = expression(XVA^"+"),
                       cols = c("#1b3a6b", "#b03030", "#2e7d32"),
                       ltys = c(1, 2, 3)) {
  ylim <- range(unlist(lapply(res_list, function(r) c(r$pde, r$dbsde, r$quant))),
                na.rm = TRUE)
  plot(NA, xlim = c(0, 1), ylim = ylim, xlab = expression(alpha),
       ylab = ylab, main = main)
  grid(col = "grey85")
  for (i in seq_along(res_list)) {
    r <- res_list[[i]]
    lines(r$alpha_fine, r$pde, col = cols[i], lty = ltys[i], lwd = 2)
    if (any(!is.na(r$dbsde))) {
      points(r$alpha_coarse, r$dbsde, col = cols[i], pch = 1, cex = 1.1)
    }
    if (!is.null(r$quant) && any(!is.na(r$quant))) {
      points(r$alpha_fine, r$quant, col = cols[i], pch = 4, cex = 0.6)
    }
  }
  legend("topleft", legend = labels, col = cols[seq_along(res_list)],
         lty = ltys[seq_along(res_list)], bty = "n", cex = 0.8)
  legend("bottomright", legend = c("PDE (line)", "deep BSDE (circle)",
         "quantisation (cross)"), bty = "n", cex = 0.7)
}

# -----------------------------------------------------------------------------

## comparison_table(): tabulates the three solvers side by side, and against the
## values read off the published figures.

# -----------------------------------------------------------------------------
comparison_table <- function(res, published = NULL) {
  idx <- match(res$alpha_coarse, res$alpha_fine)
  df <- data.frame(
    alpha  = res$alpha_coarse,
    linear = res$linear[idx],
    pde    = res$pde[idx],
    quant  = if (is.null(res$quant)) NA_real_ else res$quant[idx],
    dbsde  = res$dbsde,
    loss   = res$loss
  )
  df$rel_err_quant <- (df$quant - df$pde) / df$pde
  df$rel_err_dbsde <- (df$dbsde - df$pde) / df$pde
  if (!is.null(published)) {
    df$published <- published
    df$ratio_pub <- published / df$pde
  }
  df
}
