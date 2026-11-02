#-------------------------------------------------------------------------------

##  16_regression.R
##  Regression tests against the published record of the replication,
##  XVA_V4.pdf (the "V4" document, made to record methodological progress through study).
##
##  PURPOSE.  Every later change to the code (position recovery, experiment
##  suites, plotting) must leave the numbers of V4 unchanged.  This file turns
##  the DETERMINISTIC tables of V4 into automated checks.  The stochastic tables
##  (deep solvers: V4 Tables 2, 5, 7, 8, 9) are not regression targets, because a
##  single training run is not reproducible across platforms to the printed
##  digits; they are re-run by the experiment suites instead.
##
##  TOLERANCE.  V4 prints its numbers rounded.  A check passes when the computed
##  value agrees with the printed one to within half a unit in the last printed
##  digit (plus a small slack of 1 per cent of that half unit for rounding of
##  the rounding).  A check therefore fails only if the printed digits cannot be
##  reproduced, which is the right definition of a regression.
##
##  USAGE.  Rscript main.R regress      (about two minutes on one core)

# ------------------------------------------------------------------------------

## Half a unit in the last printed digit of a number written as a string.
half_unit <- function(s) {
  s <- sub("^[-+]", "", s)
  if (grepl("e", s, ignore.case = TRUE)) {            # scientific: mantissa digits
    mant <- sub("[eE].*$", "", s)
    ex   <- as.numeric(sub("^.*[eE]", "", s))
    dec  <- if (grepl("\\.", mant)) nchar(sub("^.*\\.", "", mant)) else 0
    return(0.5 * 10^(ex - dec))
  }
  dec <- if (grepl("\\.", s)) nchar(sub("^.*\\.", "", s)) else 0
  0.5 * 10^(-dec)
}

## One regression record: computed value against a printed target.
## 'known' marks a printed value that V4 itself gets wrong (an erratum in the
## record, not a regression of the code); such rows are reported separately.
reg_row <- function(table, item, computed, printed, known = NA_character_) {
  target <- as.numeric(printed)
  tol    <- 1.01 * half_unit(printed)
  data.frame(table = table, item = item, printed = printed,
             computed = signif(computed, 8), diff = signif(computed - target, 3),
             tol = signif(tol, 2), ok = abs(computed - target) <= tol,
             known = known, value = computed, stringsAsFactors = FALSE)
}

## V4 Table 4 prints a closed-form column that disagrees with V4 Table 10 for
## the same cases (for example 0.0101019 against 0.010108 in the crisis at
## alpha = 0).  The current closed form reproduces Table 10 and agrees with the
## PDE to about 5e-6 in relative terms; Table 4 predates the final quadrature.
V4_TABLE4_NOTE <- "V4 erratum: Table 4 closed-form column is stale (see Table 10)"

## PDE value and control at (0, S0) in the starting regime.
pde_value_control <- function(p, op) {
  res <- pde_solve(p, op)
  key <- if (identical(p$beta_mode, "frozen1")) "1" else "0"
  u   <- res$U[[key]]
  x0  <- log(p$S0)
  ux  <- d_dx(u, op$dx)
  c(value = approx(op$x, u, x0)$y, control = p$sigma * approx(op$x, ux, x0)$y)
}

# -----------------------------------------------------------------------------

## The targets, transcribed from V4.

# -----------------------------------------------------------------------------
V4_TARGETS <- list(

  ## Section 4: reference values at inception
  reference = c(Vhat0 = "0.060967", Zhat0 = "0.160956", delta = "0.53652"),

  ## Table 10: frozen regimes, r_f^- = 0.10, n = 50, N = 100
  table10 = data.frame(
    regime = c("frozen0", "frozen0", "frozen0", "frozen1", "frozen1", "frozen1"),
    alpha  = c(0, 0.5, 1, 0, 0.5, 1),
    closed = c("0.003477", "0.004395", "0.005348", "0.010108", "0.011287", "0.012466"),
    pde    = c("0.003477", "0.004396", "0.005348", "0.010108", "0.011287", "0.012466"),
    quant  = c("0.003479", "0.004397", "0.005347", "0.010103", "0.011279", "0.012455"),
    stringsAsFactors = FALSE),

  ## Table 4: the three solvers at r_f^- = 0.10 (closed form and PDE columns)
  table4 = data.frame(
    regime = c("frozen0", "frozen0", "frozen1", "frozen1"),
    alpha  = c(0, 1, 0, 1),
    closed = c("0.0034762", "0.0053474", "0.0101019", "0.0124606"),
    pde    = c("0.0034766", "0.0053477", "0.0101076", "0.0124661"),
    stringsAsFactors = FALSE),

  ## Table 6: reproduced values (the "This work" column)
  table6 = data.frame(
    regime = c("frozen0", "frozen0", "frozen1", "frozen1", "frozen1"),
    rf_m   = c(0.05, 0.05, 0.05, 0.10, 0.10),
    alpha  = c(0, 1, 0, 0, 1),
    value  = c("0.003477", "0.005314", "0.003477", "0.010108", "0.012466"),
    stringsAsFactors = FALSE),

  ## Table 11 (and the PDE column of Table 2): dynamic regime, r_f^- = 0.10
  table11 = data.frame(
    mean_U = rep(c(0.3, 1.0, 1.7), each = 3),
    alpha  = rep(c(0, 0.5, 1), 3),
    pde    = c("0.005421", "0.006412", "0.007436", "0.004167", "0.005111",
               "0.006089", "0.003895", "0.004830", "0.005798"),
    quant  = c("0.005387", "0.006375", "0.007397", "0.004155", "0.005098",
               "0.006074", "0.003889", "0.004822", "0.005788"),
    stringsAsFactors = FALSE),

  ## Table 12: joint refinement, crisis, alpha = 0.5, r_f^- = 0.10, relative
  ## errors in per cent against the PDE on a 1601 x 3000 grid
  table12 = data.frame(
    n = c(25L, 50L, 100L), N = c(50L, 100L, 200L),
    err_U = c("-0.138", "-0.070", "-0.035"),
    err_Z = c("-4.71", "-2.39", "-1.21"),
    stringsAsFactors = FALSE),

  ## Table 16: adjoint sensitivities, alpha = 0.5, r_f^- = 0.10, n = 50, N = 100
  table16 = list(
    "seller frozen0" = c(XVA0 = "0.004397", S0 = "0.02289", sigma = "-0.002538",
                         rf_m = "0.000000", hI = "-0.004361", hC = "-0.000537"),
    "seller frozen1" = c(XVA0 = "0.011279", S0 = "0.06022", sigma = "-0.003902",
                         rf_m = "0.1430", hI = "-0.005449", hC = "-0.001393"),
    "seller dynamic" = c(XVA0 = "0.004912", S0 = "0.02571", sigma = "-0.002624",
                         rf_m = "0.01056", hI = "-0.004456", hC = "-0.000620"),
    "buyer frozen0"  = c(XVA0 = "0.004424", S0 = "0.02303", sigma = "-0.002555",
                         rf_m = "-0.003307", hI = "-0.000541", hC = "-0.004389"),
    "buyer frozen1"  = c(XVA0 = "0.004424", S0 = "0.02303", sigma = "-0.002555",
                         rf_m = "-0.003307", hI = "-0.000541", hC = "-0.004389"),
    "buyer dynamic"  = c(XVA0 = "0.004424", S0 = "0.02303", sigma = "-0.002555",
                         rf_m = "-0.003307", hI = "-0.000541", hC = "-0.004389")),

  ## Table 17: frozen-grid vega
  table17 = c("seller frozen0" = "-0.001906", "seller frozen1" = "-0.005603",
              "seller dynamic" = "-0.002145", "buyer frozen0" = "-0.001916")
)

# -----------------------------------------------------------------------------

## run_regression(): computes every target and returns one row per check.

# -----------------------------------------------------------------------------
run_regression <- function(verbose = TRUE) {
  T0 <- proc.time()[["elapsed"]]
  out <- list()
  add <- function(r) out[[length(out) + 1L]] <<- r
  say <- function(...) if (verbose) cat(...)

  ## The V4 record was produced with the generator as printed ("extended" crisis
  ## convention).  The base parameter set now uses the admissible projection,
  ## which changes only the buyer in the crisis; the record is checked under the
  ## convention it was made with.
  make_params <- function(...) base::get("make_params", envir = globalenv())(crisis_short = "extended", ...)

  p0 <- make_params()
  op <- pde_operator(p0, SETTINGS$pde_nx, SETTINGS$pde_nt)

  ## reference values
  say("reference values\n")
  rt <- V4_TARGETS$reference
  add(reg_row("Sec. 4", "Vhat_0", bs_value(0, p0$S0, p0), rt[["Vhat0"]]))
  add(reg_row("Sec. 4", "Zhat_0", bs_control(0, p0$S0, p0), rt[["Zhat0"]]))
  add(reg_row("Sec. 4", "delta_0", bs_delta(0, p0$S0, p0), rt[["delta"]]))

  ## Table 10
  say("Table 10\n")
  tb <- V4_TARGETS$table10
  for (i in seq_len(nrow(tb))) {
    p <- make_params(beta_mode = tb$regime[i], rf_m = 0.10, alpha = tb$alpha[i])
    lab <- sprintf("%s, alpha = %.1f", tb$regime[i], tb$alpha[i])
    add(reg_row("Table 10", paste(lab, "closed form"), linear_benchmark(p)$value, tb$closed[i]))
    add(reg_row("Table 10", paste(lab, "PDE"),         pde_value(p, op),          tb$pde[i]))
    add(reg_row("Table 10", paste(lab, "quantisation"), quant_value(p),           tb$quant[i]))
  }

  ## Table 4
  say("Table 4\n")
  tb <- V4_TARGETS$table4
  for (i in seq_len(nrow(tb))) {
    p <- make_params(beta_mode = tb$regime[i], rf_m = 0.10, alpha = tb$alpha[i])
    lab <- sprintf("%s, alpha = %.0f", tb$regime[i], tb$alpha[i])
    add(reg_row("Table 4", paste(lab, "closed form"), linear_benchmark(p)$value, tb$closed[i],
                known = V4_TABLE4_NOTE))
    add(reg_row("Table 4", paste(lab, "PDE"),         pde_value(p, op),          tb$pde[i]))
  }

  ## Table 6
  say("Table 6\n")
  tb <- V4_TARGETS$table6
  for (i in seq_len(nrow(tb))) {
    p <- make_params(beta_mode = tb$regime[i], rf_m = tb$rf_m[i], alpha = tb$alpha[i])
    lab <- sprintf("%s, rf- = %.2f, alpha = %.0f, PDE", tb$regime[i], tb$rf_m[i], tb$alpha[i])
    add(reg_row("Table 6", lab, pde_value(p, op), tb$value[i]))
  }

  ## Table 11
  say("Table 11\n")
  tb <- V4_TARGETS$table11
  for (i in seq_len(nrow(tb))) {
    p <- make_params(beta_mode = "dynamic", rf_m = 0.10, mean_U = tb$mean_U[i],
                     alpha = tb$alpha[i])
    lab <- sprintf("lambda_U = %.1f, alpha = %.1f", tb$mean_U[i], tb$alpha[i])
    add(reg_row("Table 11", paste(lab, "PDE"),          pde_value(p, op), tb$pde[i]))
    add(reg_row("Table 11", paste(lab, "quantisation"), quant_value(p),   tb$quant[i]))
  }

  ## Table 12 (fine PDE reference, then the three refinement levels)
  say("Table 12 (fine PDE reference, about 30 s)\n")
  tb  <- V4_TARGETS$table12
  p   <- make_params(beta_mode = "frozen1", rf_m = 0.10, alpha = 0.5)
  opf <- pde_operator(p, 1601L, 3000L)
  ref <- pde_value_control(p, opf)
  for (i in seq_len(nrow(tb))) {
    q  <- quantisation_bsde(p, n = tb$n[i], N = tb$N[i])
    lab <- sprintf("n = %d, N = %d", tb$n[i], tb$N[i])
    add(reg_row("Table 12", paste(lab, "error in U0 (%)"),
                100 * (q$U0 / ref[["value"]] - 1), tb$err_U[i]))
    add(reg_row("Table 12", paste(lab, "error in Z0 (%)"),
                100 * (q$Z0 / ref[["control"]] - 1), tb$err_Z[i]))
  }

  ## Tables 16 and 17
  say("Tables 16 and 17 (adjoint)\n")
  scen <- lapply(sensitivity_scenarios(), function(p) { p$crisis_short <- "extended"; p })
  for (nm in names(V4_TARGETS$table16)) {
    tg <- V4_TARGETS$table16[[nm]]
    a  <- quantisation_adjoint(scen[[nm]], grid = "full")
    add(reg_row("Table 16", paste(nm, "XVA0"), a$U0, tg[["XVA0"]]))
    for (g in names(a$grads)) {
      ## r_f^- of the calm seller is zero on its branch; V4 prints "0".
      add(reg_row("Table 16", paste(nm, g), a$grads[[g]], tg[[g]]))
    }
  }
  for (nm in names(V4_TARGETS$table17)) {
    a <- quantisation_adjoint(scen[[nm]], grid = "frozen")
    add(reg_row("Table 17", paste(nm, "frozen vega"), a$grads[["sigma"]],
                V4_TARGETS$table17[[nm]]))
  }

  res <- do.call(rbind, out)
  attr(res, "seconds") <- proc.time()[["elapsed"]] - T0
  res
}

## Printing helper used by main.R.
report_regression <- function(res) {
  kn  <- !is.na(res$known)
  cat(sprintf("\n%d of %d regression checks pass (%.0f s); %d further rows test known V4 errata\n",
              sum(res$ok[!kn]), sum(!kn), attr(res, "seconds"), sum(kn)))
  bad <- res[!res$ok & !kn, ]
  if (nrow(bad) > 0) {
    cat("\nREGRESSIONS:\n")
    print(bad[, 1:7], row.names = FALSE)
  }
  ke <- res[kn, ]
  if (nrow(ke) > 0) {
    cat("\nKnown errata in V4 (not counted):\n")
    print(ke[, c("table", "item", "printed", "computed", "known")], row.names = FALSE)
  }
  invisible(all(res$ok[!kn]))
}

# -----------------------------------------------------------------------------

## Snapshot layer.  The printed digits of V4 cannot see a change below half a
## unit in the last digit, which is too coarse to certify that an optimisation
## leaves the numbers unchanged.  The snapshot stores every computed value of
## run_regression() at full precision; compare_snapshot() then requires
## agreement to a relative tolerance of 1e-9 (absolute floor 1e-15 for values
## that vanish on their branch).  The snapshot is written once, deliberately,
## with  Rscript main.R regress snapshot , and is not overwritten otherwise.

# -----------------------------------------------------------------------------
SNAPSHOT_FILE <- file.path("regression", "snapshot.rds")

write_snapshot <- function(res, file = SNAPSHOT_FILE) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(res[, c("table", "item", "value")], file)
  cat(sprintf("snapshot of %d values written to %s\n", nrow(res), file))
}

compare_snapshot <- function(res, file = SNAPSHOT_FILE, rtol = 1e-9, atol = 1e-15) {
  if (!file.exists(file)) {
    cat("no snapshot at ", file, "; skipping the snapshot layer\n", sep = "")
    return(invisible(NA))
  }
  ref <- readRDS(file)
  key <- function(d) paste(d$table, d$item, sep = " | ")
  m   <- match(key(ref), key(res))
  cur <- res$value[m]
  rel <- abs(cur - ref$value) / pmax(abs(ref$value), atol / rtol)
  ok  <- !is.na(m) & rel <= rtol
  cat(sprintf("snapshot: %d of %d values within %.0e relative (largest %.2e)\n",
              sum(ok), length(ok), rtol, max(rel, na.rm = TRUE)))
  if (!all(ok)) {
    print(data.frame(item = key(ref), snapshot = ref$value, current = cur,
                     rel = signif(rel, 3))[!ok, ], row.names = FALSE)
  }
  invisible(all(ok))
}
