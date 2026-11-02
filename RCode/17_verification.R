# -----------------------------------------------------------------------------

##  17_verification.R
##  Suite (a): verification of the quantisation solver (13) against the closed
##  form (08) and the PDE (07).
##
##  WHAT IS COMPARED.  For every case (side, regime, alpha, r_f^-, and lambda_U
##  for the dynamic regime) the solver is run on the jointly refined grids
##  (n, N) = (25, 50), (50, 100), (100, 200), and its value XVA_0 = Ubreve_0 and
##  control Zbreve_0 are compared with two references:
##    * the closed form, for the seller in a frozen regime, and only where its
##      assumed branch holds (linear_branch_check() in 08).  Its control is
##      Zbreve_0 = sigma S0 dXVA_0/dS0, by a central difference of the closed
##      form in S0, and both are extrapolated in the number of quadrature time
##      nodes (see cf_value_control());
##    * the PDE on a 1601 x 3000 grid, for every case.  The PDE on the standard
##      801 x 1500 grid is also recorded, so that the gap between the two PDE
##      grids (and between the fine PDE and the closed form) shows how far the
##      reference itself can be trusted.
##  The primary reference is the closed form where it is valid and the fine PDE
##  otherwise.
##
##  RICHARDSON EXTRAPOLATION.  Joint refinement doubles n, so the error of the
##  explicit scheme falls as O(Delta t) and 2 Q(2n) - Q(n) removes the leading
##  term.  The observed order log2((Q1 - Q2) / (Q2 - Q3)) checks that premise.
##
##  OUTPUT (output/data/):
##    verification_values.csv   one row per case, quantity and method (tidy)
##    verification_summary.csv  one row per case and quantity: references,
##                              errors at each level, observed order, Richardson
##    verification.rds          both, with the run metadata
##  XVA is reported in currency units and as a percentage of Vhat_0.
##
##  USAGE.  Rscript main.R verify [force]

#------------------------------------------------------------------------------- 

VERIFY_SETTINGS <- list(
  levels  = list(c(25L, 50L), c(50L, 100L), c(100L, 200L)),  # jointly refined (n, N)
  alphas  = c(0, 0.5, 1),
  rf_m    = 0.10,                  # the headline stress setting
  mean_U  = c(1.39, 0.3),          # dynamic cases: calibrated, and fast switching
  sides   = c("+", "-"),
  pde_ref = c(1601L, 3000L),       # (nx, nt) of the reference PDE
  pde_std = c(801L, 1500L),        # (nx, nt) of the standard PDE
  cf_bump = 1e-4,                  # relative bump in S0 for the closed-form control
  cf_ns   = c(1600L, 3200L, 6400L),  # closed-form time nodes, Richardson in ns
  seed    = 17L)                   # the suite is deterministic; recorded for form

VERIFY_USED_BY <- "Ch5 verification table (values and Z0 errors); convergence figure"

# -----------------------------------------------------------------------------

## The cases, built from the base parameter object p0.

# -----------------------------------------------------------------------------
verification_cases <- function(p0, s = VERIFY_SETTINGS) {
  grid <- rbind(
    expand.grid(side = s$sides, regime = c("frozen0", "frozen1"), alpha = s$alphas,
                rf_m = s$rf_m, mean_U = NA_real_, stringsAsFactors = FALSE),
    expand.grid(side = s$sides, regime = "dynamic", alpha = s$alphas,
                rf_m = s$rf_m, mean_U = s$mean_U, stringsAsFactors = FALSE))
  grid$case <- sprintf("%s %s%s alpha=%.2f rf-=%.3f",
                       ifelse(grid$side == "+", "seller", "buyer"), grid$regime,
                       ifelse(is.na(grid$mean_U), "", sprintf(" lamU=%.2f", grid$mean_U)),
                       grid$alpha, grid$rf_m)
  grid
}

case_params <- function(p0, cs) {
  p <- with_params(p0, side = cs$side, beta_mode = cs$regime, alpha = cs$alpha,
                   rf_m = cs$rf_m)
  if (!is.na(cs$mean_U)) p <- with_params(p, mean_U = cs$mean_U)
  p
}

## Closed-form value and control; NA where the closed form does not apply.
##
## The outer time quadrature of 08 is a trapezoid rule whose integrand behaves
## like 1/sqrt(T - s) near maturity once it is differentiated in S0, so the
## control converges only at first order in the number of time nodes ns: at
## the default ns = 400 its error is about 0.1 per cent, while the value is
## accurate to a few parts in 1e6.  The reference is therefore computed at
## ns = 1600, 3200 and 6400 and extrapolated with Richardson's rule of order
## one on the last two; the observed order from the three is recorded, and so
## is the size of the extrapolation as an estimate of the quadrature error.
cf_value_control <- function(p, h, ns = VERIFY_SETTINGS$cf_ns) {
  na <- c(value = NA, control = NA, branch_err = NA, ok = NA, value_quad_err = NA,
          control_quad_err = NA, value_ns_order = NA, control_ns_order = NA)
  if (!identical(p$side, "+") || identical(p$beta_mode, "dynamic")) return(na)
  lb <- linear_benchmark(p, check = TRUE)
  br <- branch_name(p, lb)
  pu <- with_params(p, S0 = p$S0 * (1 + h)); pd <- with_params(p, S0 = p$S0 * (1 - h))
  v  <- vapply(ns, function(k) linear_benchmark(p, branch = br, ns = k)$value, 0)
  z  <- vapply(ns, function(k) p$sigma * p$S0 *
          (linear_benchmark(pu, branch = br, ns = k)$value -
           linear_benchmark(pd, branch = br, ns = k)$value) / (2 * p$S0 * h), 0)
  m  <- length(ns)
  ext   <- function(x) 2 * x[m] - x[m - 1]
  order <- function(x) log2((x[m - 2] - x[m - 1]) / (x[m - 1] - x[m]))
  c(value = ext(v), control = ext(z), branch_err = lb$branch_err,
    ok = as.numeric(lb$branch_ok),
    value_quad_err = ext(v) - v[m], control_quad_err = ext(z) - z[m],
    value_ns_order = order(v), control_ns_order = order(z))
}
branch_name <- function(p, lb) if (lb$rho_f == p$rf_p) "rf_plus" else "rf_minus"

# -----------------------------------------------------------------------------

## One case: all methods.  Returns long-format rows.

# -----------------------------------------------------------------------------
verify_case <- function(cs, p0, ops, s) {
  p   <- case_params(p0, cs)
  row <- function(method, n, N, quantity, value, seconds)
    data.frame(case = cs$case, side = cs$side, regime = cs$regime, alpha = cs$alpha,
               rf_m = cs$rf_m, mean_U = cs$mean_U, method = method, n = n, N = N,
               quantity = quantity, value = value, seconds = seconds,
               stringsAsFactors = FALSE)
  el  <- function(expr) { t0 <- proc.time()[["elapsed"]]; v <- force(expr)
                          list(v = v, s = proc.time()[["elapsed"]] - t0) }
  out <- list()

  r <- el(pde_value_control(p, ops$ref))
  out[[1]] <- row("pde_ref", s$pde_ref[1], s$pde_ref[2], c("U0", "Z0"), r$v, r$s)
  r <- el(pde_value_control(p, ops$std))
  out[[2]] <- row("pde_std", s$pde_std[1], s$pde_std[2], c("U0", "Z0"), r$v, r$s)
  r <- el(cf_value_control(p, s$cf_bump))
  out[[3]] <- row("closed_form", NA, NA, c("U0", "Z0"), r$v[c("value", "control")], r$s)
  out[[3]]$branch_err <- r$v[["branch_err"]]; out[[3]]$branch_ok <- r$v[["ok"]]
  out[[3]]$quad_err <- r$v[c("value_quad_err", "control_quad_err")]
  out[[3]]$ns_order <- r$v[c("value_ns_order", "control_ns_order")]
  for (lv in s$levels) {
    r <- el(quantisation_bsde(p, n = lv[1], N = lv[2]))
    out[[length(out) + 1L]] <- row("quantisation", lv[1], lv[2], c("U0", "Z0"),
                                   c(r$v$U0, r$v$Z0), r$s)
  }
  do.call(rbind, lapply(out, function(d) {
    if (is.null(d$branch_err)) { d$branch_err <- NA_real_; d$branch_ok <- NA_real_
                                 d$quad_err <- NA_real_; d$ns_order <- NA_real_ }
    d
  }))
}

# -----------------------------------------------------------------------------

## Summary: references, errors at each level, observed order, Richardson.

# -----------------------------------------------------------------------------
verification_summary <- function(vals, Vhat0, s = VERIFY_SETTINGS) {
  key <- unique(vals[, c("case", "side", "regime", "alpha", "rf_m", "mean_U", "quantity")])
  rows <- lapply(seq_len(nrow(key)), function(i) {
    k  <- key[i, ]
    v  <- vals[vals$case == k$case & vals$quantity == k$quantity, ]
    get <- function(m) v$value[v$method == m]
    cf  <- get("closed_form")
    cf_ok <- isTRUE(v$branch_ok[v$method == "closed_form"] == 1)
    pde <- get("pde_ref")
    cfq <- v$quad_err[v$method == "closed_form"]
    ref <- if (cf_ok) cf else pde
    q   <- vapply(s$levels, function(lv) v$value[v$method == "quantisation" & v$n == lv[1]], 0)
    nl  <- length(q)
    rich <- 2 * q[-1] - q[-nl]                      # Richardson on consecutive levels
    ord  <- if (nl >= 3) log2((q[1] - q[2]) / (q[2] - q[3])) else NA_real_
    lab  <- vapply(s$levels, function(lv) sprintf("%d_%d", lv[1], lv[2]), "")
    d <- data.frame(k, reference = if (cf_ok) "closed_form" else "pde_ref",
                    ref = ref, closed_form = if (cf_ok) cf else NA_real_,
                    pde_ref = pde, pde_std = get("pde_std"),
                    cf_quad_err_rel = if (cf_ok) cfq / cf else NA_real_,
                    pde_ref_vs_cf = if (cf_ok) pde / cf - 1 else NA_real_,
                    pde_std_vs_ref = get("pde_std") / pde - 1,
                    observed_order = ord, stringsAsFactors = FALSE)
    for (j in seq_len(nl)) {
      d[[paste0("quant_", lab[j])]]   <- q[j]
      d[[paste0("relerr_", lab[j])]]  <- q[j] / ref - 1
    }
    for (j in seq_along(rich)) {
      d[[paste0("rich_", lab[j + 1])]]        <- rich[j]
      d[[paste0("relerr_rich_", lab[j + 1])]] <- rich[j] / ref - 1
    }
    if (k$quantity == "U0") {
      d$ref_pct_Vhat0 <- 100 * ref / Vhat0
    } else d$ref_pct_Vhat0 <- NA_real_
    d
  })
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------

## run_verification(): the suite.

# -----------------------------------------------------------------------------
run_verification <- function(p0 = make_params(), s = VERIFY_SETTINGS, force = FALSE) {
  run_suite("verification", "17_verification.R", p0 = p0, settings = s,
            seed = s$seed, used_by = VERIFY_USED_BY, force = force, run = function() {
    ## Shared, case-independent objects are built once here, before the cases
    ## are spread over the cores: the PDE operators (with their reference pair)
    ## and the quantisation trees and node caches of the three levels.
    t0  <- proc.time()[["elapsed"]]
    ops <- list(ref = pde_operator(p0, s$pde_ref[1], s$pde_ref[2]),
                std = pde_operator(p0, s$pde_std[1], s$pde_std[2]))
    t_ops <- proc.time()[["elapsed"]] - t0
    trees <- do.call(rbind, lapply(s$levels, function(lv) {
      t1 <- proc.time()[["elapsed"]]
      rq <- get_rmq(p0, lv[1], lv[2]); invisible(get_nodes(p0, lv[1], lv[2]))
      data.frame(n = lv[1], N = lv[2], tree_seconds = rq$seconds,
                 build_seconds = proc.time()[["elapsed"]] - t1,
                 lloyd_fallbacks = sum(rq$method == "lloyd"),
                 var_rel_err = rq$var_rel_err, mean_err = rq$mean_err)
    }))
    cases <- verification_cases(p0, s)
    vals  <- do.call(rbind, par_lapply(seq_len(nrow(cases)), function(i)
               verify_case(cases[i, ], p0, ops, s)))
    Vhat0 <- bs_value(0, p0$S0, p0)
    vals$value_pct_Vhat0 <- ifelse(vals$quantity == "U0", 100 * vals$value / Vhat0, NA_real_)
    list(values = vals, summary = verification_summary(vals, Vhat0, s),
         trees = trees,
         setup = data.frame(item = c("Vhat0", "pde_operators_seconds"),
                            value = c(Vhat0, t_ops)))
  })
}

## Console report used by main.R: one compact table per quantity.
report_verification <- function(res) {
  sm  <- res$summary
  lv  <- sprintf("%d_%d", vapply(res$meta$settings$levels, `[`, 0L, 1),
                 vapply(res$meta$settings$levels, `[`, 0L, 2))
  pct <- function(x) sprintf("%+.4f%%", 100 * x)
  for (qty in c("U0", "Z0")) {
    d <- sm[sm$quantity == qty, ]
    out <- data.frame(case = d$case, ref = ifelse(d$reference == "closed_form", "CF", "PDE"),
                      value = signif(d$ref, 6))
    for (l in lv) out[[l]] <- pct(d[[paste0("relerr_", l)]])
    out[[paste("rich", lv[length(lv)])]] <- pct(d[[paste0("relerr_rich_", lv[length(lv)])]])
    out$order <- round(d$observed_order, 2)
    cat(sprintf("\n%s: relative error of the quantisation solver against the reference\n", qty))
    print(out, row.names = FALSE)
  }
  cat("\nreference checks: max |fine PDE / closed form - 1| =",
      sprintf("%.2e", max(abs(sm$pde_ref_vs_cf), na.rm = TRUE)),
      "; max |standard PDE / fine PDE - 1| (U0) =",
      sprintf("%.2e", max(abs(sm$pde_std_vs_ref[sm$quantity == "U0"]))),
      "(Z0) =", sprintf("%.2e\n", max(abs(sm$pde_std_vs_ref[sm$quantity == "Z0"]))))
  cat("\ntrees\n"); print(res$trees, row.names = FALSE)
  cat(sprintf("\nsuite time %.0f s on %d core(s)\n", res$meta$seconds, res$meta$cores))
}
