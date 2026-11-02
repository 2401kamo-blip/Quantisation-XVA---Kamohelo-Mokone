#-------------------------------------------------------------------------------

##  15_sensitivity_validation.R
##  Validation of the adjoint sensitivities of 14_quantisation_adjoint.R.
##
##  Every article on which the adjoint is based validates it the same way: by
##  bump-and-revalue.  Two different comparisons are made here, and they answer
##  different questions.
##
##  V8  DISCRETE.  The adjoint against central finite differences of the SAME
##      quantisation scheme (same n, same N, grid rebuilt for each bump).  The
##      adjoint is the exact derivative of that discrete quantity, so the two
##      must agree to the accuracy of the finite differences.  This checks that
##      the adjoint code is correct.
##
##  V9  CONTINUOUS.  The adjoint against finite differences of the partial
##      differential equation benchmark (both sides) and of the closed form
##      (seller, frozen regimes).  The difference is the discretisation error of
##      the Greeks at n = 50, N = 100, the counterpart for sensitivities of the
##      value errors reported for 13.
##
##  Two further tables complete the picture: vega with the grid frozen against
##  vega with the grid differentiated, and the cost of the adjoint against the
##  cost of bumping.
##
##------------------------------------------------------------------------------

SENS_BUMP_DISCRETE   <- c(S0 = 1e-4, sigma = 1e-4, rf_m = 1e-4, hI = 1e-4, hC = 1e-4)
SENS_BUMP_CONTINUOUS <- c(S0 = 1e-3, sigma = 1e-3, rf_m = 1e-3, hI = 1e-3, hC = 1e-3)

## Central differences of value_fn; the bump in S0 is relative, the others absolute.
bump_sensitivities <- function(p, value_fn, h = SENS_BUMP_DISCRETE) {
  sapply(names(h), function(nm) {
    up <- p; dn <- p
    if (nm == "S0") {
      up$S0 <- p$S0 * (1 + h[[nm]]); dn$S0 <- p$S0 * (1 - h[[nm]])
      return((value_fn(up) - value_fn(dn)) / (2 * p$S0 * h[[nm]]))
    }
    up[[nm]] <- p[[nm]] + h[[nm]]; dn[[nm]] <- p[[nm]] - h[[nm]]
    (value_fn(up) - value_fn(dn)) / (2 * h[[nm]])
  })
}

## The test scenarios: both sides, the three regimes, alpha = 0.5, r_f^- = 0.10.
sensitivity_scenarios <- function(alpha = 0.5, rf_m = 0.10, mean_U = 1.39,
                                  sides = c("+", "-"),
                                  regimes = c("frozen0", "frozen1", "dynamic")) {
  out <- list()
  for (sd in sides) for (bm in regimes) {
    p <- make_params(alpha = alpha, rf_m = rf_m, beta_mode = bm, mean_U = mean_U)
    p$side <- sd
    out[[sprintf("%s %s", if (sd == "+") "seller" else "buyer", bm)]] <- p
  }
  out
}

## A relative difference with an absolute floor, so that a sensitivity that is
## zero on the branch (for example r_f^- for the calm seller) does not produce
## a meaningless ratio.
rel_diff <- function(a, b, floor = 1e-8) (a - b) / pmax(abs(b), floor)

# ---- V8: adjoint against bump-and-revalue on the same scheme ------------------
validate_adjoint_discrete <- function(scen = sensitivity_scenarios(), tol = 1e-5) {
  do.call(rbind, lapply(names(scen), function(nm) {
    p  <- scen[[nm]]
    a  <- quantisation_adjoint(p, grid = "full")
    fd <- bump_sensitivities(p, quant_value, SENS_BUMP_DISCRETE)
    rd <- rel_diff(a$grads, fd)
    data.frame(scenario = nm, parameter = names(a$grads), adjoint = a$grads,
               bump = fd, rel_diff = rd, ok = abs(rd) < tol, row.names = NULL)
  }))
}

# ---- V9: adjoint against the continuous benchmarks ----------------------------
validate_adjoint_continuous <- function(scen = sensitivity_scenarios(
                                          sides = "+", regimes = c("frozen0", "frozen1", "dynamic"))) {
  do.call(rbind, lapply(names(scen), function(nm) {
    p   <- scen[[nm]]
    a   <- quantisation_adjoint(p, grid = "full")
    pde <- bump_sensitivities(p, function(q) pde_value(q), SENS_BUMP_CONTINUOUS)
    cf  <- if (identical(p$side, "+") && !identical(p$beta_mode, "dynamic")) {
      bump_sensitivities(p, function(q) linear_benchmark(q)$value, SENS_BUMP_CONTINUOUS)
    } else rep(NA_real_, length(pde))
    data.frame(scenario = nm, parameter = names(a$grads), adjoint = a$grads,
               pde_bump = pde, rel_diff_pde = rel_diff(a$grads, pde),
               closed_bump = cf, rel_diff_closed = rel_diff(a$grads, cf),
               row.names = NULL)
  }))
}

# ---- vega: grid frozen against grid differentiated ----------------------------
vega_frozen_vs_full <- function(scen = sensitivity_scenarios()) {
  do.call(rbind, lapply(names(scen), function(nm) {
    p  <- scen[[nm]]
    a  <- quantisation_adjoint(p, grid = "full")
    fz <- quantisation_adjoint(p, grid = "frozen")
    data.frame(scenario = nm, vega_full = a$grads[["sigma"]],
               vega_frozen = fz$grads[["sigma"]],
               grid_contribution = a$vega_grid,
               frozen_error = (fz$grads[["sigma"]] - a$grads[["sigma"]]) /
                              abs(a$grads[["sigma"]]))
  }))
}

# ---- cost: one valuation, the adjoint, and bump-and-revalue ------------------

## The grid is rebuilt from scratch in every case (the cache is cleared), so
## that the comparison charges each method for the grids it needs: one for the
## valuation and the adjoint, and one per bump of S0 or sigma for bumping.
sensitivity_cost <- function(p = sensitivity_scenarios()[[2]]) {
  clear <- clear_quant_caches
  el <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0 }
  clear(); t_val  <- el(quant_value(p))
  clear(); t_adj  <- el(quantisation_adjoint(p, grid = "full"))
  clear(); t_bump <- el(bump_sensitivities(p, quant_value, SENS_BUMP_DISCRETE))
  data.frame(method = c("one valuation", "valuation + adjoint (5 Greeks)",
                        "bump-and-revalue (5 Greeks, central)"),
             seconds = c(t_val, t_adj, t_bump),
             in_valuations = c(1, t_adj / t_val, t_bump / t_val))
}
