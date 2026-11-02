# ------------------------------------------------------------------------------

##  21_regime.R
##  Suite (e): the effect of the status process.
##
##  PARTS.
##    1. heat map: XVA_0 on a 15 x 15 grid of the mean holding times
##       (lambda_U, lambda_V), log-spaced, both sides and both initial regimes
##       (one quantisation solve returns the root value in either starting
##       regime, column 1 for beta_0 = 0 and column 2 for beta_0 = 1).  The
##       convention is that of Section 4.4: lambda are mean holding times in
##       years and the switching rates are their reciprocals.
##    2. the regime effect against the funding spread r_f^- - r_r^-: crisis
##       minus calm and dynamic minus calm over a scan of r_f^-, with a linear
##       fit per 100 bp and the standing assumptions checked at every point.
##       At r_f^- = r_r^- the crisis only moves the financing of the stock from
##       the repo desk to the treasury desk at the same rate, so the effect must
##       vanish there.  The fit is a quadratic through the origin (the effect
##       is mildly convex) and the slope at zero spread.
##    3. model risk: at every heat-map point the dynamic value is compared with
##       three shortcuts a desk might use instead of the regime-switching model:
##         frozen_initial  the regime frozen at its initial value,
##         frozen_crisis   the crisis frozen throughout (the stress value),
##         mixture         the two frozen values mixed by the time-averaged
##                         probability of the crisis on [0, T],
##                         pbar = (1/T) int_0^T P(beta_t = 1 | beta_0) dt.
##    4. PDE spot checks of the heat map at its four corners and at the
##       calibrated point, both sides, both initial regimes.
##  Values are the quantisation solver extrapolated by Richardson from (50, 100)
##  and (100, 200).  The buyer's rows rest on the buyer close-out convention
##  (flagged); the buyer is regime-independent for the call (Module 18).
##
##  OUTPUT (output/data/): regime_heatmap.csv, regime_spread.csv,
##  regime_spread_fit.csv, regime_pde_check.csv.
##
##  USAGE.  Rscript main.R regime [force]
# ------------------------------------------------------------------------------

REGIME_SETTINGS <- list(
  lambda   = exp(seq(log(0.25), log(3), length.out = 15)),   # both axes, years
  sides    = c("+", "-"),
  alpha    = 0.5,
  rf_m     = c(0.075, 0.10),
  levels   = list(c(50L, 100L), c(100L, 200L)),
  spread_rf_m  = seq(0.05, 0.15, by = 0.005),
  spread_alpha = c(0, 0.5, 1),
  spread_lambda = c(U = 1.39, V = 0.99),     # the calibrated holding times
  pde = c(801L, 1500L),
  seed = 21L)

REGIME_USED_BY <- paste("Ch5 regime heat maps (seller, buyer; beta_0 = 0, 1);",
  "regime effect against the funding spread figure; frozen-versus-dynamic model-risk table")

## Richardson-extrapolated root values for every starting regime carried.
regime_root <- function(p, levels) {
  q <- lapply(levels, function(lv) quantisation_bsde(p, lv[1], lv[2], keep = TRUE))
  U <- 2 * q[[2]]$U[[1]][1, ] - q[[1]]$U[[1]][1, ]
  Z <- 2 * q[[2]]$Z[[1]][1, ] - q[[1]]$Z[[1]][1, ]
  list(U = U, Z = Z, regimes = q[[1]]$regimes)
}

## Time-averaged probability of the crisis on [0, T] given beta_0.
crisis_share <- function(p, b0) {
  s <- p$aU + p$aV; pi1 <- p$aU / s
  (pi1 * p$Tmat + (b0 - pi1) * (1 - exp(-s * p$Tmat)) / s) / p$Tmat
}

# -----------------------------------------------------------------------------

## Part 1 and 3: one (side, r_f^-) slice of the heat map.

# -----------------------------------------------------------------------------
regime_slice <- function(cs, p0, s) {
  p  <- with_params(p0, side = cs$side, alpha = s$alpha, rf_m = cs$rf_m)
  fz <- c(calm   = regime_root(with_params(p, beta_mode = "frozen0"), s$levels)$U,
          crisis = regime_root(with_params(p, beta_mode = "frozen1"), s$levels)$U)
  Vhat0 <- bs_value(0, p$S0, p)
  g <- expand.grid(lambda_U = s$lambda, lambda_V = s$lambda)
  do.call(rbind, lapply(seq_len(nrow(g)), function(i) {
    pd <- with_params(p, beta_mode = "dynamic", mean_U = g$lambda_U[i], mean_V = g$lambda_V[i])
    r  <- regime_root(pd, s$levels)
    pb <- c(crisis_share(pd, 0), crisis_share(pd, 1))
    XVA <- r$U
    froz_init <- fz[c("calm", "crisis")]
    mix <- (1 - pb) * fz[["calm"]] + pb * fz[["crisis"]]
    data.frame(side = cs$side, alpha = s$alpha, rf_m = cs$rf_m,
               lambda_U = g$lambda_U[i], lambda_V = g$lambda_V[i],
               aU = pd$aU, aV = pd$aV, beta0 = c(0L, 1L),
               XVA = XVA, XVA_pct_Vhat0 = 100 * XVA / Vhat0, Z0 = r$Z,
               crisis_share = pb, stationary_crisis = pd$aU / (pd$aU + pd$aV),
               frozen_calm = fz[["calm"]], frozen_crisis = fz[["crisis"]],
               frozen_initial = unname(froz_init), mixture = mix,
               risk_frozen_initial_pct_Vhat0 = 100 * (unname(froz_init) - XVA) / Vhat0,
               risk_frozen_crisis_pct_Vhat0  = 100 * (fz[["crisis"]] - XVA) / Vhat0,
               risk_mixture_pct_Vhat0        = 100 * (mix - XVA) / Vhat0,
               flag = if (cs$side == "-") "buyer close-out convention" else "",
               stringsAsFactors = FALSE)
  }))
}

# -----------------------------------------------------------------------------

## Part 2: the regime effect against the funding spread.

# -----------------------------------------------------------------------------
regime_spread_case <- function(cs, p0, s) {
  p <- with_params(p0, side = cs$side, alpha = cs$alpha, rf_m = cs$rf_m,
                   mean_U = s$spread_lambda[["U"]], mean_V = s$spread_lambda[["V"]])
  v <- vapply(c("frozen0", "frozen1", "dynamic"), function(bm)
    regime_root(with_params(p, beta_mode = bm), s$levels)$U[1], 0)
  Vhat0 <- bs_value(0, p$S0, p)
  data.frame(side = cs$side, alpha = cs$alpha, rf_m = cs$rf_m,
             spread = cs$rf_m - p$rr_m, calm = v[["frozen0"]], crisis = v[["frozen1"]],
             dynamic = v[["dynamic"]],
             effect_crisis_pct_Vhat0 = 100 * (v[["frozen1"]] - v[["frozen0"]]) / Vhat0,
             effect_dynamic_pct_Vhat0 = 100 * (v[["dynamic"]] - v[["frozen0"]]) / Vhat0,
             assumptions_ok = check_assumptions(p, verbose = FALSE),
             flag = if (cs$side == "-") "buyer close-out convention" else "",
             stringsAsFactors = FALSE)
}

## The effect vanishes at zero spread and is mildly convex in it: in the crisis
## r_f^- also enters the coefficient c1 of the value itself (4.41), so the
## carry compounds.  Reported: the slope at zero spread (one-sided difference
## over the first scan step), and a quadratic fit through the origin,
## effect = b1 s + b2 s^2, with its R^2; all per 100 bp.
regime_spread_fit <- function(d) {
  do.call(rbind, lapply(split(d, list(d$side, d$alpha), drop = TRUE), function(x) {
    x <- x[order(x$spread), ]
    out <- lapply(c("effect_crisis_pct_Vhat0", "effect_dynamic_pct_Vhat0"), function(y) {
      e <- x[[y]]; sp <- x$spread / 0.01                   # spread in units of 100 bp
      if (all(abs(e) < 1e-10)) return(c(slope0 = 0, b1 = 0, b2 = 0, r2 = NA))
      f <- lm(e ~ 0 + sp + I(sp^2))
      c(slope0 = (e[2] - e[1]) / (sp[2] - sp[1]), b1 = unname(coef(f)[1]),
        b2 = unname(coef(f)[2]), r2 = summary(f)$r.squared)
    })
    data.frame(side = x$side[1], alpha = x$alpha[1],
               crisis_slope_at_zero = out[[1]][["slope0"]], crisis_b1 = out[[1]][["b1"]],
               crisis_b2 = out[[1]][["b2"]], crisis_r2 = out[[1]][["r2"]],
               dynamic_slope_at_zero = out[[2]][["slope0"]], dynamic_b1 = out[[2]][["b1"]],
               dynamic_b2 = out[[2]][["b2"]], dynamic_r2 = out[[2]][["r2"]],
               at_zero_spread_crisis = x$effect_crisis_pct_Vhat0[which.min(abs(x$spread))],
               stringsAsFactors = FALSE)
  }))
}

# -----------------------------------------------------------------------------

## Part 4: PDE spot checks.

# -----------------------------------------------------------------------------
regime_pde_case <- function(cs, p0, s, op) {
  p <- with_params(p0, side = cs$side, alpha = s$alpha, rf_m = max(s$rf_m),
                   beta_mode = "dynamic", mean_U = cs$lambda_U, mean_V = cs$lambda_V)
  r <- regime_root(p, s$levels)
  v <- pde_solve(p, op)$value
  data.frame(side = cs$side, lambda_U = cs$lambda_U, lambda_V = cs$lambda_V,
             beta0 = c(0L, 1L), quant = r$U, pde = unname(v[c("0", "1")]),
             rel_diff = r$U / unname(v[c("0", "1")]) - 1, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------

## run_regime(): the suite.

# -----------------------------------------------------------------------------
run_regime <- function(p0 = make_params(), s = REGIME_SETTINGS, force = FALSE) {
  run_suite("regime", "21_regime.R", p0 = p0, settings = s, seed = s$seed,
            used_by = REGIME_USED_BY, force = force, run = function() {
    for (lv in s$levels) invisible(get_nodes(p0, lv[1], lv[2]))
    op <- pde_operator(p0, s$pde[1], s$pde[2])
    sl <- expand.grid(side = s$sides, rf_m = s$rf_m, stringsAsFactors = FALSE)
    heat <- do.call(rbind, par_lapply(seq_len(nrow(sl)), function(i) regime_slice(sl[i, ], p0, s)))
    sp <- expand.grid(side = s$sides, alpha = s$spread_alpha, rf_m = s$spread_rf_m,
                      stringsAsFactors = FALSE)
    spread <- do.call(rbind, par_lapply(seq_len(nrow(sp)), function(i) regime_spread_case(sp[i, ], p0, s)))
    L <- range(s$lambda)
    pc <- rbind(expand.grid(lambda_U = L, lambda_V = L),
                data.frame(lambda_U = p0$mean_U, lambda_V = p0$mean_V))
    pc <- merge(pc, data.frame(side = s$sides, stringsAsFactors = FALSE))
    pde <- do.call(rbind, par_lapply(seq_len(nrow(pc)), function(i) regime_pde_case(pc[i, ], p0, s, op)))
    list(heatmap = heat, spread = spread, spread_fit = regime_spread_fit(spread), pde_check = pde)
  })
}

## Console report used by main.R.
report_regime <- function(res) {
  h <- res$heatmap
  sel <- h$rf_m == 0.10 & h$side == "+"
  lam <- sort(unique(h$lambda_U)); show <- lam[c(1, 4, 8, 11, 15)]
  for (b in 0:1) {
    d <- h[sel & h$beta0 == b & h$lambda_U %in% show & h$lambda_V %in% show, ]
    cat(sprintf("\nseller XVA_0 (%% of Vhat_0), rf- = 0.10, beta_0 = %d; rows lambda_U, columns lambda_V\n", b))
    print(round(xtabs(XVA_pct_Vhat0 ~ round(lambda_U, 2) + round(lambda_V, 2), d), 2))
  }
  b <- h[h$side == "-", ]
  cat(sprintf("\nbuyer: largest range of XVA_0 over the heat map within an r_f^- slice %.1e (regime-independent)\n",
              max(tapply(b$XVA, b$rf_m, function(v) diff(range(v))))))
  f <- res$spread_fit
  cat("\nregime effect against the spread r_f^- - r_r^- (% of Vhat_0; effect = b1 s + b2 s^2, s in 100 bp)\n")
  print(f, row.names = FALSE, digits = 3)
  cat("\nmodel risk over the heat map, seller, rf- = 0.10 (shortcut - dynamic, % of Vhat_0)\n")
  m <- h[sel, ]
  print(do.call(rbind, lapply(0:1, function(b) {
    x <- m[m$beta0 == b, ]
    data.frame(beta0 = b,
               frozen_initial = sprintf("[%+.2f, %+.2f]", min(x$risk_frozen_initial_pct_Vhat0), max(x$risk_frozen_initial_pct_Vhat0)),
               frozen_crisis = sprintf("[%+.2f, %+.2f]", min(x$risk_frozen_crisis_pct_Vhat0), max(x$risk_frozen_crisis_pct_Vhat0)),
               mixture = sprintf("[%+.3f, %+.3f]", min(x$risk_mixture_pct_Vhat0), max(x$risk_mixture_pct_Vhat0)))
  })), row.names = FALSE)
  cat("\nPDE spot checks: max |quantisation / PDE - 1| =",
      sprintf("%.2e\n", max(abs(res$pde_check$rel_diff))))
  cat(sprintf("\nsuite time %.0f s on %d core(s)\n", res$meta$seconds, res$meta$cores))
}
