# ------------------------------------------------------------------------------
##  20_bands.R
##  Suite (d): the seller's and buyer's XVA, the price band, funding positions,
##  and the sensitivity to the repo rates and the credit intensities.
##
##  THE BAND.  The investor sells the claim at V^+ = Vhat + XVA^+ and buys it at
##  V^- = Vhat + XVA^-, the buyer's XVA being defined through the reflection of
##  (4.17).  The width of the band is XVA^+ - XVA^-; a negative width (the
##  buying price above the selling price) is flagged as an inverted band.  At
##  the benchmark parameters the band inverts in the calm regime for alpha < 1,
##  because the seller's DVA benefit (h_I = 0.20) exceeds the buyer's CVA cost
##  (h_C = 0.15); with L_I = L_C = 0 the two sides coincide.  The buyer's side
##  rests on the buyer close-out convention of the driver, an open
##  methodological decision, so every buyer quantity carries that flag.
##
##  PARTS.
##    1. band grid: both sides over alpha (fine grid), r_f^- and regime
##       (frozen calm, frozen crisis, dynamic at the listed lambda_U);
##    2. funding positions in currency: the six accounts at t = 0 at every band
##       point (portfolio_accounts() of 18), and expected funding profiles over
##       the tree, E[A^+] (lending to the treasury), E[A^-] (borrowing) and E[A],
##       the regime weighted by its exact marginal law, pre-default and
##       survival-weighted;
##    3. repo rates: a symmetric repo spread r_r^+/- = r0 -/+ d, then r_r^- and
##       r_r^+ alone, inside Proposition 14 (r_r^+ <= r_f^+ <= r_r^-); the
##       standing assumptions are checked at every point;
##    4. credit: h_I and h_C scaled on a grid, with mu_i = r_D + h_i moved too.
##  Values are the quantisation solver extrapolated by Richardson from
##  (50, 100) and (100, 200), as in Modules 17 to 19.
##
##  OUTPUT (output/data/): bands_grid.csv, bands_accounts_t0.csv,
##  bands_funding_profile.csv, bands_repo.csv, bands_credit.csv.
##
##  USAGE.  Rscript main.R bands [force]

# ------------------------------------------------------------------------------

BANDS_SETTINGS <- list(
  alphas   = seq(0, 1, by = 0.05),
  rf_m     = c(0.05, 0.075, 0.10),
  regimes  = list(list(regime = "frozen0", mean_U = NA), list(regime = "frozen1", mean_U = NA),
                  list(regime = "dynamic", mean_U = 1.39), list(regime = "dynamic", mean_U = 0.3)),
  levels   = list(c(50L, 100L), c(100L, 200L)),
  profile_level  = c(50L, 100L),
  profile_alphas = c(0, 0.5, 1), profile_rf_m = 0.10,
  sens_alpha = 0.5, sens_rf_m = 0.10,
  sens_regimes = list(list(regime = "frozen0", mean_U = NA), list(regime = "frozen1", mean_U = NA),
                      list(regime = "dynamic", mean_U = 1.39)),
  repo_spread = seq(0, 0.02, by = 0.0025),   # r_r^+/- = 0.05 -/+ d
  repo_minus  = seq(0.05, 0.07, by = 0.0025),
  repo_plus   = seq(0.03, 0.05, by = 0.0025),
  credit_scale = c(0.5, 0.75, 1, 1.5, 2),
  seed = 20L)

BANDS_USED_BY <- paste("Ch5 band figure (XVA+/- over alpha by r_f^- and regime);",
  "funding-position table and profile figure; repo and credit sensitivity tables")

BUYER_FLAG <- "buyer close-out convention"

## Where the two sides coincide (alpha = 1 at r_f^- = r_f^+) Richardson leaves
## noise of order 1e-17; a band counts as inverted beyond this tolerance.
INVERSION_TOL <- 1e-10

## Richardson-extrapolated (Y0, Z0) of the quantisation solver.
rich_YZ <- function(p, levels) {
  q <- lapply(levels, function(lv) quantisation_bsde(p, lv[1], lv[2]))
  c(Y0 = 2 * q[[2]]$U0 - q[[1]]$U0, Z0 = 2 * q[[2]]$Z0 - q[[1]]$Z0)
}

regime_params <- function(p, rg) {
  p <- with_params(p, beta_mode = rg$regime)
  if (!is.na(rg$mean_U)) p <- with_params(p, mean_U = rg$mean_U)
  p
}
regime_label <- function(rg) paste0(rg$regime, if (is.na(rg$mean_U)) "" else sprintf(" lamU=%.2f", rg$mean_U))

## Both sides at one parameter point: XVA, prices, width, t = 0 accounts.
band_point <- function(p, levels) {
  b0 <- if (identical(p$beta_mode, "frozen1")) 1 else 0
  Vhat0 <- bs_value(0, p$S0, p)
  sides <- lapply(c("+", "-"), function(sd) {
    ps <- with_params(p, side = sd)
    yz <- rich_YZ(ps, levels)
    a  <- portfolio_accounts(0, p$S0, b0, yz[["Y0"]], yz[["Z0"]], ps)
    list(XVA = yz[["Y0"]], Z0 = yz[["Z0"]], acc = a)
  })
  list(XVA_plus = sides[[1]]$XVA, XVA_minus = sides[[2]]$XVA, Vhat0 = Vhat0,
       acc = rbind(cbind(side = "+", sides[[1]]$acc), cbind(side = "-", sides[[2]]$acc)))
}

band_row <- function(keys, bp) {
  cbind(keys, XVA_plus = bp$XVA_plus, XVA_minus = bp$XVA_minus,
        XVA_plus_pct_Vhat0 = 100 * bp$XVA_plus / bp$Vhat0,
        XVA_minus_pct_Vhat0 = 100 * bp$XVA_minus / bp$Vhat0,
        price_sell = bp$Vhat0 + bp$XVA_plus, price_buy = bp$Vhat0 + bp$XVA_minus,
        width = bp$XVA_plus - bp$XVA_minus,
        width_pct_Vhat0 = 100 * (bp$XVA_plus - bp$XVA_minus) / bp$Vhat0,
        inverted = bp$XVA_plus - bp$XVA_minus < -INVERSION_TOL, flag = BUYER_FLAG,
        stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------

## Part 1 and the t = 0 accounts of Part 2.

# -----------------------------------------------------------------------------
bands_grid_case <- function(cs, p0, s) {
  rg <- s$regimes[[cs$rg]]
  p  <- regime_params(with_params(p0, alpha = cs$alpha, rf_m = cs$rf_m), rg)
  bp <- band_point(p, s$levels)
  keys <- data.frame(regime = regime_label(rg), alpha = cs$alpha, rf_m = cs$rf_m,
                     stringsAsFactors = FALSE)
  acc <- bp$acc
  acc_long <- do.call(rbind, lapply(ACCOUNTS, function(a)
    data.frame(keys[rep(1, nrow(acc)), ], side = acc$side, account = a,
               units = acc[[paste0("units_", a)]], value = acc[[paste0("value_", a)]],
               value_pct_Vhat0 = 100 * acc[[paste0("value_", a)]] / bp$Vhat0,
               admissible = acc$admissible, residual = acc$residual,
               flag = ifelse(acc$side == "-", BUYER_FLAG, ""),
               row.names = NULL, stringsAsFactors = FALSE)))
  list(grid = band_row(keys, bp), acc = acc_long)
}

# -----------------------------------------------------------------------------

## Part 2: expected funding profile over the tree.

# -----------------------------------------------------------------------------
funding_profile <- function(p, n, N) {
  q    <- quantisation_bsde(p, n, N, keep = TRUE)
  regs <- q$regimes
  dyn  <- identical(p$beta_mode, "dynamic")
  do.call(rbind, lapply(seq_len(n), function(k) {
    t  <- (k - 1) * p$Tmat / n
    x  <- q$grids[[k]]; pk <- q$probs[[k]]
    pr <- if (dyn) { c1 <- beta_prob_crisis(p, t, start = 0); c(1 - c1, c1) } else 1
    acc <- lapply(seq_along(regs), function(j)
      portfolio_accounts(t, exp(x), regs[j], q$U[[k]][, j], q$Z[[k]][, j], p))
    E <- function(f) sum(vapply(seq_along(regs), function(j) pr[j] * sum(pk * f(acc[[j]])), 0))
    surv <- exp(-(p$hI + p$hC) * t)
    data.frame(t = t,
               E_funding = E(function(a) a$value_funding),
               E_lending = E(function(a) pmax(a$value_funding, 0)),
               E_borrowing = E(function(a) pmax(-a$value_funding, 0)),
               E_stock = E(function(a) a$value_stock), E_repo = E(function(a) a$value_repo),
               E_collateral = E(function(a) a$value_collateral),
               P_crisis = if (dyn) pr[2] else as.numeric(regs[1] == 1),
               survival = surv)
  }))
}

bands_profile_case <- function(cs, p0, s) {
  rg <- s$regimes[[cs$rg]]
  p  <- regime_params(with_params(p0, side = cs$side, alpha = cs$alpha, rf_m = s$profile_rf_m), rg)
  pr <- funding_profile(p, s$profile_level[1], s$profile_level[2])
  cbind(data.frame(side = cs$side, regime = regime_label(rg), alpha = cs$alpha,
                   rf_m = s$profile_rf_m, flag = if (cs$side == "-") BUYER_FLAG else "",
                   stringsAsFactors = FALSE)[rep(1, nrow(pr)), ], pr, row.names = NULL)
}

# -----------------------------------------------------------------------------

## Parts 3 and 4: sensitivities.

# -----------------------------------------------------------------------------
sens_point <- function(p, keys, s) {
  ok <- check_assumptions(p, verbose = FALSE)
  bp <- band_point(p, s$levels)
  fund <- bp$acc$value_funding
  cbind(band_row(keys, bp), assumptions_ok = ok,
        funding_t0_seller = fund[bp$acc$side == "+"],
        funding_t0_buyer = fund[bp$acc$side == "-"])
}

bands_repo_case <- function(cs, p0, s) {
  rg <- s$sens_regimes[[cs$rg]]
  p  <- regime_params(with_params(p0, alpha = s$sens_alpha, rf_m = s$sens_rf_m,
                                  rr_p = cs$rr_p, rr_m = cs$rr_m), rg)
  sens_point(p, data.frame(design = cs$design, regime = regime_label(rg),
                           rr_p = cs$rr_p, rr_m = cs$rr_m, stringsAsFactors = FALSE), s)
}

bands_credit_case <- function(cs, p0, s) {
  rg <- s$sens_regimes[[cs$rg]]
  hI <- p0$hI * cs$sI; hC <- p0$hC * cs$sC
  p  <- regime_params(with_params(p0, alpha = s$sens_alpha, rf_m = s$sens_rf_m,
                                  hI = hI, hC = hC, muI = p0$rD + hI, muC = p0$rD + hC), rg)
  sens_point(p, data.frame(regime = regime_label(rg), scale_hI = cs$sI, scale_hC = cs$sC,
                           hI = hI, hC = hC, stringsAsFactors = FALSE), s)
}

# -----------------------------------------------------------------------------

## run_bands(): the suite.

# -----------------------------------------------------------------------------
run_bands <- function(p0 = make_params(), s = BANDS_SETTINGS, force = FALSE) {
  run_suite("bands", "20_bands.R", p0 = p0, settings = s, seed = s$seed,
            used_by = BANDS_USED_BY, force = force, run = function() {
    for (lv in c(s$levels, list(s$profile_level))) invisible(get_nodes(p0, lv[1], lv[2]))
    nr <- length(s$regimes); ns <- length(s$sens_regimes)
    g  <- expand.grid(alpha = s$alphas, rf_m = s$rf_m, rg = seq_len(nr))
    pa <- par_lapply(seq_len(nrow(g)), function(i) bands_grid_case(g[i, ], p0, s))
    pf <- expand.grid(side = c("+", "-"), alpha = s$profile_alphas, rg = seq_len(nr),
                      stringsAsFactors = FALSE)
    prof <- par_lapply(seq_len(nrow(pf)), function(i) bands_profile_case(pf[i, ], p0, s))
    r0 <- p0$rf_p
    rp <- rbind(
      data.frame(design = "symmetric spread", rr_p = r0 - s$repo_spread, rr_m = r0 + s$repo_spread),
      data.frame(design = "r_r^- alone", rr_p = p0$rr_p, rr_m = s$repo_minus),
      data.frame(design = "r_r^+ alone", rr_p = s$repo_plus, rr_m = p0$rr_m))
    rp <- merge(rp, data.frame(rg = seq_len(ns)))
    repo <- par_lapply(seq_len(nrow(rp)), function(i) bands_repo_case(rp[i, ], p0, s))
    cr <- expand.grid(sI = s$credit_scale, sC = s$credit_scale, rg = seq_len(ns))
    cred <- par_lapply(seq_len(nrow(cr)), function(i) bands_credit_case(cr[i, ], p0, s))
    list(grid = do.call(rbind, lapply(pa, `[[`, "grid")),
         accounts_t0 = do.call(rbind, lapply(pa, `[[`, "acc")),
         funding_profile = do.call(rbind, prof),
         repo = do.call(rbind, repo), credit = do.call(rbind, cred))
  })
}

## Console report used by main.R.
report_bands <- function(res) {
  g <- res$grid
  h <- g[g$alpha %in% c(0, 0.5, 1), ]
  cat("\nthe band, % of Vhat_0 (XVA+ / XVA- / width)\n")
  h$cell <- sprintf("%6.2f / %5.2f / %+5.2f", h$XVA_plus_pct_Vhat0, h$XVA_minus_pct_Vhat0,
                    h$width_pct_Vhat0)
  w <- reshape(h[, c("regime", "rf_m", "alpha", "cell")], idvar = c("regime", "rf_m"),
               timevar = "alpha", direction = "wide")
  names(w) <- sub("^cell\\.", "alpha=", names(w))
  print(w, row.names = FALSE)
  inv <- g[g$inverted, ]
  cat(sprintf("\ninverted band at %d of %d grid points (%s); largest inversion %.3f%% of Vhat_0\n",
              nrow(inv), nrow(g), paste(unique(inv$regime), collapse = ", "),
              if (nrow(inv)) -min(inv$width_pct_Vhat0) else 0))
  a <- res$accounts_t0
  a <- a[a$account == "funding" & a$alpha %in% c(0, 0.5, 1) & a$rf_m == 0.10, ]
  cat("\nfunding account at t = 0, rf- = 0.10 (currency; buyer reflected)\n")
  print(reshape(a[, c("regime", "side", "alpha", "value")], idvar = c("regime", "side"),
                timevar = "alpha", direction = "wide"), row.names = FALSE, digits = 4)
  pf <- res$funding_profile
  pk <- aggregate(cbind(E_lending, E_borrowing) ~ side + regime + alpha, pf, max)
  cat("\npeak expected lending / borrowing over the life (rf- = 0.10)\n")
  print(pk, row.names = FALSE, digits = 3)
  rp <- res$repo
  cat("\nrepo sensitivity (seller / buyer XVA in % of Vhat_0)\n")
  print(rp[, c("design", "regime", "rr_p", "rr_m", "XVA_plus_pct_Vhat0",
               "XVA_minus_pct_Vhat0", "assumptions_ok")][rp$regime != "dynamic lamU=1.39", ],
        row.names = FALSE, digits = 4)
  cr <- res$credit
  cat("\ncredit sensitivity, frozen crisis: XVA+ in % of Vhat_0 (rows h_I scale, columns h_C scale)\n")
  c1 <- cr[cr$regime == "frozen1", ]
  print(xtabs(XVA_plus_pct_Vhat0 ~ scale_hI + scale_hC, c1), digits = 4)
  cat("calm: band width in % of Vhat_0\n")
  c0 <- cr[cr$regime == "frozen0", ]
  print(xtabs(width_pct_Vhat0 ~ scale_hI + scale_hC, c0), digits = 3)
  cat(sprintf("\nsuite time %.0f s on %d core(s)\n", res$meta$seconds, res$meta$cores))
}
