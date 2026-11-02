# ------------------------------------------------------------------------------

##  25_buyer_convention.R
##  The buyer's crisis convention: the extended equation (default) against the
##  short-sale ban imposed on the control ("banned").  See docs/buyer_closeout.md
##  for the derivation and the case for each; the thesis decides which to report.
##
##  WHAT IS COMPARED.  The buyer's value is defined by the reflection (4.17),
##  f^-(y, z) = -f^+(-y, -z; -Vhat, -Zhat), whose close-out terms are the
##  correct ones for the holder (a CVA loss L_C (1 - alpha) Vhat on the
##  counterparty's default, no DVA).  The open point is the crisis: the
##  reflected hedge is short at every crisis node (Module 18), and the extended
##  equation finances that short by reverse repo although the repo market is
##  frozen and short sales are banned.  Under "banned" the hedge holds no stock
##  in the crisis and earns no reverse-repo carry.  The seller of a call never
##  reaches the bound, so its numbers must not move; this is checked.
##
##  PARTS.
##    1. values: XVA^- (and XVA^+) under both conventions over regime, alpha
##       and r_f^-, and the price band they imply, with inversions flagged;
##    2. the crisis effect for the buyer against the reverse-repo rate r_r^+
##       (the carry that the ban removes);
##    3. consistency: under "banned", quantisation against the PDE (both use
##       the same driver);
##    4. admissibility: under "banned" every node's stock position is
##       admissible; the nodes where the projection binds are counted;
##    5. sensitivities of the buyer under both conventions (eight, extended
##       AAD), with a same-scheme bump check.
##
##  OUTPUT (output/data/): buyerconv_values.csv, buyerconv_carry.csv,
##  buyerconv_pde_check.csv, buyerconv_admissibility.csv, buyerconv_greeks.csv.
##
##  USAGE.  Rscript main.R buyerconv [force]
# ------------------------------------------------------------------------------

BUYERCONV_SETTINGS <- list(
  conventions = c("extended", "banned"),
  regimes  = list(list(regime = "frozen0", mean_U = NA), list(regime = "frozen1", mean_U = NA),
                  list(regime = "dynamic", mean_U = 1.39), list(regime = "dynamic", mean_U = 0.3)),
  alphas   = c(0, 0.25, 0.5, 0.75, 1),
  rf_m     = c(0.05, 0.075, 0.10),
  levels   = list(c(50L, 100L), c(100L, 200L)),
  carry_rr_p = seq(0.03, 0.05, by = 0.005),
  pde      = c(801L, 1500L),
  pde_alphas = c(0, 0.5, 1),
  node_level = c(50L, 100L),
  greek_h  = c(S0 = 1e-4, sigma = 1e-4, rf_m = 1e-4, hI = 1e-4, hC = 1e-4,
               alpha = 1e-4, LI = 1e-4, LC = 1e-4),
  seed = 25L)

BUYERCONV_USED_BY <- paste("Ch5 / appendix: buyer crisis convention table (values, band,",
  "carry), consistency and admissibility checks, buyer sensitivities under both conventions")

conv_params <- function(p0, side, rg, alpha, rf_m, conv, ...) {
  p <- with_params(p0, side = side, beta_mode = rg$regime, alpha = alpha, rf_m = rf_m, ...)
  if (!is.na(rg$mean_U)) p <- with_params(p, mean_U = rg$mean_U)
  p$crisis_short <- conv
  p
}
conv_value <- function(p, levels) {
  q <- vapply(levels, function(lv) quantisation_bsde(p, lv[1], lv[2])$U0, 0)
  2 * q[2] - q[1]
}

## Part 1: one (regime, alpha, r_f^-) point, both sides, both conventions.
buyerconv_point <- function(cs, p0, s) {
  rg <- s$regimes[[cs$rg]]
  Vhat0 <- bs_value(0, p0$S0, p0)
  do.call(rbind, lapply(s$conventions, function(cv) {
    v <- vapply(c("+", "-"), function(sd)
      conv_value(conv_params(p0, sd, rg, cs$alpha, cs$rf_m, cv), s$levels), 0)
    data.frame(convention = cv, regime = regime_label(rg), alpha = cs$alpha, rf_m = cs$rf_m,
               XVA_plus = v[["+"]], XVA_minus = v[["-"]],
               XVA_minus_pct_Vhat0 = 100 * v[["-"]] / Vhat0,
               width_pct_Vhat0 = 100 * (v[["+"]] - v[["-"]]) / Vhat0,
               inverted = v[["+"]] - v[["-"]] < -INVERSION_TOL, stringsAsFactors = FALSE)
  }))
}

## Part 2: the buyer's crisis effect against r_r^+ (frozen crisis - calm).
buyerconv_carry <- function(cs, p0, s) {
  Vhat0 <- bs_value(0, p0$S0, p0)
  rg0 <- list(regime = "frozen0", mean_U = NA); rg1 <- list(regime = "frozen1", mean_U = NA)
  do.call(rbind, lapply(s$conventions, function(cv) {
    calm <- conv_value(conv_params(p0, "-", rg0, 0.5, 0.10, cv, rr_p = cs$rr_p), s$levels)
    cris <- conv_value(conv_params(p0, "-", rg1, 0.5, 0.10, cv, rr_p = cs$rr_p), s$levels)
    data.frame(convention = cv, rr_p = cs$rr_p, calm = calm, crisis = cris,
               effect_pct_Vhat0 = 100 * (cris - calm) / Vhat0, stringsAsFactors = FALSE)
  }))
}

## Part 3: quantisation against the PDE under "banned".
buyerconv_pde <- function(cs, p0, s, op) {
  rg <- s$regimes[[cs$rg]]
  p  <- conv_params(p0, cs$side, rg, cs$alpha, 0.10, "banned")
  q  <- conv_value(p, s$levels); v <- pde_value(p, op)
  ## the buyer's XVA vanishes in the frozen crisis at alpha = 1 under "banned"
  ## (no hedge, no carry, no close-out loss), so the relative difference is
  ## floored and the absolute one is reported in % of Vhat_0 as well
  data.frame(side = cs$side, regime = regime_label(rg), alpha = cs$alpha, quant = q, pde = v,
             rel_diff = rel_diff(q, v), diff_pct_Vhat0 = 100 * (q - v) / bs_value(0, p$S0, p),
             stringsAsFactors = FALSE)
}

## Part 4: admissibility of the buyer's hedge at every node under each convention.
buyerconv_admissibility <- function(cs, p0, s) {
  rg <- s$regimes[[cs$rg]]
  p  <- conv_params(p0, "-", rg, 0.5, 0.10, cs$convention)
  n  <- s$node_level[1]; N <- s$node_level[2]
  q  <- quantisation_bsde(p, n, N, keep = TRUE)
  out <- c(crisis_nodes = 0, inadmissible = 0, binding = 0, residual = 0)
  for (k in seq_len(n)) for (j in seq_along(q$regimes)) {
    b <- q$regimes[j]; if (b != 1) next
    x <- q$grids[[k]]; t <- (k - 1) * p$Tmat / n
    a <- portfolio_accounts(t, exp(x), b, q$U[[k]][, j], q$Z[[k]][, j], p)
    raw <- -(q$Z[[k]][, j] + bs_control(t, exp(x), p))      # reflected raw control
    out <- out + c(length(x), sum(!a$admissible), sum(raw < 0), 0)
    out[["residual"]] <- max(out[["residual"]], max(abs(a$residual)))
  }
  data.frame(convention = cs$convention, regime = regime_label(rg), n = n, N = N,
             crisis_nodes = out[["crisis_nodes"]], inadmissible_nodes = out[["inadmissible"]],
             projection_binds_at = if (cs$convention == "banned") out[["binding"]] else 0,
             max_abs_residual = out[["residual"]], stringsAsFactors = FALSE)
}

## Part 5: the buyer's eight sensitivities under each convention.
buyerconv_greeks <- function(cs, p0, s) {
  rg <- s$regimes[[cs$rg]]
  p  <- conv_params(p0, "-", rg, 0.5, 0.10, cs$convention)
  par <- names(s$greek_h)
  a  <- quantisation_adjoint(p, 50L, 100L, grid = "full", params = par)
  fd <- vapply(par, function(nm) central(p, nm, s$greek_h[[nm]], function(q) quant_value(q)), 0)
  data.frame(convention = cs$convention, regime = regime_label(rg), parameter = par,
             adjoint = a$grads[par], bump = fd, rel_diff = rel_diff(a$grads[par], fd),
             XVA_minus = a$U0, row.names = NULL, stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------

## run_buyer_convention(): the suite.

# -----------------------------------------------------------------------------
run_buyer_convention <- function(p0 = make_params(), s = BUYERCONV_SETTINGS, force = FALSE) {
  run_suite("buyerconv", "25_buyer_convention.R", p0 = p0, settings = s, seed = s$seed,
            used_by = BUYERCONV_USED_BY, force = force, run = function() {
    for (lv in c(s$levels, list(s$node_level))) invisible(get_nodes(p0, lv[1], lv[2]))
    nr <- length(s$regimes)
    g1 <- expand.grid(rg = seq_len(nr), alpha = s$alphas, rf_m = s$rf_m)
    vals <- do.call(rbind, par_lapply(seq_len(nrow(g1)), function(i) buyerconv_point(g1[i, ], p0, s)))
    g2 <- data.frame(rr_p = s$carry_rr_p)
    carry <- do.call(rbind, par_lapply(seq_len(nrow(g2)), function(i) buyerconv_carry(g2[i, , drop = FALSE], p0, s)))
    op <- pde_operator(p0, s$pde[1], s$pde[2])
    g3 <- expand.grid(side = c("+", "-"), rg = seq_len(nr), alpha = s$pde_alphas,
                      stringsAsFactors = FALSE)
    pde <- do.call(rbind, par_lapply(seq_len(nrow(g3)), function(i) buyerconv_pde(g3[i, ], p0, s, op)))
    g4 <- expand.grid(convention = s$conventions, rg = seq_len(nr), stringsAsFactors = FALSE)
    adm <- do.call(rbind, par_lapply(seq_len(nrow(g4)), function(i) buyerconv_admissibility(g4[i, ], p0, s)))
    grk <- do.call(rbind, par_lapply(seq_len(nrow(g4)), function(i) buyerconv_greeks(g4[i, ], p0, s)))
    list(values = vals, carry = carry, pde_check = pde, admissibility = adm, greeks = grk)
  })
}

## Console report used by main.R.
report_buyer_convention <- function(res) {
  v <- res$values
  seller_moved <- max(abs(v$XVA_plus[v$convention == "banned"] - v$XVA_plus[v$convention == "extended"]))
  cat(sprintf("\nseller unchanged by the convention: max |difference| = %.1e\n", seller_moved))
  h <- v[v$rf_m == 0.10 & v$alpha %in% c(0, 0.5, 1), ]
  h$cell <- sprintf("%5.2f / %+6.2f", h$XVA_minus_pct_Vhat0, h$width_pct_Vhat0)
  cat("\nbuyer XVA^- / band width, % of Vhat_0, rf- = 0.10\n")
  print(reshape(h[, c("convention", "regime", "alpha", "cell")], idvar = c("convention", "regime"),
                timevar = "alpha", direction = "wide"), row.names = FALSE)
  cat("\ninverted band points (of", sum(v$convention == "extended"), "per convention):",
      paste(sprintf("%s %d", c("extended", "banned"),
                    c(sum(v$inverted[v$convention == "extended"]), sum(v$inverted[v$convention == "banned"]))),
            collapse = "; "), "\n")
  cat("\nbuyer crisis effect (frozen crisis - calm, % of Vhat_0) against r_r^+\n")
  print(reshape(res$carry[, c("convention", "rr_p", "effect_pct_Vhat0")], idvar = "rr_p",
                timevar = "convention", direction = "wide"), row.names = FALSE, digits = 4)
  cat(sprintf("\n'banned': quantisation against PDE, max |rel diff| = %.1e, max |diff| = %.1e%% of Vhat_0\n",
              max(abs(res$pde_check$rel_diff)), max(abs(res$pde_check$diff_pct_Vhat0))))
  cat("\nadmissibility of the buyer's hedge at the crisis nodes\n")
  print(res$admissibility, row.names = FALSE)
  g <- res$greeks
  cat(sprintf("\nbuyer sensitivities: adjoint against same-scheme bumps, max |rel diff| = %.1e\n",
              max(abs(g$rel_diff))))
  w <- reshape(g[, c("convention", "regime", "parameter", "adjoint")],
               idvar = c("regime", "parameter"), timevar = "convention", direction = "wide")
  print(w[w$regime %in% c("frozen1", "dynamic lamU=1.39"), ], row.names = FALSE, digits = 4)
  cat(sprintf("\nsuite time %.0f s on %d core(s)\n", res$meta$seconds, res$meta$cores))
}
