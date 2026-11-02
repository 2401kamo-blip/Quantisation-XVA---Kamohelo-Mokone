# ------------------------------------------------------------------------------

##  19_shapley.R
##  Suite (c): Shapley decomposition of XVA_0 into its drivers.
##
##  THE GAME.  The players are features of the model.  For a coalition S of
##  features switched ON (the others OFF) the worth v(S) is XVA_0 of that
##  configuration.  The Shapley value of feature i is its marginal contribution
##  v(S + i) - v(S) averaged over all orderings of the features: 4! = 24
##  orderings for the main decomposition, 5! = 120 for the extended one.  The
##  orderings are enumerated literally.  By efficiency
##      XVA_0 (all on) = v(empty) + sum_i phi_i,
##  where the baseline v(empty) is what remains with every feature off: the
##  carry of the funding, repo and collateral accounts at the base rates and
##  the survival discounting at h_I, h_C (the default intensities are not
##  features; they stay on throughout).
##
##  WHAT "OFF" MEANS.
##    funding_spread  r_f^- = r_f^+ : borrowing from the treasury costs what
##                    lending earns; the asymmetry of the funding account is
##                    removed and nothing else is changed.
##    dva             L_I = 0 : the hedger (the investor I) still defaults at
##                    h_I, but its close-out carries no loss, theta_I = Vhat.
##    cva             L_C = 0 : the same for the counterparty, theta_C = Vhat.
##    regime          beta == 0 : the market stays in the normal regime; repo
##                    and short sales are available throughout.  For a dynamic
##                    target this removes the possibility of a crisis; for the
##                    frozen-crisis target it replaces the crisis by the calm
##                    regime.
##    collateral      alpha = 0 (extended decomposition only): no collateral
##                    account; the close-out exposure is the full Vhat.
##  ON means the value of the target case.
##
##  NOTES.  At r_f^- = r_f^+ = r_r^- the crisis only moves the financing of
##  the stock from the repo desk to the treasury desk at the same rate, so the
##  regime has no effect without the funding spread (V4, Table 6); Shapley
##  shares that interaction equally between the two features, and the spread
##  of the marginal contributions across orderings records it.  For the seller
##  of a call Vhat > 0, so theta_C = Vhat whatever L_C, and the CVA share is
##  zero; for the buyer the DVA share is zero likewise.  The buyer's values
##  rest on the buyer close-out convention of the driver (flagged).
##
##  VALUES.  Each v(S) is the quantisation solver extrapolated by Richardson
##  from (50, 100) and (100, 200), as in Modules 17 and 18.
##
##  OUTPUT (output/data/):
##    shapley_values.csv      phi per case and feature, in currency, in % of
##                            Vhat_0 and in % of XVA_0 - v(empty)
##    shapley_coalitions.csv  v(S) for every coalition of every case
##    shapley_orderings.csv   the marginal contributions in every ordering
##    shapley_dividends.csv   the Harsanyi dividends d(S) of every coalition
##
##  HARSANYI DIVIDENDS.  d(S) = sum over T in S of (-1)^(|S| - |T|) v(T) is the
##  part of v(S) not explained by its sub-coalitions: d({i}) is the stand-alone
##  effect of feature i, d({i, j}) the pure interaction of i and j, and so on.
##  Shapley splits each dividend equally among its members,
##      phi_i = sum over S containing i of d(S) / |S|,
##  which is checked.  Where two features are complements (the regime acts only
##  through the funding spread), Shapley's equal split is an axiom, not an
##  economic attribution; the dividend shows the joint term explicitly.
##
##  USAGE.  Rscript main.R shapley [force]
# ------------------------------------------------------------------------------

SHAPLEY_SETTINGS <- list(
  features   = c("funding_spread", "dva", "cva", "regime"),
  extended   = c("funding_spread", "dva", "cva", "regime", "collateral"),
  sides      = c("+", "-"),
  targets    = list(list(regime = "dynamic", mean_U = 1.39),
                    list(regime = "dynamic", mean_U = 0.3),
                    list(regime = "frozen1", mean_U = NA)),
  alphas     = c(0, 0.5, 1),
  ext_alphas = c(0.5, 1),           # alpha = 0 leaves the collateral feature empty
  rf_m       = c(0.075, 0.10),
  levels     = list(c(50L, 100L), c(100L, 200L)),
  seed       = 19L)

SHAPLEY_USED_BY <- paste("Ch5 Shapley table (headline decomposition);",
  "decomposition across alpha and r_f^- figure; ordering-dispersion table")

## The parameter object of a coalition: features not in 'on' are switched off.
shapley_params <- function(p_on, on) {
  p <- p_on
  if (!"funding_spread" %in% on) p <- with_params(p, rf_m = p$rf_p)
  if (!"dva" %in% on)            p <- with_params(p, LI = 0)
  if (!"cva" %in% on)            p <- with_params(p, LC = 0)
  if (!"regime" %in% on)         p <- with_params(p, beta_mode = "frozen0")
  if (!"collateral" %in% on)     p <- with_params(p, alpha = 0)
  p
}

## XVA_0 by Richardson on the two levels.
shapley_worth <- function(p, levels) {
  q <- vapply(levels, function(lv) quantisation_bsde(p, lv[1], lv[2])$U0, 0)
  2 * q[2] - q[1]
}

## All orderings of 1..k (k! of them), one per row.
all_orderings <- function(k) {
  if (k == 1L) return(matrix(1L, 1, 1))
  sub <- all_orderings(k - 1L)
  do.call(rbind, lapply(seq_len(k), function(first)
    cbind(first, matrix(setdiff(seq_len(k), first)[sub], nrow(sub)))))
}

# -----------------------------------------------------------------------------

## One decomposition: the target case p_on and a set of features.

# -----------------------------------------------------------------------------
shapley_case <- function(cs, p0, s) {
  feats <- if (cs$version == "main") s$features else s$extended
  if (cs$version == "main") feats <- setdiff(feats, "collateral")
  k <- length(feats)
  p_on <- with_params(p0, side = cs$side, beta_mode = cs$regime, alpha = cs$alpha,
                      rf_m = cs$rf_m)
  if (!is.na(cs$mean_U)) p_on <- with_params(p_on, mean_U = cs$mean_U)
  ## the main decomposition keeps collateral at the target's alpha throughout
  base_on <- if (cs$version == "main") "collateral" else character(0)

  ## worth of every coalition, indexed by bitmask
  masks <- 0:(2^k - 1)
  members <- function(mask) feats[bitwAnd(mask, 2^(seq_len(k) - 1)) > 0]
  worth <- vapply(masks, function(mk)
    shapley_worth(shapley_params(p_on, c(base_on, members(mk))), s$levels), 0)

  ## marginal contributions in every ordering
  ord <- all_orderings(k)
  marg <- t(apply(ord, 1, function(o) {
    mk <- 0; out <- numeric(k)
    for (i in o) { new <- mk + 2^(i - 1); out[i] <- worth[new + 1] - worth[mk + 1]; mk <- new }
    out
  }))
  phi <- colMeans(marg)

  ## Harsanyi dividends by Moebius inversion over the subsets of each mask
  size <- vapply(masks, function(mk) sum(bitwAnd(mk, 2^(seq_len(k) - 1)) > 0), 0)
  div  <- vapply(masks, function(mk) {
    sub <- masks[bitwAnd(masks, mk) == masks]
    sum((-1)^(size[mk + 1] - size[sub + 1]) * worth[sub + 1])
  }, 0)
  phi_div <- vapply(seq_len(k), function(i) {
    has <- bitwAnd(masks, 2^(i - 1)) > 0
    sum(div[has] / size[has])
  }, 0)

  Vhat0 <- bs_value(0, p0$S0, p0)
  total <- worth[2^k]; base <- worth[1]
  keys <- data.frame(case = cs$case, version = cs$version, side = cs$side,
                     regime = cs$regime, mean_U = cs$mean_U, alpha = cs$alpha,
                     rf_m = cs$rf_m, stringsAsFactors = FALSE)
  values <- cbind(keys[rep(1, k + 1), ], row.names = NULL,
    feature = c("baseline", feats),
    value = c(base, phi),
    pct_Vhat0 = 100 * c(base, phi) / Vhat0,
    pct_of_adjustment = c(NA, 100 * phi / (total - base)),
    min_marginal = c(NA, apply(marg, 2, min)),
    max_marginal = c(NA, apply(marg, 2, max)),
    XVA0 = total, XVA0_pct_Vhat0 = 100 * total / Vhat0,
    efficiency_residual = total - base - sum(phi),
    dividend_check = c(NA, phi - phi_div),
    orderings = nrow(ord),
    convention_flag = if (cs$side == "-") "buyer close-out convention" else "")
  coal <- cbind(keys[rep(1, length(masks)), ], row.names = NULL, mask = masks,
                on = vapply(masks, function(mk) paste(members(mk), collapse = "+"), ""),
                size = vapply(masks, function(mk) length(members(mk)), 0L),
                worth = worth, worth_pct_Vhat0 = 100 * worth / Vhat0)
  divs <- cbind(keys[rep(1, length(masks)), ], row.names = NULL, mask = masks,
                on = vapply(masks, function(mk) paste(members(mk), collapse = "+"), ""),
                size = size, dividend = div, dividend_pct_Vhat0 = 100 * div / Vhat0)
  orders <- cbind(keys[rep(1, nrow(ord) * k), ], row.names = NULL,
                  ordering = rep(seq_len(nrow(ord)), each = k),
                  sequence = rep(apply(ord, 1, function(o) paste(feats[o], collapse = ">")), each = k),
                  feature = rep(feats, nrow(ord)),
                  position = as.vector(t(apply(ord, 1, order))),
                  marginal = as.vector(t(marg)))
  list(values = values, coalitions = coal, orderings = orders, dividends = divs)
}

shapley_cases <- function(p0, s = SHAPLEY_SETTINGS) {
  one <- function(version, alphas) do.call(rbind, lapply(s$targets, function(tg)
    expand.grid(version = version, side = s$sides, regime = tg$regime, mean_U = tg$mean_U,
                alpha = alphas, rf_m = s$rf_m, stringsAsFactors = FALSE)))
  g <- rbind(one("main", s$alphas), one("extended", s$ext_alphas))
  g$case <- sprintf("%s %s %s%s alpha=%.2f rf-=%.3f", g$version,
                    ifelse(g$side == "+", "seller", "buyer"), g$regime,
                    ifelse(is.na(g$mean_U), "", sprintf(" lamU=%.2f", g$mean_U)),
                    g$alpha, g$rf_m)
  g
}

# -----------------------------------------------------------------------------

## run_shapley(): the suite.

# -----------------------------------------------------------------------------
run_shapley <- function(p0 = make_params(), s = SHAPLEY_SETTINGS, force = FALSE) {
  run_suite("shapley", "19_shapley.R", p0 = p0, settings = s, seed = s$seed,
            used_by = SHAPLEY_USED_BY, force = force, run = function() {
    for (lv in s$levels) invisible(get_nodes(p0, lv[1], lv[2]))   # before forking
    cases <- shapley_cases(p0, s)
    parts <- par_lapply(seq_len(nrow(cases)), function(i) shapley_case(cases[i, ], p0, s))
    pick  <- function(nm) do.call(rbind, lapply(parts, `[[`, nm))
    list(values = pick("values"), coalitions = pick("coalitions"),
         orderings = pick("orderings"), dividends = pick("dividends"))
  })
}

## Console report used by main.R: the headline decompositions.
report_shapley <- function(res) {
  v <- res$values
  for (ver in c("main", "extended")) {
    h <- v[v$version == ver & v$rf_m == 0.10 & v$alpha == 0.5, ]
    w <- reshape(h[, c("case", "feature", "pct_Vhat0")], idvar = "case",
                 timevar = "feature", direction = "wide")
    names(w) <- sub("^pct_Vhat0\\.", "", names(w))
    w$XVA0 <- h$XVA0_pct_Vhat0[match(w$case, h$case)]
    w$case <- sub("^(main|extended) ", "", w$case)
    cat(sprintf("\n%s decomposition (%d orderings), alpha = 0.5, rf- = 0.10, %% of Vhat_0\n",
                ver, h$orderings[1]))
    print(w, row.names = FALSE, digits = 3)
  }
  cat(sprintf("\nmax |efficiency residual| = %.1e; max |phi - dividend form| = %.1e\n",
              max(abs(v$efficiency_residual)), max(abs(v$dividend_check), na.rm = TRUE)))
  d <- res$dividends
  d <- d[d$version == "main" & d$rf_m == 0.10 & d$alpha == 0.5 & d$size > 0 &
         abs(d$dividend_pct_Vhat0) > 1e-6, ]
  cat("\nmain, alpha = 0.5, rf- = 0.10: non-zero Harsanyi dividends (% of Vhat_0)\n")
  print(data.frame(case = sub("^main ", "", d$case), coalition = d$on,
                   dividend = d$dividend_pct_Vhat0), row.names = FALSE, digits = 3)
  h <- v[v$version == "main" & v$rf_m == 0.10 & v$alpha == 0.5 & v$feature != "baseline" &
         v$side == "+", ]
  cat("\nseller, main: spread of the marginal contribution across the 24 orderings (% of Vhat_0)\n")
  p0 <- res$meta$params; V0 <- 100 / bs_value(0, p0$S0, p0)
  print(data.frame(case = sub("^main ", "", h$case), feature = h$feature,
                   phi = h$pct_Vhat0, min = h$min_marginal * V0, max = h$max_marginal * V0),
        row.names = FALSE, digits = 3)
  cat(sprintf("\nsuite time %.0f s on %d core(s)\n", res$meta$seconds, res$meta$cores))
}
