#-------------------------------------------------------------------------------

##  18_replication.R
##  Suite (b): the replicating portfolio.
##
##  PART 1 -- RECOVERY.  From the solver's (Ybreve, Zbreve) the six accounts of
##  Table 4.1 are recovered in units and in currency, by (4.12), (4.18), (4.9)
##  and (4.10):
##      total control       Z   = Zbreve + Zhat
##      jump controls       Z^I = theta_I(Vhat) - Ybreve,  Z^C = theta_C(Vhat) - Ybreve
##      stock               xi^S = Z / (sigma S),      value (1 - beta 1{xi^S < 0}) xi^S S
##      investor bond       xi^I = -Z^I / P^I,          value -Z^I
##      counterparty bond   xi^C = -Z^C / P^C,          value -Z^C
##      repo                value -(1 - beta) xi^S S    (4.9)
##      collateral          value C = alpha Vhat
##      funding             value A of (4.16), the balancing item (4.10)
##  with P^i_t = exp(-(r_D + h_i)(T - t)) the pre-default bond price.  The
##  accounts must sum to the wealth V = Vhat + Ybreve (the balance check; it is
##  an identity of the model, so it checks that the recovery and the driver's
##  funding requirement agree).  A node with beta = 1 and Z < 0 requires a short
##  sale in a crisis and is flagged INADMISSIBLE (Section 4.3.1).
##
##  Units of the three cash accounts are the currency value divided by the
##  account's numeraire B_t = exp(r t), with r the rate for the sign of the
##  position; this is exact when the sign has not changed on [0, t].
##
##  THE BUYER.  The buyer's hedge portfolio has wealth -V^- and is obtained by
##  the reflection of (4.17): the seller's recovery applied to
##  (-Ybreve^-, -Zbreve^-, -Vhat, -Zhat).  Its stock position is short, and in
##  a crisis therefore inadmissible.  This rests on the buyer's close-out
##  convention of the driver, which is an open methodological decision; the
##  buyer's rows are reported but flagged.
##
##  PART 2 -- THE SELF-FINANCING SIMULATOR (seller only).  The discretised
##  strategy is run along exact paths of (S, beta) and of the two default times
##  under Q, on a fine grid of nt steps.  At each rebalancing date the stock,
##  bond, repo and collateral accounts are set to the model's targets at the
##  current state and the funding account takes the difference between the
##  actual wealth W and their sum.  Between dates the units are held and each
##  cash account accrues at the rate for its sign; when beta switches into the
##  crisis, the repo balance is moved into the funding account (wealth neutral:
##  the repo market freezes, the treasury desk takes over the financing).  A
##  short stock position is worth nothing in the crisis (Table 4.1).  At the
##  first default, the replication error is W - theta_i(Vhat); otherwise it is
##  W_T - (S_T - K)^+.  The controls come either from the PDE solution (any
##  rebalancing step, which isolates the error of discrete rebalancing) or from
##  the quantisation tree (rebalancing on its own dates).  All runs share the
##  same random numbers, so differences across steps are not sampling noise.
##
##  OUTPUT (output/data/):
##    replication_t0.csv             the six accounts at t = 0, all cases
##    replication_nodes.csv          the accounts at every tree node (headline)
##    replication_admissibility.csv  inadmissible nodes and balance residuals
##    replication_sim.csv            replication error by step and control source
##    replication_decay.csv          fitted decay rate of the error in the step
##
##  USAGE.  Rscript main.R replicate [force]

#-------------------------------------------------------------------------------

REPLICATION_SETTINGS <- list(
  ## part 1
  sides   = c("+", "-"),
  regimes = c("frozen0", "frozen1", "dynamic"),
  alphas  = c(0, 0.5, 1),
  rf_m    = c(0.05, 0.075, 0.10),
  t0_levels = list(c(50L, 100L), c(100L, 200L)),  # t = 0 values by Richardson
  node_level = c(50L, 100L),                       # tree for the node tables
  node_csv  = list(alpha = 0.5, rf_m = 0.10),      # headline slice written as CSV
  ## part 2
  sim_regimes = list(list(regime = "frozen0", mean_U = NA),
                     list(regime = "frozen1", mean_U = NA),
                     list(regime = "dynamic", mean_U = 1.39),
                     list(regime = "dynamic", mean_U = 0.3)),
  sim_alpha = 0.5, sim_rf_m = 0.10,
  sim_M     = 20000L,                              # paths
  sim_nt    = 1600L,                               # fine steps (and PDE steps)
  sim_pde_nx = 801L,
  sim_m     = c(25L, 50L, 100L, 200L, 400L, 800L, 1600L),  # rebalancing dates
  sim_quant = list(c(25L, 50L), c(50L, 100L), c(100L, 200L)),
  seed      = 18L)

REPLICATION_USED_BY <- paste("Ch5 replicating-portfolio table (t = 0 accounts);",
  "positions over (t, S) figure; admissibility table; replication-error decay figure")

ACCOUNTS <- c("stock", "funding", "investor_bond", "counterparty_bond", "repo", "collateral")

# -----------------------------------------------------------------------------

## portfolio_accounts(): the six accounts at nodes (vectors of equal length).
##   t, S, beta   time, stock price, regime
##   y, z         the solver's Ybreve, Zbreve for the side p$side
## Returns a data frame with units and values of every account, the wealth, the
## balance residual and the admissibility flag.

# -----------------------------------------------------------------------------
portfolio_accounts <- function(t, S, beta, y, z, p) {
  s   <- if (identical(p$side, "+")) 1 else -1        # buyer: reflection (4.17)
  Vh  <- s * bs_value(t, S, p)
  Zh  <- s * bs_control(t, S, p)
  d   <- driver_terms(s * y, s * z, beta, Vh, Zh, p)
  V   <- Vh + s * y                                   # wealth of the hedge portfolio
  xiS <- d$Z / (p$sigma * S)
  PI  <- exp(-(p$rD + p$hI) * (p$Tmat - t))
  PC  <- exp(-(p$rD + p$hC) * (p$Tmat - t))
  val <- data.frame(
    stock             = (1 - beta * (xiS < 0)) * xiS * S,
    funding           = d$A,
    investor_bond     = -d$zI,
    counterparty_bond = -d$zC,
    repo              = -(1 - beta) * xiS * S,
    collateral        = p$alpha * Vh)
  rate <- function(x, rp, rm) ifelse(x >= 0, rp, rm)
  units <- data.frame(
    stock             = xiS,
    funding           = val$funding / exp(rate(val$funding, p$rf_p, p$rf_m) * t),
    investor_bond     = -d$zI / PI,
    counterparty_bond = -d$zC / PC,
    repo              = val$repo / exp(rate(val$repo, p$rr_p, p$rr_m) * t),
    ## Table 4.1: the collateral account holds psi^c with -psi^c B^{r_c} = C
    collateral        = -val$collateral / exp(rate(val$collateral, p$rc_p, p$rc_m) * t))
  names(val)   <- paste0("value_", names(val))
  names(units) <- paste0("units_", names(units))
  data.frame(t = t, S = S, beta = beta, Vhat = Vh, Ybreve = s * y, wealth = V,
             Z_total = d$Z, ZI = d$zI, ZC = d$zC, units, val,
             residual = rowSums(val) - V,
             admissible = !(beta == 1 & d$Z < 0))
}

## The cases of part 1, built from the base parameter object.
replication_cases <- function(p0, s = REPLICATION_SETTINGS) {
  g <- expand.grid(side = s$sides, regime = s$regimes, alpha = s$alphas,
                   rf_m = s$rf_m, stringsAsFactors = FALSE)
  g$mean_U <- ifelse(g$regime == "dynamic", p0$mean_U, NA_real_)
  g$case <- sprintf("%s %s alpha=%.2f rf-=%.3f",
                    ifelse(g$side == "+", "seller", "buyer"), g$regime, g$alpha, g$rf_m)
  g
}

# -----------------------------------------------------------------------------

## Part 1 for one case: the t = 0 accounts (Richardson-extrapolated Y0, Z0) and
## the accounts at every node of the tree.

# -----------------------------------------------------------------------------
replication_case <- function(cs, p0, s) {
  p  <- with_params(p0, side = cs$side, beta_mode = cs$regime, alpha = cs$alpha,
                    rf_m = cs$rf_m)
  b0 <- if (identical(cs$regime, "frozen1")) 1 else 0
  ## t = 0: Richardson on the two levels (first order, Module 17)
  q  <- lapply(s$t0_levels, function(lv) quantisation_bsde(p, lv[1], lv[2]))
  Y0 <- 2 * q[[2]]$U0 - q[[1]]$U0
  Z0 <- 2 * q[[2]]$Z0 - q[[1]]$Z0
  a0 <- portfolio_accounts(0, p$S0, b0, Y0, Z0, p)
  Vhat0 <- bs_value(0, p$S0, p)
  t0 <- data.frame(case = cs$case, side = cs$side, regime = cs$regime, alpha = cs$alpha,
                   rf_m = cs$rf_m, mean_U = cs$mean_U, account = ACCOUNTS,
                   units = unlist(a0[paste0("units_", ACCOUNTS)]),
                   value = unlist(a0[paste0("value_", ACCOUNTS)]),
                   stringsAsFactors = FALSE, row.names = NULL)
  t0$value_pct_Vhat0 <- 100 * t0$value / Vhat0
  t0$wealth   <- a0$wealth
  t0$XVA0     <- Y0
  t0$XVA0_pct_Vhat0 <- 100 * Y0 / Vhat0
  t0$residual <- a0$residual
  t0$admissible <- a0$admissible

  ## every node of the tree, every regime column carried
  n <- s$node_level[1]; N <- s$node_level[2]
  qk <- quantisation_bsde(p, n, N, keep = TRUE)
  regs <- qk$regimes
  nodes <- do.call(rbind, lapply(seq_len(n), function(k) {
    x <- qk$grids[[k]]; pk <- qk$probs[[k]]
    do.call(rbind, lapply(seq_along(regs), function(j) {
      a <- portfolio_accounts((k - 1) * p$Tmat / n, exp(x), regs[j],
                              qk$U[[k]][, j], qk$Z[[k]][, j], p)
      cbind(case = cs$case, side = cs$side, regime = cs$regime, alpha = cs$alpha,
            rf_m = cs$rf_m, k = k - 1L, node = seq_along(x), weight = pk, a,
            stringsAsFactors = FALSE)
    }))
  }))
  crisis <- nodes$beta == 1
  adm <- data.frame(case = cs$case, side = cs$side, regime = cs$regime, alpha = cs$alpha,
                    rf_m = cs$rf_m, n = n, N = N, nodes = nrow(nodes),
                    crisis_nodes = sum(crisis),
                    inadmissible_nodes = sum(!nodes$admissible),
                    ## tree probability of an inadmissible crisis node, averaged
                    ## over dates (conditional on the crisis column)
                    inadmissible_mass = if (any(crisis))
                      sum(nodes$weight[!nodes$admissible]) / n else 0,
                    min_crisis_Z = if (any(crisis)) min(nodes$Z_total[crisis]) else NA_real_,
                    max_abs_residual = max(abs(nodes$residual)),
                    t0_admissible = a0$admissible,
                    convention_flag = if (cs$side == "-")
                      "buyer: reflected strategy, depends on the buyer close-out convention" else "",
                    stringsAsFactors = FALSE)
  list(t0 = t0, nodes = nodes, adm = adm)
}

# -----------------------------------------------------------------------------

## Part 2: the self-financing simulator.
##
## ctrl(j, S, beta) returns list(Y, Zb) at fine step j for vectors S, beta.

# -----------------------------------------------------------------------------
pde_controls <- function(p, s) {
  op  <- pde_operator(p, s$sim_pde_nx, s$sim_nt)
  sol <- pde_solve(p, op, keep = TRUE)
  x   <- sol$x; dx <- x[2] - x[1]
  ## the solution in regime r at fine step j, for the regime-switch diagnostic
  value <- function(j, S, r) {
    r <- rep_len(r, length(S)); out <- numeric(length(S))
    for (rr in unique(r)) {
      i <- r == rr
      out[i] <- approx(x, sol$hist[[as.character(rr)]][, j + 1L], log(S[i]), rule = 2)$y
    }
    out
  }
  structure(function(j, S, beta) {
    Y <- Zb <- numeric(length(S)); lx <- log(S)
    for (r in unique(beta)) {
      i <- beta == r; u <- sol$hist[[as.character(r)]][, j + 1L]
      Y[i]  <- approx(x, u, lx[i], rule = 2)$y
      Zb[i] <- p$sigma * approx(x, d_dx(u, dx), lx[i], rule = 2)$y
    }
    list(Y = Y, Zb = Zb)
  }, value = if (identical(p$beta_mode, "dynamic")) value else NULL)
}

quant_controls <- function(p, n, N, nt) {
  q <- quantisation_bsde(p, n, N, keep = TRUE)
  every <- nt %/% n
  col <- function(r) match(r, q$regimes)
  function(j, S, beta) {
    k <- j %/% every                                # the tree date
    x <- q$grids[[k + 1L]]; lx <- log(S)
    Y <- Zb <- numeric(length(S))
    for (r in unique(beta)) {
      i <- beta == r; cl <- col(r)
      if (length(x) == 1L) { Y[i] <- q$U[[k + 1L]][1, cl]; Zb[i] <- q$Z[[k + 1L]][1, cl]; next }
      Y[i]  <- approx(x, q$U[[k + 1L]][, cl], lx[i], rule = 2)$y
      Zb[i] <- approx(x, q$Z[[k + 1L]][, cl], lx[i], rule = 2)$y
    }
    list(Y = Y, Zb = Zb)
  }
}

simulate_replication <- function(p, ctrl, m, s) {
  stopifnot(identical(p$side, "+"), s$sim_nt %% m == 0L)
  set.seed(s$seed)                                  # the same paths for every run
  M <- s$sim_M; nt <- s$sim_nt; ds <- p$Tmat / nt; every <- nt %/% m
  dyn <- identical(p$beta_mode, "dynamic")
  tauI <- rexp(M, p$hI); tauC <- rexp(M, p$hC)
  if (dyn) tr <- beta_transition(p, ds)
  S    <- rep(p$S0, M)
  beta <- rep(if (identical(p$beta_mode, "frozen1")) 1 else 0, M)
  rate <- function(x, rp, rm) ifelse(x >= 0, rp, rm)
  alive <- rep(TRUE, M); err <- rep(NA_real_, M); event <- rep("maturity", M)
  inadm <- 0L
  xiS <- xiI <- xiC <- R <- Cl <- F <- numeric(M)
  ## Regime-switch martingale (diagnostic, PDE controls, dynamic regime):
  ##   Jm = sum over switches of [u(t, S, beta_new) - u(t, S, beta_old)]
  ##        - int a(beta) [u(t, S, 1 - beta) - u(t, S, beta)] dt,
  ## the component of the XVA that the reduced equation (4.23) leaves out and
  ## that no account of the portfolio can hedge.  The error is predicted to be
  ## -Jm plus the error of discrete rebalancing.
  uval <- attr(ctrl, "value"); track <- dyn && !is.null(uval)
  Jm <- numeric(M)

  c0 <- ctrl(0L, S[1], beta[1])
  W  <- rep(bs_value(0, p$S0, p) + c0$Y, M)          # initial wealth V_0 = Vhat_0 + XVA_0

  for (j in 0:(nt - 1L)) {
    t <- j * ds
    if (j %% every == 0L) {                         # rebalance the living paths
      i  <- which(alive)
      cc <- ctrl(j, S[i], beta[i])
      Vh <- bs_value(t, S[i], p); Zh <- bs_control(t, S[i], p)
      d  <- driver_terms(cc$Y, cc$Zb, beta[i], Vh, Zh, p)
      x  <- d$Z / (p$sigma * S[i])
      short_crisis <- beta[i] == 1 & x < 0
      inadm <- inadm + sum(short_crisis)
      x[short_crisis] <- 0                          # no short sale in a crisis
      PI <- exp(-(p$rD + p$hI) * (p$Tmat - t)); PC <- exp(-(p$rD + p$hC) * (p$Tmat - t))
      xiS[i] <- x
      xiI[i] <- -d$zI / PI; xiC[i] <- -d$zC / PC
      R[i]   <- -(1 - beta[i]) * x * S[i]
      Cl[i]  <- p$alpha * Vh
      F[i]   <- W[i] - x * S[i] - xiI[i] * PI - xiC[i] * PC - R[i] - Cl[i]
    }
    ## one fine step
    S  <- S * exp((p$rD - 0.5 * p$sigma^2) * ds + p$sigma * sqrt(ds) * rnorm(M))
    if (track) {
      i  <- which(alive)
      u0 <- uval(j, S[i], 0); u1 <- uval(j, S[i], 1)
      a  <- ifelse(beta[i] == 0, p$aU, p$aV)
      Jm[i] <- Jm[i] - a * ifelse(beta[i] == 0, u1 - u0, u0 - u1) * ds
      beta_old <- beta
    }
    if (dyn) {
      u  <- runif(M)
      sw <- u < ifelse(beta == 0, tr$p01, tr$p10)
      into_crisis <- sw & beta == 0
      F[into_crisis] <- F[into_crisis] + R[into_crisis]; R[into_crisis] <- 0
      beta <- ifelse(sw, 1 - beta, beta)
    }
    if (track) {
      i <- which(alive & beta != beta_old)
      if (length(i)) Jm[i] <- Jm[i] + uval(j + 1L, S[i], beta[i]) - uval(j + 1L, S[i], beta_old[i])
    }
    R  <- R  * exp(rate(R,  p$rr_p, p$rr_m) * ds)
    Cl <- Cl * exp(rate(Cl, p$rc_p, p$rc_m) * ds)
    F  <- F  * exp(rate(F,  p$rf_p, p$rf_m) * ds)
    t1 <- t + ds
    PI <- exp(-(p$rD + p$hI) * (p$Tmat - t1)); PC <- exp(-(p$rD + p$hC) * (p$Tmat - t1))
    stock <- (1 - beta * (xiS < 0)) * xiS * S
    W <- stock + xiI * PI + xiC * PC + R + Cl + F

    ## defaults in (t, t1]: the defaulted party's bond is worth zero
    dI <- alive & tauI <= t1 & tauI < tauC
    dC <- alive & tauC <= t1 & tauC < tauI
    if (any(dI)) {
      Vh <- bs_value(t1, S[dI], p)
      err[dI] <- (W[dI] - xiI[dI] * PI) - (Vh - p$LI * pmax((1 - p$alpha) * Vh, 0))
      event[dI] <- "investor_default"; alive[dI] <- FALSE
    }
    if (any(dC)) {
      Vh <- bs_value(t1, S[dC], p)
      err[dC] <- (W[dC] - xiC[dC] * PC) - (Vh + p$LC * pmax(-(1 - p$alpha) * Vh, 0))
      event[dC] <- "counterparty_default"; alive[dC] <- FALSE
    }
  }
  err[alive] <- W[alive] - pmax(S[alive] - p$K, 0)
  list(err = err, event = event, inadmissible_rebalances = inadm,
       err_adj = if (track) err + Jm else rep(NA_real_, M))
}

replication_sim_case <- function(sc, p0, s) {
  p <- with_params(p0, side = "+", beta_mode = sc$regime, alpha = s$sim_alpha,
                   rf_m = s$sim_rf_m)
  if (!is.na(sc$mean_U)) p <- with_params(p, mean_U = sc$mean_U)
  Vhat0 <- bs_value(0, p$S0, p)
  lab <- sprintf("seller %s%s alpha=%.2f rf-=%.3f", sc$regime,
                 if (is.na(sc$mean_U)) "" else sprintf(" lamU=%.2f", sc$mean_U),
                 s$sim_alpha, s$sim_rf_m)
  runs <- c(lapply(s$sim_m, function(m) list(source = "pde", m = m, n = NA, N = NA)),
            lapply(s$sim_quant, function(l) list(source = "quantisation", m = l[1],
                                                 n = l[1], N = l[2])))
  pc <- pde_controls(p, s)
  do.call(rbind, lapply(runs, function(r) {
    ctrl <- if (r$source == "pde") pc else quant_controls(p, r$n, r$N, s$sim_nt)
    t0  <- proc.time()[["elapsed"]]
    out <- simulate_replication(p, ctrl, r$m, s)
    e   <- out$err
    ev  <- function(k) if (any(out$event == k)) mean(e[out$event == k]) else NA_real_
    data.frame(case = lab, regime = sc$regime, mean_U = sc$mean_U, alpha = s$sim_alpha,
               rf_m = s$sim_rf_m, source = r$source, n = r$n, N = r$N,
               rebalances = r$m, step = p$Tmat / r$m, paths = length(e),
               mean = mean(e), sd = sd(e), rmse = sqrt(mean(e^2)),
               mean_abs = mean(abs(e)),
               q01 = quantile(e, 0.01, names = FALSE), q99 = quantile(e, 0.99, names = FALSE),
               rmse_pct_Vhat0 = 100 * sqrt(mean(e^2)) / Vhat0,
               mean_pct_Vhat0 = 100 * mean(e) / Vhat0,
               se_mean = sd(e) / sqrt(length(e)),
               frac_investor_default = mean(out$event == "investor_default"),
               frac_counterparty_default = mean(out$event == "counterparty_default"),
               mean_err_maturity = ev("maturity"),
               mean_err_investor_default = ev("investor_default"),
               mean_err_counterparty_default = ev("counterparty_default"),
               inadmissible_rebalances = out$inadmissible_rebalances,
               ## error with the unhedgeable regime-switch martingale added back
               rmse_adj_pct_Vhat0 = 100 * sqrt(mean(out$err_adj^2)) / Vhat0,
               rmse_regime_pct_Vhat0 = 100 * sqrt(mean((out$err_adj - e)^2)) / Vhat0,
               seconds = proc.time()[["elapsed"]] - t0, stringsAsFactors = FALSE)
  }))
}

## Least-squares slope of log RMSE (and log |mean|) on log step, PDE controls.
replication_decay <- function(sim) {
  d <- sim[sim$source == "pde", ]
  do.call(rbind, lapply(split(d, d$case), function(x) {
    fr <- lm(log(rmse) ~ log(step), data = x)
    fm <- lm(log(abs(mean)) ~ log(step), data = x)
    adj <- if (all(is.finite(x$rmse_adj_pct_Vhat0)))
      unname(coef(lm(log(rmse_adj_pct_Vhat0) ~ log(step), data = x))[2]) else NA_real_
    data.frame(case = x$case[1], rmse_slope = unname(coef(fr)[2]),
               rmse_slope_se = summary(fr)$coefficients[2, 2],
               rmse_adj_slope = adj,
               mean_slope = unname(coef(fm)[2]),
               rmse_at_finest_pct_Vhat0 = x$rmse_pct_Vhat0[which.min(x$step)],
               stringsAsFactors = FALSE)
  }))
}

# -----------------------------------------------------------------------------

## run_replication(): the suite.

# -----------------------------------------------------------------------------
run_replication <- function(p0 = make_params(), s = REPLICATION_SETTINGS, force = FALSE) {
  run_suite("replication", "18_replication.R", p0 = p0, settings = s,
            seed = s$seed, used_by = REPLICATION_USED_BY, force = force, run = function() {
    for (lv in c(s$t0_levels, list(s$node_level), s$sim_quant)) {
      invisible(get_nodes(p0, lv[1], lv[2]))         # trees built once, before forking
    }
    cases <- replication_cases(p0, s)
    parts <- par_lapply(seq_len(nrow(cases)), function(i) replication_case(cases[i, ], p0, s))
    nodes <- do.call(rbind, lapply(parts, `[[`, "nodes"))
    keep  <- nodes$alpha == s$node_csv$alpha & nodes$rf_m == s$node_csv$rf_m
    sim   <- do.call(rbind, par_lapply(s$sim_regimes, replication_sim_case, p0 = p0, s = s))
    list(t0 = do.call(rbind, lapply(parts, `[[`, "t0")),
         nodes = nodes[keep, ],
         admissibility = do.call(rbind, lapply(parts, `[[`, "adm")),
         sim = sim, decay = replication_decay(sim))
  })
}

## Console report used by main.R.
report_replication <- function(res) {
  t0 <- res$t0
  h  <- t0[t0$alpha == 0.5 & t0$rf_m == 0.10, ]
  w  <- reshape(h[, c("case", "account", "value")], idvar = "case", timevar = "account",
                direction = "wide")
  names(w) <- sub("^value\\.", "", names(w))
  w$XVA0  <- h$XVA0[match(w$case, h$case)]
  w$resid <- h$residual[match(w$case, h$case)]
  cat("\nthe six accounts at t = 0 (currency), alpha = 0.5, rf- = 0.10\n")
  print(w, row.names = FALSE, digits = 4)
  cat("\nadmissibility over the tree\n")
  a <- res$admissibility
  print(aggregate(cbind(inadmissible_nodes, crisis_nodes, max_abs_residual) ~ side + regime,
                  a, max), row.names = FALSE, digits = 3)
  cat("\nreplication error (seller), PDE controls: RMSE in % of Vhat_0 by rebalancing count\n")
  sm <- res$sim
  print(reshape(sm[sm$source == "pde", c("case", "rebalances", "rmse_pct_Vhat0")],
                idvar = "case", timevar = "rebalances", direction = "wide"),
        row.names = FALSE, digits = 3)
  dd <- sm[sm$source == "pde" & sm$regime == "dynamic", ]
  if (nrow(dd)) {
    cat("\ndynamic regime: RMSE after adding back the regime-switch martingale (% of Vhat_0)\n")
    print(reshape(dd[, c("case", "rebalances", "rmse_adj_pct_Vhat0")],
                  idvar = "case", timevar = "rebalances", direction = "wide"),
          row.names = FALSE, digits = 3)
    cat("size of the regime-switch martingale alone (% of Vhat_0):",
        paste(sprintf("%s %.2f", dd$case[dd$rebalances == max(dd$rebalances)],
                      dd$rmse_regime_pct_Vhat0[dd$rebalances == max(dd$rebalances)]),
              collapse = "; "), "\n")
  }
  cat("\nquantisation controls (rebalancing on the tree dates)\n")
  print(sm[sm$source == "quantisation", c("case", "n", "N", "rmse_pct_Vhat0",
                                         "mean_pct_Vhat0", "inadmissible_rebalances")],
        row.names = FALSE, digits = 3)
  cat("\ndecay of the RMSE in the rebalancing step (slope of log RMSE on log step)\n")
  print(res$decay, row.names = FALSE, digits = 3)
  cat(sprintf("\nsuite time %.0f s on %d core(s)\n", res$meta$seconds, res$meta$cores))
}
