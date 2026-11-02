
##  00_parameters.R
##  Central parameter store for the replication of
##      Pang, W. and Sturm, S. (2020) "XVA Valuation under Market Illiquidity",
## -----------------------------------------------------------------------------

## make_params()
## Arguments override the benchmark values of Section 7.2 of the paper.

## -----------------------------------------------------------------------------
make_params <- function(...) {

  p <- list(

    ## ---------------- contract --------------------------------------------
    S0    = 1.00,      # initial stock price
    K     = 1.00,      # strike
    Tmat  = 0.25,      # maturity in years

    ## ---------------interest rates--------------------------------------- 
    rr_p  = 0.05,      # r_r^+  repo rate received when lending cash
    rr_m  = 0.05,      # r_r^-  repo rate paid when borrowing cash
    rc_p  = 0.01,      # r_c^+  collateral rate received by the provider
    rc_m  = 0.01,      # r_c^-  collateral rate paid by the receiver
    rf_p  = 0.05,      # r_f^+  funding rate received when lending to treasury
    rf_m  = 0.05,      # r_f^-  funding rate paid when borrowing from treasury
    rD    = 0.01,      # r_D    the investor's chosen discount rate

    ## --------------------credit-------------------------------
    muI   = 0.21,      # mu^I
    muC   = 0.16,      # mu^C
    hI    = 0.20,      # h_I^Q = mu^I - r_D
    hC    = 0.15,      # h_C^Q = mu^C - r_D
    LI    = 0.50,      # L_I loss rate on the investor's default
    LC    = 0.50,      # L_C loss rate on the counterparty's default

    ## ---- market -------------------------------------------------------------
    sigma = 0.30,      # stock volatility

    ## ---- collateral ---------------------------------------------------------
    ##C_t = alpha * Vhat_t with alpha in [0, 1].
    alpha = 0.00,

    ## ---- status process beta (Section 3, Table 2) ---------------------------
    ## Table 2 reports "Average Length (days) 507 / 361" together with
    ## "Estimates of lambda_U and lambda_V  1.39 / 0.99".  Since 507/365 = 1.389
    ## and 361/365 = 0.989, the reported figures are MEAN HOLDING TIMES IN YEARS
    ## and not exponential rates, even though Definition 1 calls them
    ## "parameters" of the exponential law.  This project therefore stores the
    ## means explicitly and derives the rates.  Set lambda_convention = "rate" to
    ## reinstate the literal reading of Definition 1.
    mean_U = 1.39,     # expected length of a calm regime, in years
    mean_V = 0.99,     # expected length of a crisis regime, in years
    lambda_convention = "mean",

    ## beta_mode selects which of the three experiments of Section 7.2 is run:
    ##   "frozen0"  beta_t == 0 for all t   (Figure 4a)
    ##   "frozen1"  beta_t == 1 for all t   (Figure 4b)
    ##   "dynamic"  beta is the alternating renewal process (Figure 5)
    beta_mode = "frozen0",

    ## ---- side of the trade --------------------------------------------------
    ## "+" is the seller's XVA of Definition 22, "-" is the buyer's.
    side = "+",

    ## ---- crisis convention for a short stock position (04_driver.R) ---------
    ## "banned": the control is projected onto the admissible set, Z >= 0 when
    ## beta = 1 (adopted for the thesis, docs/buyer_closeout.md); "extended":
    ## the generator as printed.  Only the buyer is affected.
    crisis_short = "banned"
  )

  ov <- list(...)
  unknown <- setdiff(names(ov), names(p))
  if (length(unknown) > 0L) {
    stop("make_params(): unknown parameter(s): ", paste(unknown, collapse = ", "))
  }
  p[names(ov)] <- ov
  derive_rates(p)
}

## Derived switching rates of the two-state Markov chain representation of
## beta.  a_U is the rate of leaving the calm state, a_V the rate of leaving
## the crisis state.
derive_rates <- function(p) {
  if (identical(p$lambda_convention, "mean")) {
    p$aU <- 1 / p$mean_U
    p$aV <- 1 / p$mean_V
  } else {
    p$aU <- p$mean_U
    p$aV <- p$mean_V
  }
  p
}

# -----------------------------------------------------------------------------

## with_params(): a copy of an existing parameter object with some entries
## changed and the switching rates re-derived.  Every suite starts from one
## parameter object and builds its cases with this, so a change to the base set
## propagates to every experiment.  The rates aU, aV are derived, not set.

# -----------------------------------------------------------------------------
with_params <- function(p, ...) {
  ov <- list(...)
  bad <- setdiff(names(ov), setdiff(names(p), c("aU", "aV")))
  if (length(bad) > 0L) {
    stop("with_params(): unknown or derived parameter(s): ", paste(bad, collapse = ", "))
  }
  p[names(ov)] <- ov
  derive_rates(p)
}

# -----------------------------------------------------------------------------

## check_assumptions()
##
## Verifies the standing hypotheses of the paper for a given parameter set.
## The function prints a report and returns TRUE when every condition holds.
## 

# -----------------------------------------------------------------------------
check_assumptions <- function(p, verbose = TRUE) {

  chk <- list()

  ## Assumption 12.1: r_f^+ <= r_f^-
  chk[["A12.1  rf+ <= rf-"]]          <- p$rf_p <= p$rf_m
  ## Assumption 12.2: max(r_f^+, r_D) < min(mu_I, mu_C)
  chk[["A12.2  max(rf+,rD) < min(muI,muC)"]] <- max(p$rf_p, p$rD) < min(p$muI, p$muC)
  ## Assumption 12.3: r_r^+ <= r_f^-  in the calm regime
  chk[["A12.3  rr+ <= rf-"]]          <- p$rr_p <= p$rf_m
  ## Proposition 14: r_r^+ <= r_f^+ <= r_r^-
  chk[["P14    rr+ <= rf+ <= rr-"]]   <- (p$rr_p <= p$rf_p) && (p$rf_p <= p$rr_m)

  ## Assumption 17(i)
  chk[["A17.i  rf- < 1/(5 T^{3/2})"]] <- p$rf_m < 1 / (5 * sqrt(p$Tmat^3))
  ## Assumption 17(ii)
  K17 <- (max(p$rf_m + p$rD, abs(p$rD - p$rr_m))) / min(p$sigma, 1)
  K17 <- max(K17, abs(p$rD - p$rr_p) / min(p$sigma, 1))
  chk[["A17.ii Lipschitz constant < 1/(5T)"]] <- K17 < 1 / (5 * p$Tmat)
  ## Assumption 17(iii).  The paper writes lambda_I, lambda_C for the default
  ## intensities that appear in the norm of Assumption 42; these are h_I^Q, h_C^Q.
  chk[["A17.iii rf- - rD < sqrt(h)/ (5T)"]] <-
    (p$rf_m - p$rD) < min(sqrt(p$hI), sqrt(p$hC)) / (5 * p$Tmat)

  if (verbose) {
    cat("Standing assumptions\n")
    for (nm in names(chk)) {
      cat(sprintf("  %-38s %s\n", nm, ifelse(chk[[nm]], "OK", "VIOLATED")))
    }
  }
  invisible(all(unlist(chk)))
}
