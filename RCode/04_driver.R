#-------------------------------------------------------------------------------
##  04_driver.R
##  The driver of the reduced XVA BSDE, assembled in the three layers used by
##  the paper.  Everything in this file is a pointwise, fully vectorised
##  function of numeric vectors of equal length.
##
##  LAYER 1 -- closeout adjustments (Equation (14))
##  LAYER 2 -- the XVA generator in the full filtration (Equation (15))
##  LAYER 3 -- the reduced generator (Equation (18))
##
##  The buyer's driver is obtained from the symmetry of Equation (15),
##      gbreve^-(t, u, z; beta, Vhat, Zhat) = -gbreve^+(t, -u, -z; beta, -Vhat, -Zhat).



# ----------------------- Layer 1 -------------------------------------------
pos_part <- function(x) pmax(x, 0)
neg_part <- function(x) pmax(-x, 0)

theta_tilde_I <- function(Vh, p) -p$LI * pos_part((1 - p$alpha) * Vh)
theta_tilde_C <- function(Vh, p)  p$LC * neg_part((1 - p$alpha) * Vh)

## ---------------Crisis convention for a short stock position -----------------
## How the equation treats a short stock position in the crisis (beta = 1),
## which Table 4.1 bans.  Documented in docs/buyer_closeout.md.
##   "extended" (p$crisis_short absent or "extended"): the generator is
##        defined for every Z; for beta = 1 and Z < 0 the indicator iota is zero,
##        so the short is financed by reverse repo as in the calm regime
##        (Section 4.3.1: the equation then solves an extended problem).
##   "banned" (the base parameter set, make_params()): the total control is
##        projected onto the admissible set, Z >= 0
##        when beta = 1; the hedger holds no stock in the crisis and earns no
##        reverse-repo carry.  The seller of a call never reaches the bound, so
##        only the buyer is affected.
crisis_banned  <- function(p) identical(p$crisis_short, "banned")
## (written with arithmetic, not ifelse(), so that a scalar beta with a vector of
## controls keeps the length of Z)
crisis_control <- function(Z, beta, p) if (crisis_banned(p)) Z + (beta == 1) * (pmax(Z, 0) - Z) else Z
## derivative of the projection with respect to the raw control
crisis_dZ <- function(Z, beta, p) if (crisis_banned(p)) as.numeric(!(beta == 1 & Z < 0)) else 1

## ---------------Layers 2 and 3, seller's side --------------------------------


g_plus <- function(u, z, beta, Vh, Zh, p, grad = FALSE) {

  Zraw <- z + Zh
  Z    <- crisis_control(Zraw, beta, p)
  iota <- as.numeric(Z > 0 & beta == 1)

  tI <- theta_tilde_I(Vh, p)
  tC <- theta_tilde_C(Vh, p)
  zI <- tI - u
  zC <- tC - u

  A <- u - (iota / p$sigma) * Z + zI + zC + (1 - p$alpha) * Vh

  bracket <-
      p$rf_p * pos_part(A) - p$rf_m * neg_part(A) +
      (p$rD - p$rr_m * (1 - iota)) / p$sigma * pos_part(Z) -
      (p$rD - p$rr_p * (1 - iota)) / p$sigma * neg_part(Z) +
      p$rc_p * pos_part(p$alpha * Vh) - p$rc_m * neg_part(p$alpha * Vh) -
      p$rD * zI - p$rD * zC

  g <- -bracket + p$rD * Vh + p$hI * (tI - u) + p$hC * (tC - u)

  if (!grad) return(g)

  rho <- p$rf_p * as.numeric(A > 0) + p$rf_m * as.numeric(A < 0)
  dgdu <- rho - 2 * p$rD - p$hI - p$hC
  dgdz <- -( rho * (-iota / p$sigma) +
             (p$rD - p$rr_m * (1 - iota)) / p$sigma * as.numeric(Z > 0) +
             (p$rD - p$rr_p * (1 - iota)) / p$sigma * as.numeric(Z < 0) ) *
          crisis_dZ(Zraw, beta, p)

  list(g = g, dgdu = dgdu, dgdz = dgdz)
}

# ---------------------- Buyer's side ------------------------------------------
## By the chain rule applied to gbreve^-(u,z) = -gbreve^+(-u,-z,-Vhat,-Zhat),
## the partial derivatives of gbreve^- equal those of gbreve^+ evaluated at the
## reflected arguments; the two sign changes cancel.
g_minus <- function(u, z, beta, Vh, Zh, p, grad = FALSE) {
  if (!grad) return(-g_plus(-u, -z, beta, -Vh, -Zh, p))
  r <- g_plus(-u, -z, beta, -Vh, -Zh, p, grad = TRUE)
  list(g = -r$g, dgdu = r$dgdu, dgdz = r$dgdz)
}

## Dispatch on p$side.
get_driver <- function(p) if (identical(p$side, "+")) g_plus else g_minus

# -----------------------------------------------------------------------------

## driver_terms(): the intermediate quantities of g_plus() -- the total control
## Z, the crisis indicator iota, the close-out adjustments, the jump controls
## Z^I = theta_I - u, Z^C = theta_C - u of (4.22) and the funding requirement A
## of (4.16) -- for the recovery of the replicating portfolio in 18.  The same
## expressions as g_plus(), seller's convention; the buyer's are obtained by
## the reflection (u, z, Vhat, Zhat) -> -(u, z, Vhat, Zhat).

## -----------------------------------------------------------------------------
driver_terms <- function(u, z, beta, Vh, Zh, p) {
  Z    <- crisis_control(z + Zh, beta, p)     # as in g_plus()
  iota <- as.numeric(Z > 0 & beta == 1)
  tI   <- theta_tilde_I(Vh, p)
  tC   <- theta_tilde_C(Vh, p)
  zI   <- tI - u
  zC   <- tC - u
  A    <- u - (iota / p$sigma) * Z + zI + zC + (1 - p$alpha) * Vh
  list(Z = Z, iota = iota, tI = tI, tC = tC, zI = zI, zC = zC, A = A)
}

# -----------------------------------------------------------------------------

## linear_coefficients()

# -----------------------------------------------------------------------------
linear_coefficients <- function(p, beta, rho_f) {
  iota <- if (beta == 1) 1 else 0
  list(iota = iota,
       c1 = rho_f - 2 * p$rD - p$hI - p$hC,
       c2 = (rho_f * iota - p$rD + p$rr_m * (1 - iota)) / p$sigma)
}

linear_c0 <- function(t, S, p, rho_f) {
  Vh <- bs_value(t, S, p)
  tI <- theta_tilde_I(Vh, p)
  tC <- theta_tilde_C(Vh, p)
  -rho_f * (tI + tC + (1 - p$alpha) * Vh) -
    p$rc_p * pos_part(p$alpha * Vh) + p$rc_m * neg_part(p$alpha * Vh) +
    (p$rD + p$hI) * tI + (p$rD + p$hC) * tC + p$rD * Vh
}

# -----------------------------------------------------------------------------

## Adjoint (reverse-mode) pass of the driver, for 14_quantisation_adjoint.R.
##
## Written by hand in the adjoint design paradigm of Capriotti and Giles:
## g_plus() is a straight-line sequence of assignments, and the adjoint visits
## them in reverse order, accumulating gbar * dg/d(input).  Given the incoming
## adjoint gbar (one value per node) it returns the adjoints of the node inputs
## u, z, Vh, Zh and, summed later over nodes, of the parameters sigma, r_f^-,
## h_I and h_C.  The kinks of the positive and negative parts are treated as in
## g_plus(grad = TRUE): the indicators are held locally constant.  As a check,
## the u-adjoint reproduces dgdu = rho - 2 r_D - h_I - h_C exactly.

# -----------------------------------------------------------------------------
g_plus_adjoint <- function(u, z, beta, Vh, Zh, p, gbar) {
  Zraw <- z + Zh
  Z    <- crisis_control(Zraw, beta, p)
  iota <- as.numeric(Z > 0 & beta == 1)
  a1   <- (1 - p$alpha) * Vh
  tI   <- theta_tilde_I(Vh, p)
  tC   <- theta_tilde_C(Vh, p)
  zI   <- tI - u
  zC   <- tC - u
  A    <- u - (iota / p$sigma) * Z + zI + zC + (1 - p$alpha) * Vh
  kp   <- (p$rD - p$rr_m * (1 - iota)) / p$sigma        # coefficient of Z^+
  km   <- (p$rD - p$rr_p * (1 - iota)) / p$sigma        # coefficient of Z^-
  rho  <- p$rf_p * as.numeric(A > 0) + p$rf_m * as.numeric(A < 0)

  ## g = -bracket + r_D Vh + h_I (tI - u) + h_C (tC - u)
  bb  <- -gbar
  Vhb <- p$rD * gbar
  hIb <- gbar * (tI - u)
  hCb <- gbar * (tC - u)
  tIb <- p$hI * gbar
  tCb <- p$hC * gbar
  ub  <- -(p$hI + p$hC) * gbar

  ## bracket
  Ab   <- bb * rho
  rfmb <- bb * (-neg_part(A))
  Zb   <- bb * (kp * as.numeric(Z > 0) + km * as.numeric(Z < 0))
  sgb  <- bb * (-(kp * pos_part(Z) - km * neg_part(Z)) / p$sigma)
  Vhb  <- Vhb + bb * p$alpha * (p$rc_p * as.numeric(p$alpha * Vh > 0) +
                                p$rc_m * as.numeric(p$alpha * Vh < 0))
  zIb  <- -p$rD * bb
  zCb  <- -p$rD * bb
  ## collateral interest r_c^+ (alpha Vh)^+ - r_c^- (alpha Vh)^- in alpha
  alb  <- bb * Vh * (p$rc_p * as.numeric(p$alpha * Vh > 0) +
                     p$rc_m * as.numeric(p$alpha * Vh < 0))

  ## A = u - (iota/sigma) Z + zI + zC + (1 - alpha) Vh
  ub  <- ub + Ab
  Zb  <- Zb - (iota / p$sigma) * Ab
  sgb <- sgb + (iota / p$sigma^2) * Z * Ab
  zIb <- zIb + Ab
  zCb <- zCb + Ab
  Vhb <- Vhb + (1 - p$alpha) * Ab
  alb <- alb - Vh * Ab

  ## zI = tI - u, zC = tC - u
  tIb <- tIb + zIb; ub <- ub - zIb
  tCb <- tCb + zCb; ub <- ub - zCb

  ## tI = -L_I (a1)^+,  tC = L_C (a1)^-,  a1 = (1 - alpha) Vh
  Vhb <- Vhb + tIb * (-p$LI * (1 - p$alpha) * as.numeric(a1 > 0)) +
               tCb * (-p$LC * (1 - p$alpha) * as.numeric(a1 < 0))
  ## ... and in the loss rates and in alpha (through a1 = (1 - alpha) Vh)
  LIb <- tIb * (-pos_part(a1))
  LCb <- tCb * neg_part(a1)
  a1b <- tIb * (-p$LI * as.numeric(a1 > 0)) + tCb * (-p$LC * as.numeric(a1 < 0))
  alb <- alb - Vh * a1b

  ## Z = crisis_control(z + Zh): the projection passes the adjoint through
  ## where it is inactive and stops it where it binds
  Zb <- Zb * crisis_dZ(Zraw, beta, p)
  list(u = ub, z = Zb, Vh = Vhb, Zh = Zb,
       sigma = sgb, rf_m = rfmb, hI = hIb, hC = hCb,
       alpha = alb, LI = LIb, LC = LCb)
}

## Buyer: g^-(u,z;Vh,Zh) = -g^+(-u,-z;-Vh,-Zh).  The adjoint of g^+ is taken
## at the reflected point with adjoint -gbar; the state adjoints change sign
## again on the way back, the parameter adjoints do not.
g_minus_adjoint <- function(u, z, beta, Vh, Zh, p, gbar) {
  a <- g_plus_adjoint(-u, -z, beta, -Vh, -Zh, p, -gbar)
  list(u = -a$u, z = -a$z, Vh = -a$Vh, Zh = -a$Zh,
       sigma = a$sigma, rf_m = a$rf_m, hI = a$hI, hC = a$hC,
       alpha = a$alpha, LI = a$LI, LC = a$LC)
}

get_driver_adjoint <- function(p) {
  if (identical(p$side, "+")) g_plus_adjoint else g_minus_adjoint
}
