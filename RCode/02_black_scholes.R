
##  02_black_scholes.R
##  The third-party valuation (Vhat, Zhat) of Equation (12).
#-------------------------------------------------------------------------------

## d1 and d2, guarding the degenerate case t = T.
bs_d1d2 <- function(t, S, p) {
  tau <- pmax(p$Tmat - t, 0)
  if (length(tau) == 1L) tau <- rep(tau, length(S))
  d1  <- rep(NA_real_, length(S))
  d2  <- rep(NA_real_, length(S))
  pos <- tau > 0
  if (any(pos)) {
    tt <- tau[pos]
    Sp <- S[pos]
    d1[pos] <- (log(Sp / p$K) + (p$rD + 0.5 * p$sigma^2) * tt) / (p$sigma * sqrt(tt))
    d2[pos] <- d1[pos] - p$sigma * sqrt(tt)
  }
  list(d1 = d1, d2 = d2, tau = tau)
}

## Vhat(t, S): the third-party (Black-Scholes) value of the European call.
bs_value <- function(t, S, p) {
  dd  <- bs_d1d2(t, S, p)
  out <- pmax(S - p$K, 0)                       # terminal payoff, used when tau = 0
  pos <- dd$tau > 0
  if (any(pos)) {
    out[pos] <- S[pos] * pnorm(dd$d1[pos]) -
      p$K * exp(-p$rD * dd$tau[pos]) * pnorm(dd$d2[pos])
  }
  out
}

## Zhat(t, S) = sigma S N(d1): the martingale representation coefficient of Vhat.
bs_control <- function(t, S, p) {
  dd  <- bs_d1d2(t, S, p)
  out <- p$sigma * S * as.numeric(S > p$K)      # limit as t -> T
  pos <- dd$tau > 0
  if (any(pos)) out[pos] <- p$sigma * S[pos] * pnorm(dd$d1[pos])
  out
}

## Delta, provided only for reporting and for the sanity checks.
bs_delta <- function(t, S, p) {
  dd  <- bs_d1d2(t, S, p)
  out <- as.numeric(S > p$K)
  pos <- dd$tau > 0
  if (any(pos)) out[pos] <- pnorm(dd$d1[pos])
  out
}

# -----------------------------------------------------------------------------

## bs_partials(): first derivatives of the reference pair with respect to the
## stock and the volatility, used by the adjoint of 14_quantisation_adjoint.R.

# -----------------------------------------------------------------------------
bs_partials <- function(t, S, p) {
  dd <- bs_d1d2(t, S, p)
  if (any(dd$tau <= 0)) stop("bs_partials(): defined for t < T only.")
  st <- sqrt(dd$tau)
  N1 <- pnorm(dd$d1); f1 <- dnorm(dd$d1)
  list(dV_dS = N1, dV_dsigma = S * f1 * st,
       dZ_dS = p$sigma * N1 + f1 / st,
       dZ_dsigma = S * N1 - S * f1 * dd$d2)
}
