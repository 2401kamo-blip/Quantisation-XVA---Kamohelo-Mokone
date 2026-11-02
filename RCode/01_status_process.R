
##  01_status_process.R
##  The status process beta of that controls the regime switching in the process.

##  Two constructions are provided:
##    (a) sim_beta_renewal()  -- literal transcription of Definition 1;
##    (b) sim_beta_chain()    -- exact simulation on a fixed grid using the
##        transition matrix of the equivalent two-state Markov chain.
##  Construction (b) is exact at the grid points, is O(N) per path with no
##  rejection, and is what the BSDE solvers use.  Construction (a) is retained
##  because it is the definition, and the two are cross-checked in
##  10_validation.R.
#-------------------------------------------------------------------------------

## -----------------------------------------------------------------------------

## Generator of the two-state chain.

# -----------------------------------------------------------------------------
beta_transition <- function(p, dt) {
  s <- p$aU + p$aV
  e <- exp(-s * dt)
  list(p01 = (p$aU / s) * (1 - e),
       p10 = (p$aV / s) * (1 - e),
       pi1 = p$aU / s)
}

# -----------------------------------------------------------------------------

## sim_beta_chain(): exact grid simulation of beta for M independent paths.
##
## Returns an M x (N+1) matrix B with B[, n] = beta at time t_{n-1}, and an
## M x N matrix dJ counting the switches observed between consecutive grid points.


# -----------------------------------------------------------------------------
sim_beta_chain <- function(p, M, N, start = 0) {
  dt <- p$Tmat / N
  tr <- beta_transition(p, dt)
  B  <- matrix(0, nrow = M, ncol = N + 1)
  dJ <- matrix(0, nrow = M, ncol = N)
  B[, 1] <- start
  for (n in seq_len(N)) {
    cur <- B[, n]
    u   <- runif(M)
    ## probability of being in the other state at the end of the step
    psw <- ifelse(cur == 0, tr$p01, tr$p10)
    sw  <- u < psw
    B[, n + 1] <- ifelse(sw, 1 - cur, cur)
    dJ[, n]    <- as.numeric(sw)
  }
  list(B = B, dJ = dJ)
}

# -----------------------------------------------------------------------------

## compensated_jumps(): the increments of the compensated jump counting process
##
##      Jtilde_t = J_t - int_0^t lambda^J_s ds,
##
## where J counts every switch of beta (Definition 5) and, for the two-state
## chain, the intensity of the next switch is the rate of leaving the current
## state:
##      lambda^J_t = aU 1{beta_t = 0} + aV 1{beta_t = 1}.

# -----------------------------------------------------------------------------
compensated_jumps <- function(p, B, dJ, N) {
  dt  <- p$Tmat / N
  lam <- ifelse(B[, seq_len(N), drop = FALSE] == 0, p$aU, p$aV)
  dJ - lam * dt
}

# -----------------------------------------------------------------------------

## sim_beta_renewal(): literal Definition 1, for validation only.
## Returns beta sampled on the grid t_0, ..., t_N of [0, Tmat].

# -----------------------------------------------------------------------------
sim_beta_renewal <- function(p, M, N, start = 0) {
  tg <- seq(0, p$Tmat, length.out = N + 1)
  B  <- matrix(0, nrow = M, ncol = N + 1)
  for (m in seq_len(M)) {
    state <- start
    tcur  <- 0
    times <- numeric(0)
    states <- numeric(0)
    while (tcur < p$Tmat) {
      rate  <- if (state == 0) p$aU else p$aV
      tcur  <- tcur + rexp(1, rate)
      state <- 1 - state
      times  <- c(times, tcur)
      states <- c(states, state)
    }
    path <- rep(start, N + 1)
    for (k in seq_along(times)) path[tg >= times[k]] <- states[k]
    B[m, ] <- path
  }
  B
}

# -----------------------------------------------------------------------------

## Closed-form benchmarks for the two-state chain, used in 10_validation.R.
##   P(beta_t = 1 | beta_0 = 0) = (aU/s)(1 - e^{-s t})
##   E[number of switches on [0,t] | beta_0 = 0] is obtained by integrating the
##   state-dependent intensity against the transition law.

# -----------------------------------------------------------------------------
beta_prob_crisis <- function(p, t, start = 0) {
  s <- p$aU + p$aV
  if (start == 0) (p$aU / s) * (1 - exp(-s * t)) else 1 - (p$aV / s) * (1 - exp(-s * t))
}

beta_expected_jumps <- function(p, t, start = 0) {
  ## E[J_t] = int_0^t ( aU P(beta_u = 0) + aV P(beta_u = 1) ) du,
  ## with P(beta_u = 1) = pi1 + (b0 - pi1) e^{-s u} and b0 the initial state.
  s   <- p$aU + p$aV
  pi1 <- p$aU / s
  b0  <- as.numeric(start == 1)
  int_q1 <- pi1 * t + (b0 - pi1) * (1 - exp(-s * t)) / s
  p$aU * (t - int_q1) + p$aV * int_q1
}
