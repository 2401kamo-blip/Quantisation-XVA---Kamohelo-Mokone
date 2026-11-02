#-------------------------------------------------------------------------------
##  05_network.R

##  The neural-network machinery of E, Han and Jentzen (2017), Section 4.1,
##  implemented in base R with hand-coded backpropagation.  No external package
##  is used, so no automatic differentiation is available and every gradient in
##  this file is derived by hand.  This is deliberate: the same derivation is
##  what will later be replaced by adjoint algorithmic differentiation.
##
##  ARCHITECTURE (E, Han and Jentzen 2017, Section 4.1 and Remark 4.1)
##  Each of the N-1 subnetworks has four layers: an input layer of dimension
##  d_in, two hidden layers of dimension H = d_in + 10, and an output layer of
##  dimension d_out.  Batch normalisation is applied immediately BEFORE the
##  first affine map and immediately AFTER each of the three affine maps; the
##  rectifier x -> max(x, 0) is the activation of the two hidden layers.  There
##  are no bias vectors, because the shift parameter of each batch-normalisation
##  layer already plays that role.  The parameter count is therefore
##      2 d_in                                  (input normalisation)
##    + d_in H + 2 H                            (first affine map + its BN)
##    + H^2 + 2 H                               (second affine map + its BN)
##    + H d_out + 2 d_out                       (third affine map + its BN)
##  which for d_in = d_out = d and H = d + 10 reproduces Equation (4.5) of E,
##  Han and Jentzen exactly.
##
##  In this application d_in is 1 when beta is frozen (the state is the stock
##  price alone) and 2 when beta is dynamic (the state is the pair (S, beta)).

#------------------------------------------------------------------------------

BN_EPS <- 1e-6   # numerical floor inside the normalisation
BN_MOM <- 0.99   # momentum of the moving averages used at evaluation time

## Names of the trainable arrays of a subnetwork.  m0..v3 are the moving
## averages, which are updated by the forward pass but never by the optimiser.
TRAIN_KEYS <- c("W1", "W2", "W3",
                "g0", "b0", "g1", "b1", "g2", "b2", "g3", "b3")

# -----------------------------------------------------------------------------

## init_subnet(): He initialisation of the affine maps, identity initialisation
## of the normalisation layers.
##
## out_scale sets the initial scale parameter of the FINAL normalisation layer,
## and therefore the order of magnitude at which the network output starts.  The
## default of one reproduces E, Han and Jentzen.  A smaller value is needed when
## the quantity being learnt is much smaller than one and its training signal is
## rare, which is the case for the jump controls of 11_deep_bsde_jumps.R; see
## the discussion there.

# -----------------------------------------------------------------------------
init_subnet <- function(d_in, d_out, H = d_in + 10L, out_scale = 1) {
  rW <- function(a, b) matrix(rnorm(a * b, sd = sqrt(2 / a)), nrow = a, ncol = b)
  list(
    W1 = rW(d_in, H), W2 = rW(H, H), W3 = rW(H, d_out),
    g0 = rep(1, d_in), b0 = rep(0, d_in),
    g1 = rep(1, H),    b1 = rep(0, H),
    g2 = rep(1, H),    b2 = rep(0, H),
    g3 = rep(out_scale, d_out), b3 = rep(0, d_out),
    m0 = rep(0, d_in), v0 = rep(1, d_in),
    m1 = rep(0, H),    v1 = rep(1, H),
    m2 = rep(0, H),    v2 = rep(1, H),
    m3 = rep(0, d_out), v3 = rep(1, d_out)
  )
}

# -----------------------------------------------------------------------------

## Batch normalisation, forward.
##   xhat = (x - mu) / sqrt(var + eps),   out = gamma * xhat + delta
## The batch statistics mu and var are the POPULATION moments of the mini-batch,
## matching the convention of the reference TensorFlow implementation.

# -----------------------------------------------------------------------------
bn_forward <- function(x, gamma, delta, mv, vv, train = TRUE) {
  m <- nrow(x)
  if (train) {
    mu  <- colMeans(x)
    va  <- colMeans(x^2) - mu^2
    mv  <- BN_MOM * mv + (1 - BN_MOM) * mu
    vv  <- BN_MOM * vv + (1 - BN_MOM) * va
  } else {
    mu <- mv; va <- vv
  }
  iv   <- 1 / sqrt(va + BN_EPS)
  ## Column-wise operations by recycling rep(v, each = m) rather than sweep():
  ## the same elementwise arithmetic, without the aperm() overhead of sweep(),
  ## which dominated the training time.
  xhat <- (x - rep(mu, each = m)) * rep(iv, each = m)
  out  <- xhat * rep(gamma, each = m) + rep(delta, each = m)
  list(out = out, cache = list(xhat = xhat, iv = iv, gamma = gamma),
       mv = mv, vv = vv)
}

# -----------------------------------------------------------------------------

## Batch normalisation, backward.  With m the batch size,
##   dgamma = sum_i dy_i xhat_i,   ddelta = sum_i dy_i,
##   dx = (gamma iv / m) [ m dy - sum_i dy_i - xhat sum_i (dy_i xhat_i) ].

# -----------------------------------------------------------------------------
bn_backward <- function(dy, cache) {
  m    <- nrow(dy)
  k    <- ncol(dy)
  xhat <- cache$xhat
  dgam <- colSums(dy * xhat)
  ddel <- colSums(dy)
  s1   <- matrix(ddel, nrow = m, ncol = k, byrow = TRUE)
  s2   <- matrix(dgam, nrow = m, ncol = k, byrow = TRUE)
  dx   <- (m * dy - s1 - xhat * s2) * rep(cache$gamma * cache$iv / m, each = m)
  list(dx = dx, dgamma = dgam, ddelta = ddel)
}

# -----------------------------------------------------------------------------

## subnet_forward(): returns the network output together with everything the
## backward pass needs, and the updated moving averages.

# -----------------------------------------------------------------------------
subnet_forward <- function(P, x, train = TRUE) {
  f0 <- bn_forward(x, P$g0, P$b0, P$m0, P$v0, train); h0 <- f0$out
  a1 <- h0 %*% P$W1
  f1 <- bn_forward(a1, P$g1, P$b1, P$m1, P$v1, train); r1 <- pmax(f1$out, 0)
  a2 <- r1 %*% P$W2
  f2 <- bn_forward(a2, P$g2, P$b2, P$m2, P$v2, train); r2 <- pmax(f2$out, 0)
  a3 <- r2 %*% P$W3
  f3 <- bn_forward(a3, P$g3, P$b3, P$m3, P$v3, train)

  if (train) {
    P$m0 <- f0$mv; P$v0 <- f0$vv; P$m1 <- f1$mv; P$v1 <- f1$vv
    P$m2 <- f2$mv; P$v2 <- f2$vv; P$m3 <- f3$mv; P$v3 <- f3$vv
  }
  list(y = f3$out, P = P,
       cache = list(h0 = h0, c0 = f0$cache, c1 = f1$cache,
                    r1 = r1, c2 = f2$cache, r2 = r2, c3 = f3$cache))
}

# -----------------------------------------------------------------------------

## subnet_backward(): gradients of the loss with respect to every trainable
## array, given the gradient dy with respect to the network output.

# -----------------------------------------------------------------------------
subnet_backward <- function(P, cache, dy) {
  G <- list()

  b3 <- bn_backward(dy, cache$c3)
  G$g3 <- b3$dgamma; G$b3 <- b3$ddelta
  G$W3 <- t(cache$r2) %*% b3$dx
  dr2  <- b3$dx %*% t(P$W3)

  dh2 <- dr2 * (cache$r2 > 0)                       # rectifier mask
  b2  <- bn_backward(dh2, cache$c2)
  G$g2 <- b2$dgamma; G$b2 <- b2$ddelta
  G$W2 <- t(cache$r1) %*% b2$dx
  dr1  <- b2$dx %*% t(P$W2)

  dh1 <- dr1 * (cache$r1 > 0)
  b1  <- bn_backward(dh1, cache$c1)
  G$g1 <- b1$dgamma; G$b1 <- b1$ddelta
  G$W1 <- t(cache$h0) %*% b1$dx
  dh0  <- b1$dx %*% t(P$W1)

  b0 <- bn_backward(dh0, cache$c0)
  G$g0 <- b0$dgamma; G$b0 <- b0$ddelta
  G
}

# -----------------------------------------------------------------------------

## Adam (Kingma and Ba 2015), exactly as in Example 5.2 of E, Han and Jentzen.
##   m_k = b1 m_{k-1} + (1-b1) grad
##   v_k = b2 v_{k-1} + (1-b2) grad^2
##   theta_k = theta_{k-1} - lr * (m_k / (1-b1^k)) / ( sqrt(v_k/(1-b2^k)) + eps )

# -----------------------------------------------------------------------------
adam_init <- function(par) {
  list(m = lapply(par, function(a) a * 0),
       v = lapply(par, function(a) a * 0),
       b1 = 0.9, b2 = 0.999, eps = 1e-8, k = 0L)
}

adam_step <- function(par, grad, state, lr) {
  state$k <- state$k + 1L
  for (nm in names(par)) {
    state$m[[nm]] <- state$b1 * state$m[[nm]] + (1 - state$b1) * grad[[nm]]
    state$v[[nm]] <- state$b2 * state$v[[nm]] + (1 - state$b2) * grad[[nm]]^2
    mh <- state$m[[nm]] / (1 - state$b1^state$k)
    vh <- state$v[[nm]] / (1 - state$b2^state$k)
    par[[nm]] <- par[[nm]] - lr * mh / (sqrt(vh) + state$eps)
  }
  list(par = par, state = state)
}
