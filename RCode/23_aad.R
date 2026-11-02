# ------------------------------------------------------------------------------

##  23_aad.R
##  Suite (g): the adjoint (AAD) sensitivities of the quantisation solver.
##
##  The five sensitivities of 14_quantisation_adjoint.R -- S0 (delta), sigma
##  (vega), r_f^-, h_I, h_C -- on the six scenarios of 15 (both sides; calm,
##  crisis, dynamic; alpha = 0.5, r_f^- = 0.10), built from the base parameter
##  object.
##
##  PARTS.
##    1. discrete: the full adjoint against central bump-and-revalue of the SAME
##       scheme (grid rebuilt for every bump of S0 and sigma) at the three
##       jointly refined levels; agreement to the accuracy of the differences
##       shows the adjoint code is right at every level, not only at (50, 100).
##    2. continuous: the adjoint at each level and its Richardson extrapolation
##       against central bumps of the PDE on 1601 x 3000 (each bump of S0 or
##       sigma gets its own centred operator) and, for the seller in a frozen
##       regime, of the closed form extrapolated in its quadrature nodes (the
##       default quadrature of 08 costs about 0.1 per cent in a derivative,
##       Module 17; 15 is left as it was).
##    3. vega with the grid frozen against the full adjoint, across alpha, with
##       PDE bumps of vega as the reference (the record's Figure 7, all regimes).
##    4. cost: one cold valuation, valuation plus adjoint (five sensitivities),
##       and central bump-and-revalue, at every level; caches cleared before
##       each measurement and run serially so the timings do not compete.
##    5. delta cross-check: by the Feynman-Kac representation Zbreve = sigma S
##       d_S u, so dXVA_0/dS0 = Zbreve_0 / (sigma S0), the XVA part of the stock
##       position.  The adjoint delta is compared with the solver's own
##       Zbreve_0 / (sigma S0) at every level, after Richardson, and with the
##       PDE bump.
##  The buyer's rows rest on the buyer close-out convention (flagged).
##  Extended (approved): alpha, L_I and L_C are differentiated as well, so every
##  part covers eight sensitivities; like r_f^-, h_I and h_C they enter the
##  driver only and reuse the grid when bumped.
##
##  OUTPUT (output/data/): aad_discrete.csv, aad_continuous.csv,
##  aad_vega_alpha.csv, aad_cost.csv, aad_delta.csv.
##
##  USAGE.  Rscript main.R aad [force]
# ------------------------------------------------------------------------------

AAD_SETTINGS <- list(
  sides   = c("+", "-"),
  regimes = c("frozen0", "frozen1", "dynamic"),
  alpha   = 0.5, rf_m = 0.10,
  levels  = list(c(25L, 50L), c(50L, 100L), c(100L, 200L)),
  h_discrete   = c(S0 = 1e-4, sigma = 1e-4, rf_m = 1e-4, hI = 1e-4, hC = 1e-4,
                   alpha = 1e-4, LI = 1e-4, LC = 1e-4),
  h_continuous = c(S0 = 1e-3, sigma = 1e-3, rf_m = 1e-3, hI = 1e-3, hC = 1e-3,
                   alpha = 1e-3, LI = 1e-3, LC = 1e-3),
  pde     = c(1601L, 3000L),
  cf_ns   = c(3200L, 6400L),
  vega_alphas = seq(0, 1, by = 0.1),
  vega_cases  = list(c(side = "+", regime = "frozen0"), c(side = "+", regime = "frozen1"),
                     c(side = "+", regime = "dynamic"), c(side = "-", regime = "frozen0")),
  vega_level  = c(50L, 100L),
  seed = 23L)

AAD_USED_BY <- paste("Ch5 AAD tables: adjoint against same-scheme and PDE/closed-form bumps",
  "by level; frozen-versus-full vega figure across alpha; cost table; delta cross-check")

aad_scenarios <- function(p0, s) {
  g <- expand.grid(side = s$sides, regime = s$regimes, stringsAsFactors = FALSE)
  g$name <- sprintf("%s %s", ifelse(g$side == "+", "seller", "buyer"), g$regime)
  g
}
aad_params <- function(p0, side, regime, s, alpha = s$alpha)
  with_params(p0, side = side, beta_mode = regime, alpha = alpha, rf_m = s$rf_m)

## The bumped copies of a parameter object, S0 relative, the others absolute.
bumped <- function(p, nm, h, sgn) {
  if (nm == "S0") with_params(p, S0 = p$S0 * (1 + sgn * h)) else
    do.call(with_params, c(list(p), setNames(list(p[[nm]] + sgn * h), nm)))
}
central <- function(p, nm, h, f) {
  d <- if (nm == "S0") 2 * p$S0 * h else 2 * h
  (f(bumped(p, nm, h, 1)) - f(bumped(p, nm, h, -1))) / d
}

## PDE operators: the base one and one centred on each bumped S0 and sigma.
pde_ops <- function(p0, s) {
  h <- s$h_continuous
  op <- function(p) pde_operator(p, s$pde[1], s$pde[2])
  list(base = op(p0),
       S0_up = op(bumped(p0, "S0", h[["S0"]], 1)), S0_dn = op(bumped(p0, "S0", h[["S0"]], -1)),
       sigma_up = op(bumped(p0, "sigma", h[["sigma"]], 1)),
       sigma_dn = op(bumped(p0, "sigma", h[["sigma"]], -1)))
}
pde_op_for <- function(ops, p, p0) {
  if (p$S0 > p0$S0) return(ops$S0_up)
  if (p$S0 < p0$S0) return(ops$S0_dn)
  if (p$sigma > p0$sigma) return(ops$sigma_up)
  if (p$sigma < p0$sigma) return(ops$sigma_dn)
  ops$base
}
pde_greeks <- function(p, p0, ops, h) {
  f <- function(q) pde_value(q, pde_op_for(ops, q, p0))
  vapply(names(h), function(nm) central(p, nm, h[[nm]], f), 0)
}

## Closed form extrapolated in the quadrature nodes, on a fixed branch.
cf_value_ext <- function(p, br, ns) {
  v <- vapply(ns, function(k) linear_benchmark(p, branch = br, ns = k)$value, 0)
  2 * v[2] - v[1]
}
cf_greeks <- function(p, h, ns) {
  lb <- linear_benchmark(p)
  br <- if (lb$rho_f == p$rf_p) "rf_plus" else "rf_minus"
  vapply(names(h), function(nm) central(p, nm, h[[nm]], function(q) cf_value_ext(q, br, ns)), 0)
}

# -----------------------------------------------------------------------------

## Parts 1, 2 and 5 for one scenario.

# -----------------------------------------------------------------------------
aad_scenario <- function(sc, p0, s, ops) {
  p <- aad_params(p0, sc$side, sc$regime, s)
  flag <- if (sc$side == "-") "buyer close-out convention" else ""
  par <- names(s$h_discrete)
  adj <- lapply(s$levels, function(lv) quantisation_adjoint(p, lv[1], lv[2], grid = "full",
                                                             params = par))
  fwd <- lapply(s$levels, function(lv) quantisation_bsde(p, lv[1], lv[2]))
  lab <- vapply(s$levels, function(lv) sprintf("%d_%d", lv[1], lv[2]), "")

  ## part 1: same-scheme bumps at every level
  disc <- do.call(rbind, lapply(seq_along(s$levels), function(j) {
    lv <- s$levels[[j]]
    fd <- vapply(par, function(nm) central(p, nm, s$h_discrete[[nm]],
                 function(q) quantisation_bsde(q, lv[1], lv[2])$U0), 0)
    a <- adj[[j]]$grads[par]
    data.frame(scenario = sc$name, side = sc$side, regime = sc$regime, n = lv[1], N = lv[2],
               parameter = par, adjoint = a, bump = fd, rel_diff = rel_diff(a, fd),
               flag = flag, row.names = NULL, stringsAsFactors = FALSE)
  }))

  ## part 2: continuous references
  pg <- pde_greeks(p, p0, ops, s$h_continuous)
  cg <- if (sc$side == "+" && sc$regime != "dynamic") cf_greeks(p, s$h_continuous, s$cf_ns) else
          setNames(rep(NA_real_, length(par)), par)
  A  <- sapply(adj, function(a) a$grads[par])                 # parameters x levels
  rich <- 2 * A[, 3] - A[, 2]
  ref  <- ifelse(is.na(cg), pg, cg)
  cont <- data.frame(scenario = sc$name, side = sc$side, regime = sc$regime, parameter = par,
                     pde_bump = pg, closed_form_bump = cg,
                     reference = ifelse(is.na(cg), "pde", "closed_form"),
                     ## NA where the closed-form sensitivity vanishes on its branch
                     pde_vs_closed_form = ifelse(abs(cg) > 1e-10, pg / cg - 1, NA_real_),
                     flag = flag,
                     row.names = NULL, stringsAsFactors = FALSE)
  for (j in seq_along(lab)) {
    cont[[paste0("adjoint_", lab[j])]] <- A[, j]
    cont[[paste0("relerr_", lab[j])]]  <- rel_diff(A[, j], ref)
  }
  cont$adjoint_rich <- rich
  cont$relerr_rich  <- rel_diff(rich, ref)
  cont$adjoint_frozen_grid_50_100 <- quantisation_adjoint(p, 50L, 100L, grid = "frozen",
                                                          params = par)$grads[par]

  ## part 5: delta against the solver's own Zbreve_0 / (sigma S0)
  sS <- p$sigma * p$S0
  dl <- vapply(adj, function(a) a$grads[["S0"]], 0)
  zs <- vapply(fwd, function(q) q$Z0 / sS, 0)
  delta <- data.frame(scenario = sc$name, side = sc$side, regime = sc$regime,
                      level = c(lab, "richardson"),
                      adjoint_delta = c(dl, 2 * dl[3] - dl[2]),
                      Z0_over_sigmaS0 = c(zs, 2 * zs[3] - zs[2]),
                      delta_root = c(vapply(adj, function(a) a$delta_root, 0), NA),
                      pde_delta = pg[["S0"]], flag = flag, stringsAsFactors = FALSE)
  delta$adjoint_vs_Z <- delta$adjoint_delta - delta$Z0_over_sigmaS0
  delta$adjoint_vs_pde_rel <- delta$adjoint_delta / delta$pde_delta - 1
  delta$Z_vs_pde_rel <- delta$Z0_over_sigmaS0 / delta$pde_delta - 1
  ## the hedge's stock position xi^S_0 = (Zbreve_0 + Zhat_0)/(sigma S0) for the
  ## seller; the buyer's reflected hedge holds the negative of the buyer's total
  sg <- if (sc$side == "+") 1 else -1
  delta$total_shares <- sg * (delta$Z0_over_sigmaS0 + bs_control(0, p$S0, p) / sS)
  list(discrete = disc, continuous = cont, delta = delta)
}

# -----------------------------------------------------------------------------

## Part 3: frozen-grid against full vega across alpha.

# -----------------------------------------------------------------------------
aad_vega_case <- function(cs, p0, s, ops) {
  p  <- aad_params(p0, cs$side, cs$regime, s, alpha = cs$alpha)
  lv <- s$vega_level
  full <- quantisation_adjoint(p, lv[1], lv[2], grid = "full")
  frz  <- quantisation_adjoint(p, lv[1], lv[2], grid = "frozen")
  pv   <- central(p, "sigma", s$h_continuous[["sigma"]],
                  function(q) pde_value(q, pde_op_for(ops, q, p0)))
  data.frame(side = cs$side, regime = cs$regime, alpha = cs$alpha,
             vega_full = full$grads[["sigma"]], vega_frozen = frz$grads[["sigma"]],
             grid_contribution = full$vega_grid, vega_pde = pv,
             full_vs_pde = full$grads[["sigma"]] / pv - 1,
             frozen_vs_pde = frz$grads[["sigma"]] / pv - 1,
             ## absolute errors: the relative ones blow up where vega crosses zero
             full_minus_pde = full$grads[["sigma"]] - pv,
             frozen_minus_pde = frz$grads[["sigma"]] - pv,
             flag = if (cs$side == "-") "buyer close-out convention" else "",
             stringsAsFactors = FALSE)
}

# -----------------------------------------------------------------------------

## Part 4: cost at every level (cold caches, serial).

# -----------------------------------------------------------------------------
aad_cost <- function(p, s) {
  el <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr); proc.time()[["elapsed"]] - t0 }
  do.call(rbind, lapply(s$levels, function(lv) {
    v  <- function(q) quantisation_bsde(q, lv[1], lv[2])$U0
    clear_quant_caches(); t_val <- el(v(p))
    clear_quant_caches(); t_adj <- el(quantisation_adjoint(p, lv[1], lv[2], grid = "full",
                                                           params = names(s$h_discrete)))
    clear_quant_caches(); t_bmp <- el(vapply(names(s$h_discrete), function(nm)
                                       central(p, nm, s$h_discrete[[nm]], v), 0))
    clear_quant_caches()
    data.frame(n = lv[1], N = lv[2],
               method = c("one valuation",
                          sprintf("valuation + adjoint (%d sensitivities)", length(s$h_discrete)),
                          sprintf("central bump-and-revalue (%d sensitivities)", length(s$h_discrete))),
               seconds = c(t_val, t_adj, t_bmp), in_valuations = c(t_val, t_adj, t_bmp) / t_val)
  }))
}

# -----------------------------------------------------------------------------

## run_aad(): the suite.

# -----------------------------------------------------------------------------
run_aad <- function(p0 = make_params(), s = AAD_SETTINGS, force = FALSE) {
  run_suite("aad", "23_aad.R", p0 = p0, settings = s, seed = s$seed,
            used_by = AAD_USED_BY, force = force, run = function() {
    sc0  <- aad_scenarios(p0, s)
    cost <- aad_cost(aad_params(p0, "+", "frozen1", s), s)       # serial, cold
    ## trees for the base and every bumped S0 and sigma, before forking
    for (lv in s$levels) for (nm in c("S0", "sigma")) for (h in c(s$h_discrete[[nm]], s$h_continuous[[nm]]))
      for (sg in c(-1, 1)) invisible(get_nodes(bumped(p0, nm, h, sg), lv[1], lv[2]))
    for (lv in s$levels) invisible(get_nodes(p0, lv[1], lv[2]))
    ops  <- pde_ops(p0, s)
    parts <- par_lapply(seq_len(nrow(sc0)), function(i) aad_scenario(sc0[i, ], p0, s, ops))
    vg <- do.call(rbind, lapply(s$vega_cases, function(v)
            data.frame(side = v[["side"]], regime = v[["regime"]], alpha = s$vega_alphas,
                       stringsAsFactors = FALSE)))
    vega <- do.call(rbind, par_lapply(seq_len(nrow(vg)), function(i) aad_vega_case(vg[i, ], p0, s, ops)))
    pick <- function(nm) do.call(rbind, lapply(parts, `[[`, nm))
    list(discrete = pick("discrete"), continuous = pick("continuous"),
         vega_alpha = vega, cost = cost, delta = pick("delta"))
  })
}

## Console report used by main.R.
report_aad <- function(res) {
  d <- res$discrete
  cat("\npart 1: adjoint against same-scheme bumps, largest |relative difference| per level\n")
  a <- aggregate(list(max_rel = abs(d$rel_diff)),
                 by = list(parameter = d$parameter, level = paste(d$n, d$N, sep = "_")), max)
  print(reshape(a, idvar = "parameter", timevar = "level", direction = "wide"),
        row.names = FALSE, digits = 2)
  c2 <- res$continuous
  cat("\npart 2: relative error against the continuous reference (closed form where available, else PDE)\n")
  show <- c2[, c("scenario", "parameter", "reference", grep("^relerr_", names(c2), value = TRUE))]
  show[, -(1:3)] <- lapply(show[, -(1:3)], function(x) sprintf("%+.3f%%", 100 * x))
  ## sensitivities that vanish identically on their branch (seller L_C, buyer
  ## L_I, calm seller r_f^-) have no meaningful relative error
  refv <- ifelse(is.na(c2$closed_form_bump), c2$pde_bump, c2$closed_form_bump)
  print(show[abs(refv) > 1e-9, ], row.names = FALSE)
  cat(sprintf("max |PDE bump / closed-form bump - 1| = %.2e\n",
              max(abs(c2$pde_vs_closed_form), na.rm = TRUE)))
  v <- res$vega_alpha
  cat("\npart 3: vega across alpha (error against the PDE bump)\n")
  print(aggregate(cbind(full_abs = abs(full_minus_pde), frozen_abs = abs(frozen_minus_pde),
                        max_abs_vega = abs(vega_pde)) ~ side + regime, v, max),
        row.names = FALSE, digits = 3)
  sg <- v[sign(v$vega_frozen) != sign(v$vega_pde), c("side", "regime", "alpha")]
  if (nrow(sg)) { cat("frozen vega has the wrong sign at:\n"); print(sg, row.names = FALSE) }
  cat("\npart 4: cost (seller, crisis; cold caches)\n")
  print(res$cost, row.names = FALSE, digits = 3)
  dl <- res$delta
  cat("\npart 5: delta cross-check, adjoint delta against Zbreve_0 / (sigma S0)\n")
  print(dl[, c("scenario", "level", "adjoint_delta", "Z0_over_sigmaS0", "adjoint_vs_pde_rel",
               "Z_vs_pde_rel")], row.names = FALSE, digits = 4)
  cat(sprintf("\nsuite time %.0f s on %d core(s)\n", res$meta$seconds, res$meta$cores))
}
