# ------------------------------------------------------------------------------
##  26_plots.R
##  Figures and tables for Chapter 5, one function per item of the one-pager
##  "Chapter 5: Results. Revised structure, figures and tables" (14 Oct 2026).
##
##  This module only READS the suite outputs in output/data/ (Modules 17-25);
##  it never computes a value, so any item redraws in about a second.
##
##  RUNNING.
##    Rscript main.R plots                  every item of the one-pager
##    Rscript main.R plots fig02 tab03      only these items
##    Rscript main.R plots extras           the items no longer in the thesis
##  or, in RStudio, the block "Figures and tables, one at a time" in main.R:
##  run the lines that load the modules, then any single line such as fig02().
##
##  OUTPUT.
##    figures/figNN_<name>.pdf  vector figures at the thesis text width (16 cm =
##                              6.30 in), drawn with cairo_pdf so that the labels
##                              carry the thesis symbols; the proportions follow
##                              the 11 x 4.2 layout of the main.R figures, scaled
##                              to the text width so that the type prints at its
##                              true size;
##    tables/tabNN_<name>.tex   booktabs tabular fragments, without a float or a
##                              caption (the thesis wraps them);
##    figures/figure_map.csv    item, file, suite, data file, section, purpose.
##
##  STYLE.  Panel titles centred above the panels and a key for the solvers'
##  marks, as in the main.R figures.  Colour: the Okabe-Ito palette in a fixed
##  order, never cycled, with a line type and marker per series as well, so no
##  identity rests on colour alone; a series keeps its styling in every panel.
##  Heat maps use one hue from light to dark.
# ------------------------------------------------------------------------------

PLOT_SETTINGS <- list(
  width   = 6.30,          # inches: the text width of the thesis
  height1 = 2.70,          # one row of two panels (the main.R proportions, 11 : 4.2)
  height2 = 4.90,          # two rows
  pointsize = 9,
  family  = "sans",        # "Carlito" matches the thesis caption font if installed
  fig_dir = "figures",
  tab_dir = "tables")

PAL <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
LTY <- c(1, 2, 4, 5, 6, 3)
PCH <- c(16, 17, 15, 18, 1, 2)
INK <- "#222222"; MUTED <- "#666666"; GRIDC <- "#E3E3E3"

## Fixed roles: the regimes and the funding rates always take these slots.
REGIME_SLOT <- c("frozen0" = 1, "frozen1" = 2, "dynamic lamU=1.39" = 3, "dynamic lamU=0.30" = 4)
RF_LEVELS   <- c(0.05, 0.075, 0.10)                      # slots 1, 2, 3
regime_name <- function(r) unname(c("frozen0" = "normal", "frozen1" = "crisis",
                             "dynamic lamU=1.39" = "dynamic*','~lambda[U]=='1.39'",
                             "dynamic lamU=0.30" = "dynamic*','~lambda[U]=='0.30'")[r])
regime_tex  <- function(r) unname(c("frozen0" = "normal", "frozen1" = "crisis",
                             "dynamic lamU=1.39" = "dynamic, $\\lambda_U=1.39$",
                             "dynamic lamU=0.30" = "dynamic, $\\lambda_U=0.30$")[r])
lab_expr <- function(s) parse(text = s)
lam_lab  <- function(l) sprintf("%.2f", l)                       # lambda_U always to two decimals
rf_lab   <- function(r) ifelse(abs(r * 100 - round(r * 100)) < 1e-9, sprintf("%.2f", r), sprintf("%.3f", r))
rf_expr  <- function(r) lab_expr(paste0("r[f]^{phantom()-phantom()}=='", rf_lab(r), "'"))

# ---- data access ------------------------------------------------------------------
.suite_cache <- new.env()
suite_data <- function(name) {
  if (!is.null(.suite_cache[[name]])) return(.suite_cache[[name]])
  f <- file.path(OUTPUT_DATA, paste0(name, ".rds"))
  if (!file.exists(f)) stop("no data for suite '", name, "' (", f, ")")
  .suite_cache[[name]] <- readRDS(f)
}
base_params <- function() suite_data("verification")$meta$params
pct_Vhat0   <- function() { p <- base_params(); 100 / bs_value(0, p$S0, p) }   # currency -> % of Vhat_0

# ---- the item map -------------------------------------------------------------------
.fig_map <- new.env()
record_output <- function(item, file, suite, where, description) {
  assign(item, data.frame(item = item, file = file, suite = paste(suite, collapse = "+"),
                          data = paste(file.path(OUTPUT_DATA, paste0(suite, ".rds")), collapse = "+"),
                          where = where, description = description, drawn = format(Sys.time()),
                          stringsAsFactors = FALSE), envir = .fig_map)
}

# ---- figure helpers -------------------------------------------------------------------
fig_open <- function(file, rows = 1, cols = 2, height = NULL, s = PLOT_SETTINGS) {
  dir.create(s$fig_dir, showWarnings = FALSE, recursive = TRUE)
  h <- if (is.null(height)) (if (rows == 1) s$height1 else s$height2) else height
  cairo_pdf(file.path(s$fig_dir, file), width = s$width, height = h,
            pointsize = s$pointsize, family = s$family)
  par(mfrow = c(rows, cols), mar = c(3.5, 4.3, 2.2, 1.0), mgp = c(2.3, 0.6, 0),
      tcl = -0.3, las = 1, col.axis = INK, col.lab = INK, col.main = INK, fg = MUTED,
      cex.main = 0.95, cex.axis = 0.9, cex.lab = 1.0, font.main = 1)
}
fig_close <- function() invisible(dev.off())

## An empty panel with recessive grid and a centred title above it.
## Log-axis labels: plain numbers over a moderate range; over three decades or
## below 0.001, the powers of ten only, written as 10^k.
log_ticks <- function(at) {
  if (diff(range(log10(at))) >= 3 || min(at) < 1e-3) {
    at <- at[abs(log10(at) - round(log10(at))) < 1e-9]
    list(at = at, labels = parse(text = paste0("10^", round(log10(at)))))
  } else list(at = at, labels = formatC(at, format = "fg", digits = 3))
}
panel <- function(xlim, ylim, xlab, ylab, log = "", main = NULL, xat = NULL) {
  plot(NA, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, log = log, axes = FALSE)
  lab <- function(side) { if (grepl(c("x", "y")[side], log)) { tk <- log_ticks(axTicks(side))
    axis(side, at = tk$at, labels = tk$labels, col = MUTED) } else axis(side, col = MUTED) }
  if (is.null(xat)) lab(1) else axis(1, at = xat, labels = xat, col = MUTED)
  lab(2); box(col = MUTED)
  if (log == "") grid(col = GRIDC, lty = 1) else abline(h = axTicks(2), v = axTicks(1), col = GRIDC)
  if (!is.null(main)) title(main = main, line = 0.7)
}
series <- function(x, y, k, type = "l", lwd = 1.6, pch = NULL, cex = 0.8, lty = NULL, ...) {
  if (type %in% c("l", "b", "o")) lines(x, y, col = PAL[k], lty = if (is.null(lty)) LTY[k] else lty, lwd = lwd, ...)
  if (type %in% c("p", "b", "o")) points(x, y, col = PAL[k], pch = if (is.null(pch)) PCH[k] else pch, cex = cex)
}
leg <- function(pos, labels, k, type = "l", y = NULL, cex = 0.8, ...) {
  legend(pos, y, legend = labels, col = PAL[k], lty = if (type == "p") NA else LTY[k],
         pch = if (type == "l") NA else PCH[k], lwd = 1.6, bty = "n",
         seg.len = 2.4, text.col = INK, cex = cex, y.intersp = 1.15, ...)
}
## The key for the solvers' marks, as in the main.R figures.
solver_key <- function(pos, deep = TRUE, tied = FALSE) {
  lab <- c("PDE (line)", "quantisation (cross)", if (deep) "deep BSDE, mean of 4 seeds (circle)",
           if (tied) "deep, jump-tied \u03b6 = 10 (triangle)")
  legend(pos, legend = lab, pch = c(NA, 4, if (deep) 1, if (tied) 2), lty = c(1, NA, if (deep) NA, if (tied) NA),
         col = MUTED, text.col = MUTED, bty = "n", cex = 0.72, seg.len = 1.6)
}
errbars <- function(x, lo, hi, k) arrows(x, lo, x, hi, angle = 90, code = 3, length = 0.02,
                                         col = PAL[k], lwd = 1)

## ---- table helpers ----------------------------------------------------------------------
fmt <- function(x, d = 3, pct = FALSE, sign = FALSE) {
  v <- if (pct) 100 * x else x
  v[!is.na(v) & round(v, d) == 0] <- 0              # no signed zeros
  s <- formatC(v, format = "f", digits = d, flag = if (sign) "+" else "")
  s[!is.na(v) & v == 0] <- formatC(0, format = "f", digits = d)
  s <- gsub("^-", "$-$", gsub("^\\+", "$+$", s))
  s[is.na(x)] <- "--"
  s
}
fmte <- function(x, d = 2) {                    # 1.2e-05 -> $1.2\times10^{-5}$
  unname(vapply(x, function(v) {
    if (is.na(v)) return("--"); if (v == 0) return("0")
    e <- floor(log10(abs(signif(v, d)))); m <- signif(v, d) / 10^e
    sprintf("$%s\\times10^{%d}$", formatC(m, format = "f", digits = d - 1), e)
  }, ""))
}
## A booktabs tabular fragment.  A row whose first cell is "MIDRULE" becomes a
## \midrule; 'header2' (optional) is a first header line, e.g. \multicolumn groups.
tex_table <- function(file, df, align, header, note = NULL, header2 = NULL, s = PLOT_SETTINGS) {
  dir.create(s$tab_dir, showWarnings = FALSE, recursive = TRUE)
  body <- apply(df, 1, function(r) if (identical(unname(r[1]), "MIDRULE")) "\\midrule" else
    paste(paste(r, collapse = " & "), "\\\\"))
  lines <- c(sprintf("%% %s -- generated by R/26_plots.R on %s; do not edit by hand.",
                     file, format(Sys.Date())),
             if (!is.null(note)) paste("%", note),
             "\\setlength{\\tabcolsep}{4pt}% local to the enclosing table environment",
             sprintf("\\begin{tabular}{%s}", align), "\\toprule",
             if (!is.null(header2)) header2,
             paste(paste(header, collapse = " & "), "\\\\"), "\\midrule",
             body, "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(s$tab_dir, file))
}

# ------------------------------------------------------------------------------

##  FIGURES

# ------------------------------------------------------------------------------

# ---- Figure 1 (5.2): convergence of value and control -----------------------------
fig01 <- function() {
  sm <- suite_data("verification")$summary
  sm$rk <- ifelse(sm$regime == "dynamic", sprintf("dynamic lamU=%.2f", sm$mean_U), sm$regime)
  d <- sm[sm$side == "+" & sm$alpha == 0.5, ]
  n <- c(25, 50, 100)
  fig_open("fig01_convergence.pdf")
  for (q in c("U0", "Z0")) {
    x <- d[d$quantity == q, ]
    e <- abs(as.matrix(x[, c("relerr_25_50", "relerr_50_100", "relerr_100_200")])) * 100
    rr <- abs(x$relerr_rich_100_200) * 100
    yl <- range(c(e, rr[rr > 0]), na.rm = TRUE) * c(0.5, 2)
    panel(c(22, 130), yl, expression(n~"(with"~N==2*n*")"), "|relative error| (%)", log = "xy",
          xat = c(25, 50, 100),
          main = if (q == "U0") expression("(a) value"~XVA[0]^"+") else expression("(b) control"~"Z\u0306"[0]))
    for (i in seq_len(nrow(x))) {
      k <- REGIME_SLOT[[x$rk[i]]]
      series(n, e[i, ], k, type = "b")
      points(125, max(rr[i], yl[1] * 1.01), col = PAL[k], pch = c(1, 2, 0, 5)[k], cex = 0.9)
    }
    axis(1, at = 125, labels = "R", col = MUTED, font = 2)
    if (q == "U0") leg("bottomleft", lab_expr(regime_name(x$rk)), REGIME_SLOT[x$rk], type = "b")
  }
  fig_close()
  record_output("fig01", "fig01_convergence.pdf", "verification", "5.2",
                "relative error of value and control under joint refinement (seller, alpha = 0.5); R: Richardson")
}

# ---- Figure 2 (5.3): XVA against alpha, normal | crisis, lines by r_f^- -------------
fig02 <- function() {
  r <- suite_data("ps"); f4 <- r$figure4; dp <- r$deep_points; sc <- pct_Vhat0()
  f4$pde <- sc * f4$pde; f4$quantisation <- sc * f4$quantisation; dp$U0_mean <- sc * dp$U0_mean
  fig_open("fig02_xva_alpha_rf.pdf")
  yl <- range(c(f4$pde, dp$U0_mean[dp$figure == "4"])); yl[1] <- yl[1] - 0.16 * diff(yl)
  for (rg in c("frozen0", "frozen1")) {
    x <- f4[f4$regime == rg, ]
    panel(c(0, 1), yl, expression(alpha), expression(XVA[0]^"+"~"(% of"~hat(V)[0]*")"),
          main = if (rg == "frozen0") expression("(a) normal regime,"~beta[t]==0) else
                   expression("(b) crisis regime,"~beta[t]==1))
    for (j in 1:3) {
      y <- x[x$rf_m == RF_LEVELS[j], ]; series(y$alpha, y$pde, j)
      i <- seq(1, nrow(y), 2); points(y$alpha[i], y$quantisation[i], pch = 4, col = PAL[j], cex = 0.6)
      d <- dp[dp$figure == "4" & dp$regime == rg & dp$rf_m == RF_LEVELS[j], ]
      points(d$alpha, d$U0_mean, pch = 1, col = PAL[j], cex = 1)
    }
    if (rg == "frozen0") leg("topleft", rf_expr(RF_LEVELS), 1:3) else solver_key("bottomright")
  }
  fig_close()
  record_output("fig02", "fig02_xva_alpha_rf.pdf", "ps", "5.3",
                "seller XVA in % of Vhat_0 against alpha in the normal and crisis regimes, one line per r_f^-; at r_f^- = 0.05 the regimes coincide")
}

# ---- Figure 3 (5.3): XVA against alpha, calm against dynamic, lines by lambda_U -----
fig03 <- function() {
  r <- suite_data("ps"); f5 <- r$figure5; dp <- r$deep_points; sc <- pct_Vhat0()
  f5 <- f5[f5$convention == "mean", ]
  f5$pde <- sc * f5$pde; f5$quantisation <- sc * f5$quantisation; dp$U0_mean <- sc * dp$U0_mean
  fig_open("fig03_xva_alpha_lambda.pdf")
  yl <- range(c(f5$pde, f5$quantisation, dp$U0_mean[dp$figure == "5"]), na.rm = TRUE); yl[1] <- yl[1] - 0.22 * diff(yl)
  lam <- c(0.3, 1, 1.7)
  g <- suite_data("bands")$grid; g <- g[g$regime == "dynamic lamU=1.39", ]     # calibrated lambda_U, already in %
  for (rf in c(0.05, 0.10)) {
    x <- f5[f5$rf_m == rf, ]
    panel(c(0, 1), yl, expression(alpha), expression(XVA[0]^"+"~"(% of"~hat(V)[0]*")"),
          main = if (rf == 0.05) expression("(a)"~r[f]^{phantom()-phantom()}=="0.05"~"(benchmark)") else
                   expression("(b)"~r[f]^{phantom()-phantom()}=="0.10"))
    y <- x[x$regime == "frozen0", ]; series(y$alpha, y$pde, 1)
    for (j in 1:3) {
      y <- x[x$regime == "dynamic" & x$lambda_U == lam[j], ]; series(y$alpha, y$pde, j + 1)
      i <- seq(1, nrow(y), 2); points(y$alpha[i], y$quantisation[i], pch = 4, col = PAL[j + 1], cex = 0.6)
      d <- dp[dp$figure == "5" & dp$rf_m == rf & dp$lambda_U == lam[j], ]
      z <- d[d$solver == "plain", ]; points(z$alpha, z$U0_mean, pch = 1, col = PAL[j + 1])
      z <- d[d$solver == "tied", ];  points(z$alpha, z$U0_mean, pch = 2, col = PAL[j + 1])
    }
    y <- g[abs(g$rf_m - rf) < 1e-9, ]; y <- y[order(y$alpha), ]
    series(y$alpha, y$XVA_plus_pct_Vhat0, 5); i <- seq(1, nrow(y), 2)
    points(y$alpha[i], y$XVA_plus_pct_Vhat0[i], pch = 4, col = PAL[5], cex = 0.6)
    if (rf == 0.05) {
      leg("topleft", c(expression(beta[t]==0), lab_expr(paste0("lambda[U]=='", lam_lab(c(0.3, 1, 1.39, 1.7)), "'",
                                                             c("", "", "~'(calibrated)'", "")))), c(1, 2, 3, 5, 4))
      text(0.98, yl[1] + 0.03 * diff(yl), "all curves coincide", adj = c(1, 0), cex = 0.75, col = MUTED)
    } else solver_key("bottomright", tied = TRUE)
  }
  fig_close()
  record_output("fig03", "fig03_xva_alpha_lambda.pdf", c("ps", "bands"), "5.3",
                "seller XVA in % of Vhat_0 against alpha, frozen normal regime against the dynamic regime, one line per mean holding time lambda_U; lambda_U = 1.39 (calibrated) from the bands suite, quantisation only")
}

# ---- Figure 4 (5.3): (lambda_U, lambda_V) heat map ------------------------------------
fig04 <- function() {
  h <- suite_data("regime")$heatmap
  ramp <- colorRampPalette(c("#F2F7FC", "#9ECAE1", "#2171B5", "#08306B"))(64)
  fig_open("fig04_heatmap.pdf", height = 3.0)
  layout(matrix(1:4, 1), widths = c(1, 0.19, 1, 0.19)); par(cex = 0.85)
  lam <- sort(unique(h$lambda_U))
  for (b in 0:1) {
    x <- h[h$side == "+" & h$rf_m == 0.10 & h$beta0 == b, ]
    z <- matrix(x$XVA_pct_Vhat0[order(x$lambda_V, x$lambda_U)], length(lam), length(lam))
    zl <- range(z)
    par(mar = c(3.5, 4.0, 2.2, 0.4))
    image(log(lam), log(lam), z, col = ramp, axes = FALSE, zlim = zl, useRaster = TRUE,
          xlab = expression(lambda[U]~"(years)"), ylab = expression(lambda[V]~"(years)"))
    contour(log(lam), log(lam), z, add = TRUE, col = INK, labcex = 0.7, lwd = 0.7, nlevels = 8)
    at <- c(0.25, 0.5, 1, 2, 3); axis(1, log(at), at, col = MUTED); axis(2, log(at), at, col = MUTED); box(col = MUTED)
    points(log(1.39), log(0.99), pch = 4, col = "#000000", cex = 1.1, lwd = 1.5)
    title(main = if (b == 0) expression("(a) starting normal,"~beta[0]==0) else
                   expression("(b) starting in crisis,"~beta[0]==1), line = 0.7)
    par(mar = c(3.5, 0.5, 2.2, 2.4))                     # colour bar, % of Vhat_0
    zz <- seq(zl[1], zl[2], length.out = length(ramp) + 1)
    image(1, zz[-1] - diff(zz) / 2, matrix(zz[-1] - diff(zz) / 2, 1), col = ramp, zlim = zl,
          axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
    axis(4, col = MUTED, cex.axis = 0.8); box(col = MUTED)
    mtext(expression("% of"~hat(V)[0]), side = 3, line = 0.4, cex = 0.7, col = INK)
  }
  fig_close()
  record_output("fig04", "fig04_heatmap.pdf", "regime", "5.3",
                "seller XVA_0 in % of Vhat_0 over the mean holding times, r_f^- = 0.10, alpha = 0.5; cross: calibrated point; one colour scale per panel")
}

# ---- Figure 5 (App. C): regime effect against the spread ------------------------------
fig05 <- function() {
  sp <- suite_data("regime")$spread; sp <- sp[sp$side == "+", ]
  fig_open("fig05_regime_spread.pdf")
  for (col in c("effect_crisis_pct_Vhat0", "effect_dynamic_pct_Vhat0")) {
    panel(range(sp$spread) * 1e4, range(sp[[col]]),
          expression(r[f]^{phantom()-phantom()}-r[r]^{phantom()-phantom()}~"(bp)"),
          expression("regime effect (% of"~hat(V)[0]*")"),
          main = if (col == "effect_crisis_pct_Vhat0") "(a) crisis minus normal" else "(b) dynamic minus normal")
    for (j in 1:3) { a <- c(0, 0.5, 1)[j]; x <- sp[sp$alpha == a, ]; series(x$spread * 1e4, x[[col]], j) }
    if (col == "effect_crisis_pct_Vhat0") leg("topleft", lab_expr(paste0("alpha==", c(0, 0.5, 1))), 1:3)
  }
  fig_close()
  record_output("fig05", "fig05_regime_spread.pdf", "regime", "App. C",
                "regime effect against the funding spread, seller, calibrated holding times")
}

# ---- Figure 6 (5.4): seller and buyer against alpha by r_f^- (the band) ----------------
fig06 <- function() {
  g <- suite_data("bands")$grid
  fig_open("fig06_band.pdf")
  yl <- range(c(g$XVA_plus_pct_Vhat0, g$XVA_minus_pct_Vhat0)[g$regime %in% c("frozen0", "frozen1")])
  for (rg in c("frozen0", "frozen1")) {
    d <- g[g$regime == rg, ]
    panel(c(0, 1), yl, expression(alpha), expression(XVA[0]^"\u00b1"~"(% of"~hat(V)[0]*")"),
          main = if (rg == "frozen0") "(a) normal regime" else "(b) crisis regime")
    for (j in 1:3) {
      x <- d[d$rf_m == RF_LEVELS[j], ]
      lines(x$alpha, x$XVA_plus_pct_Vhat0, col = PAL[j], lty = 1, lwd = 1.6)
      lines(x$alpha, x$XVA_minus_pct_Vhat0, col = PAL[j], lty = 2, lwd = 1.6)
      i <- seq(1, nrow(x), 5); points(x$alpha[i], x$XVA_plus_pct_Vhat0[i], col = PAL[j], pch = PCH[j], cex = 0.7)
      points(x$alpha[i], x$XVA_minus_pct_Vhat0[i], col = PAL[j], pch = PCH[j], cex = 0.7)
    }
    if (rg == "frozen0") {
      legend("topleft", legend = c(rf_expr(RF_LEVELS), expression("seller,"~XVA^"+"), expression("buyer,"~XVA^"\u2212")),
             col = c(PAL[1:3], INK, INK), pch = c(PCH[1:3], NA, NA), lty = c(NA, NA, NA, 1, 2),
             lwd = 1.6, bty = "n", cex = 0.78, seg.len = 2.4)
    }
  }
  fig_close()
  record_output("fig06", "fig06_band.pdf", "bands", "5.4",
                "seller (solid) and buyer (dashed) XVA against alpha, one colour per r_f^-; buyer under the admissible crisis convention (short sale banned)")
}

# ---- Figure 7 (App. C): repo-rate and credit-intensity panels --------------------------
fig07 <- function() {
  r <- suite_data("bands")
  rp <- r$repo[r$repo$design == "symmetric spread", ]
  fig_open("fig07_repo_credit.pdf")
  d <- 1e4 * (rp$rr_m - 0.05)
  panel(range(d), range(c(rp$XVA_plus_pct_Vhat0, rp$XVA_minus_pct_Vhat0)),
        expression("repo spread"~d~"(bp):"~~r[r]^{phantom()-phantom()}==0.05+d*','~~r[r]^"+"==0.05-d),
        expression(XVA[0]~"(% of"~hat(V)[0]*")"), main = "(a) repo rates")
  x <- rp[rp$regime == "frozen0", ]; series(d[rp$regime == "frozen0"], x$XVA_plus_pct_Vhat0, 1, type = "b")
  x <- rp[rp$regime == "frozen1", ]; series(d[rp$regime == "frozen1"], x$XVA_plus_pct_Vhat0, 2, type = "b")
  x <- rp[rp$regime == "frozen0", ]; series(d[rp$regime == "frozen0"], x$XVA_minus_pct_Vhat0, 5, type = "b")
  x <- rp[rp$regime == "frozen1", ]; series(d[rp$regime == "frozen1"], x$XVA_minus_pct_Vhat0, 6, type = "b")
  leg(par("usr")[1], y = 17.2, c("seller, normal", "seller, crisis", "buyer, normal", "buyer, crisis"), c(1, 2, 5, 6), type = "b")
  cr <- r$credit[r$credit$regime == "frozen0", ]
  sc <- c(0.5, 1, 2)
  panel(range(cr$scale_hI), range(cr$width_pct_Vhat0),
        expression(h[I]~"as a multiple of its benchmark"),
        expression("band width"~XVA^"+"-XVA^"\u2212"~"(% of"~hat(V)[0]*")"), main = "(b) credit intensities, normal regime")
  abline(h = 0, col = MUTED, lwd = 0.8)
  for (j in 1:3) { x <- cr[cr$scale_hC == sc[j], ]; x <- x[order(x$scale_hI), ]; series(x$scale_hI, x$width_pct_Vhat0, j, type = "b") }
  leg("topright", lab_expr(paste0("h[C]==", sc, "%*%h[C]^{bench}")), 1:3, type = "b")
  text(par("usr")[1], par("usr")[3], "negative width: inverted band", adj = c(-0.05, -0.8), cex = 0.72, col = MUTED)
  fig_close()
  record_output("fig07", "fig07_repo_credit.pdf", "bands", "App. C",
                "alpha = 0.5, r_f^- = 0.10. (a) seller and buyer XVA, normal and crisis, against a symmetric repo spread d from 0 to 200 bp (r_r^- = 0.05 + d, r_r^+ = 0.05 - d); (b) normal-regime band width against h_I scaled 0.5 to 2 times its benchmark, one line per h_C at 0.5, 1 and 2 times its benchmark (negative: inverted band)")
}

# ---- Figure 8 (5.5): stock units against S and time, normal against crisis ------------
fig08 <- function() {
  nd <- suite_data("replication")$nodes; p <- base_params()
  ks <- c(10, 25, 45)                                   # tree dates t = k T / 50
  fig_open("fig08_stock_units.pdf")
  for (rg in c("frozen0", "frozen1")) {
    panel(c(0.6, 1.6), c(-0.02, 1.05), expression(S[t]), expression("stock position"~xi[t]^S~"(shares)"),
          main = if (rg == "frozen0") "(a) normal regime" else "(b) crisis regime")
    for (j in seq_along(ks)) {
      x <- nd[nd$side == "+" & nd$regime == rg & nd$k == ks[j] & nd$S > 0.6 & nd$S < 1.6, ]
      x <- x[order(x$S), ]; x <- x[-c(1, nrow(x)), ]     # outermost cells: tail conditional means
      lines(x$S, x$units_stock, col = PAL[j], lty = 1, lwd = 1.6)
      lines(x$S, bs_delta(x$t[1], x$S, p), col = PAL[j], lty = 3, lwd = 1.1)
    }
    if (rg == "frozen0")
      legend("topleft", legend = c(sprintf("t = %.3f", ks * p$Tmat / 50), "hedge (solid)", "Black-Scholes delta (dotted)"),
             col = c(PAL[seq_along(ks)], INK, INK), lty = c(1, 1, 1, 1, 3), lwd = c(1.6, 1.6, 1.6, 1.6, 1.1),
             bty = "n", cex = 0.75, seg.len = 2.2)
  }
  fig_close()
  record_output("fig08", "fig08_stock_units.pdf", "replication", "5.5",
                "seller's stock position against S at three dates, with the Black-Scholes delta; alpha = 0.5, r_f^- = 0.10; (50,100) tree without the two outermost cells")
}

# ---- Figure 9 (5.5): positions over alpha and r_f^- ----------------------------------
fig09 <- function() {
  a <- suite_data("bands")$accounts_t0; a <- a[a$side == "+" & a$regime %in% c("frozen0", "frozen1"), ]
  fig_open("fig09_positions.pdf", rows = 2, cols = 3, height = PLOT_SETTINGS$height2)
  par(mar = c(3.9, 4.6, 2.0, 0.6), mgp = c(2.9, 0.6, 0), cex = 0.85)
  get <- function(rg, rf, acc, what = "value") { x <- a[a$regime == rg & a$rf_m == rf & a$account %in% acc, ]
    tapply(x[[what]], x$alpha, sum) }
  ttl <- c(stock = "stock~(shares)", bonds = "bonds~P^I+P^C", funding = "funding~account")
  for (rg in c("frozen0", "frozen1")) for (cc in names(ttl)) {
    acc <- switch(cc, stock = "stock", bonds = c("investor_bond", "counterparty_bond"), funding = "funding")
    what <- if (cc == "stock") "units" else "value"
    ys <- sapply(RF_LEVELS, function(rf) get(rg, rf, acc, what)); al <- as.numeric(rownames(ys))
    yl <- range(ys); if (diff(yl) < 1e-3 * max(abs(yl))) yl <- yl + c(-1, 1) * 1e-3 * max(abs(yl))
    panel(c(0, 1), yl, expression(alpha), if (cc == "stock") expression(xi[0]^S) else "account value",
          main = lab_expr(paste0("'(", letters[(rg == "frozen1") * 3 + match(cc, names(ttl))], ") ",
                                 if (rg == "frozen0") "normal: " else "crisis: ", "'~", ttl[[cc]])))
    if (cc == "funding") abline(h = 0, col = MUTED, lwd = 0.8)
    for (j in 1:3) { series(al, ys[, j], j); i <- seq(1, length(al), 5); points(al[i], ys[i, j], col = PAL[j], pch = PCH[j], cex = 0.7) }
    if (rg == "frozen0" && cc == "stock") leg("topleft", rf_expr(RF_LEVELS), 1:3, type = "b", cex = 0.8)
  }
  fig_close()
  record_output("fig09", "fig09_positions.pdf", "bands", "5.5",
                "seller's stock units (shares) and the values of the bond and funding accounts at t = 0, per unit notional, over alpha, one line per r_f^-, normal and crisis")
}

# ---- Figure 10 (5.5): replication error and its decay ---------------------------------
fig10 <- function() {
  sm <- suite_data("replication")$sim
  d <- sm[sm$source == "pde", ]
  rk <- ifelse(d$regime == "dynamic", sprintf("dynamic lamU=%.2f", d$mean_U), d$regime)
  fig_open("fig10_replication_decay.pdf")
  for (pnl in 1:2) {
    yv <- if (pnl == 1) d$rmse_pct_Vhat0 else ifelse(d$regime == "dynamic", d$rmse_adj_pct_Vhat0, d$rmse_pct_Vhat0)
    panel(range(d$step), range(sm$rmse_pct_Vhat0, na.rm = TRUE) * c(0.8, 1.2),
          expression("rebalancing step"~Delta*t), expression("RMSE (% of"~hat(V)[0]*")"), log = "xy",
          main = if (pnl == 1) "(a) as hedged" else "(b) regime-switch martingale added back")
    for (g in unique(rk)) { i <- rk == g; o <- order(d$step[i]); series(d$step[i][o], yv[i][o], REGIME_SLOT[[g]], type = "b") }
    if (pnl == 1) leg("topleft", lab_expr(regime_name(unique(rk))), REGIME_SLOT[unique(rk)], type = "b")
  }
  fig_close()
  record_output("fig10", "fig10_replication_decay.pdf", "replication", "5.5",
                "RMSE of the seller's replication error against the rebalancing step; alpha = 0.5, r_f^- = 0.10; 20,000 paths on 1,600 fine steps, 25 to 1,600 rebalancing dates, hedge ratios from the PDE (801 space nodes); slopes in Table 11")
}

# ---- Figure 11 (5.6): stacked Shapley contributions -----------------------------------
fig11 <- function() {
  v <- suite_data("shapley")$values
  v <- v[v$version == "main" & v$alpha == 0.5 & v$rf_m == 0.10, ]
  rows <- list(c("+", "dynamic", "1.39", "seller*','~dynamic*','~lambda[U]=='1.39'"),
               c("+", "dynamic", "0.3", "seller*','~dynamic*','~lambda[U]=='0.30'"),
               c("+", "frozen1", NA, "seller*','~crisis"),
               c("-", "dynamic", "1.39", "buyer*','~dynamic*','~lambda[U]=='1.39'"))
  feats <- c("baseline", "funding_spread", "dva", "cva", "regime")
  cols  <- c("#BDBDBD", PAL[1:4])
  M <- sapply(rows, function(rw) { x <- v[v$side == rw[1] & v$regime == rw[2] &
                                         (is.na(rw[3]) | abs(v$mean_U - as.numeric(rw[3])) < 1e-9), ]
                                   x$pct_Vhat0[match(feats, x$feature)] })
  tot <- colSums(M)
  fig_open("fig11_shapley.pdf", cols = 1, height = 3.1)
  par(mar = c(3.4, 11.5, 2.6, 1.0))
  ny <- ncol(M); yy <- rev(seq_len(ny))
  plot(NA, xlim = c(min(0, colSums(pmin(M, 0))) - 0.5, max(colSums(pmax(M, 0))) * 1.08), ylim = c(0.4, ny + 0.6),
       axes = FALSE, xlab = expression("contribution (% of"~hat(V)[0]*")"), ylab = "")
  abline(v = axTicks(1), col = GRIDC); axis(1, col = MUTED); abline(v = 0, col = MUTED)
  axis(2, at = yy, labels = lab_expr(sapply(rows, `[`, 4)), tick = FALSE, las = 1)
  for (j in seq_len(ny)) {
    pos <- 0; neg <- 0
    for (f in seq_along(feats)) {
      val <- M[f, j]; if (abs(val) < 1e-9) next
      if (val > 0) { rect(pos, yy[j] - 0.32, pos + val, yy[j] + 0.32, col = cols[f], border = "white", lwd = 1.5); pos <- pos + val }
      else         { rect(neg + val, yy[j] - 0.32, neg, yy[j] + 0.32, col = cols[f], border = "white", lwd = 1.5); neg <- neg + val }
    }
    points(tot[j], yy[j], pch = 23, bg = INK, col = INK, cex = 0.9)
    text(max(pos, tot[j]), yy[j], sprintf("  %.2f", tot[j]), adj = 0, cex = 0.75, col = INK, xpd = NA)
  }
  legend(par("usr")[1], ny + 0.6 + 0.12 * ny, yjust = 0,
         legend = c("baseline", "funding spread", "DVA", "CVA", "regime", expression(XVA[0]~"(total)")),
         fill = c(cols, NA), border = c(rep(NA, 5), NA), pch = c(rep(NA, 5), 23), pt.bg = INK,
         col = c(rep(NA, 5), INK), bty = "n", cex = 0.8, horiz = TRUE, xpd = NA,
         x.intersp = 0.6, text.width = NA)
  fig_close()
  record_output("fig11", "fig11_shapley.pdf", "shapley", "5.6",
                "stacked Shapley contributions over 24 orderings, alpha = 0.5, r_f^- = 0.10; spread and regime are complements")
}

# ---- Figure 12 (5.7): frozen-grid against full-grid vega ---------------------------
fig12 <- function() {
  v <- suite_data("aad")$vega_alpha
  fig_open("fig12_vega.pdf")
  for (rg in c("frozen0", "frozen1")) {
    x <- v[v$side == "+" & v$regime == rg, ]
    panel(c(0, 1), range(c(x$vega_full, x$vega_frozen, x$vega_pde)) * 1e3, expression(alpha),
          expression(partialdiff*XVA[0]^"+"/partialdiff*sigma~"("%*%10^-3*")"),
          main = if (rg == "frozen0") "(a) normal regime" else "(b) crisis regime")
    abline(h = 0, col = MUTED, lwd = 0.8)
    series(x$alpha, 1e3 * x$vega_full, 1); series(x$alpha, 1e3 * x$vega_frozen, 2)
    points(x$alpha, 1e3 * x$vega_pde, pch = 1, col = INK, cex = 1)
    if (rg == "frozen0") legend("bottomright", c("adjoint, grid differentiated", "adjoint, grid frozen", "PDE bump"),
                                col = c(PAL[1:2], INK), lty = c(LTY[1:2], NA), pch = c(NA, NA, 1), lwd = 1.6, bty = "n", cex = 0.78)
  }
  fig_close()
  record_output("fig12", "fig12_vega.pdf", "aad", "5.7", "vega with the grid differentiated and frozen, adjoint at (50,100) (not extrapolated), against PDE bumps; seller, r_f^- = 0.10")
}

# ---- Figure 13 (5.7): AAD cost scaling and delta cross-check ------------------------
fig13 <- function() {
  r <- suite_data("aad"); ct <- r$cost; dl <- r$delta
  fig_open("fig13_aad_cost_delta.pdf")
  meth <- unique(ct$method)
  panel(c(22, 115), range(ct$seconds) * c(0.7, 1.4), expression(n~"(with"~N==2*n*")"), "seconds (cold caches)",
        log = "xy", xat = c(25, 50, 100), main = "(a) cost of eight sensitivities")
  for (j in seq_along(meth)) { x <- ct[ct$method == meth[j], ]; series(x$n, x$seconds, c(1, 3, 2)[j], type = "b") }
  leg("topleft", c("one valuation", "valuation + adjoint", "central bumping"), c(1, 3, 2), type = "b")
  d <- dl[dl$scenario == "seller frozen1", ]; lv <- d[d$level != "richardson", ]; ri <- d[d$level == "richardson", ]
  n <- c(25, 50, 100)
  e1 <- abs(lv$adjoint_vs_pde_rel) * 100; e2 <- abs(lv$Z_vs_pde_rel) * 100
  yl <- range(c(e1, e2, abs(ri$adjoint_vs_pde_rel) * 100, abs(ri$Z_vs_pde_rel) * 100)) * c(0.5, 2)
  panel(c(22, 130), yl, expression(n~"(with"~N==2*n*")"), "|relative error| against PDE (%)", log = "xy",
        xat = c(25, 50, 100), main = "(b) delta, seller in crisis")
  series(n, e1, 1, type = "b"); series(n, e2, 2, type = "b")
  points(125, abs(ri$adjoint_vs_pde_rel) * 100, pch = 1, col = PAL[1]); points(125, abs(ri$Z_vs_pde_rel) * 100, pch = 2, col = PAL[2])
  axis(1, at = 125, labels = "R", col = MUTED, font = 2)
  leg("bottomleft", c(expression("adjoint"~partialdiff*XVA[0]/partialdiff*S[0]), expression("Z\u0306"[0]/(sigma*S[0]))), 1:2, type = "b",
      inset = c(0.02, 0.02))
  fig_close()
  record_output("fig13", "fig13_aad_cost_delta.pdf", "aad", "5.7",
                "cost of the eight sensitivities by level; delta by the adjoint and from the control against the PDE bump (R: Richardson)")
}

# ---- Figure 14 (5.8): accuracy against run time, by method --------------------------
method_performance <- function() {
  vf <- suite_data("verification"); dp <- suite_data("deep")
  case <- "seller dynamic lamU=0.30 alpha=0.00 rf-=0.100"
  sm <- vf$summary[vf$summary$case == case, ]; va <- vf$values[vf$values$case == case, ]
  tr <- vf$trees
  q <- function(n, N, qty) sm[sm$quantity == qty, ][[sprintf("relerr_%d_%d", n, N)]]
  qsec <- function(n) va$seconds[va$method == "quantisation" & va$n == n & va$quantity == "U0"] +
                      tr$build_seconds[tr$n == n]
  rows <- list(
    data.frame(method = "quantisation (25,50)", group = "quant", U0 = q(25, 50, "U0"), Z0 = q(25, 50, "Z0"), sec = qsec(25)),
    data.frame(method = "quantisation (50,100)", group = "quant", U0 = q(50, 100, "U0"), Z0 = q(50, 100, "Z0"), sec = qsec(50)),
    data.frame(method = "quantisation (100,200)", group = "quant", U0 = q(100, 200, "U0"), Z0 = q(100, 200, "Z0"), sec = qsec(100)),
    data.frame(method = "quantisation, Richardson", group = "rich",
               U0 = sm$relerr_rich_100_200[sm$quantity == "U0"], Z0 = sm$relerr_rich_100_200[sm$quantity == "Z0"],
               sec = qsec(50) + qsec(100)),
    data.frame(method = "PDE 801 x 1500", group = "pde", U0 = sm$pde_std_vs_ref[sm$quantity == "U0"],
               Z0 = sm$pde_std_vs_ref[sm$quantity == "Z0"], sec = va$seconds[va$method == "pde_std"][1]))
  ru <- dp$runs[dp$runs$lambda_U == 0.3, ]
  dsel <- list(c("plain", NA, "deep, plain"), c("tied", 1000, "deep, tied (\u03b6 = 1000)"),
               c("untied", 1000, "deep, untied (\u03b6 = 1000)"), c("tied", 10, "deep, tied (\u03b6 = 10)"))
  for (s in dsel) {
    x <- ru[ru$solver == s[1] & (is.na(s[2]) | ru$K == as.numeric(s[2])), ]
    rows[[length(rows) + 1L]] <- data.frame(method = s[3], group = "deep", U0 = mean(x$err_U0_rel), Z0 = mean(x$err_Z0_rel),
                                            sec = mean(x$seconds), U0_lo = min(x$err_U0_rel), U0_hi = max(x$err_U0_rel),
                                            Z0_lo = min(x$err_Z0_rel), Z0_hi = max(x$err_Z0_rel), seeds = nrow(x))
  }
  out <- do.call(rbind, lapply(rows, function(d) { for (cc in c("U0_lo", "U0_hi", "Z0_lo", "Z0_hi", "seeds"))
    if (is.null(d[[cc]])) d[[cc]] <- NA; d }))
  out$case <- "seller, dynamic, lambda_U = 0.30, alpha = 0, r_f^- = 0.10"
  out
}
fig14 <- function() {
  m <- method_performance()
  slot <- c(1, 1, 1, 3, 5, 2, 4, 6, 4); pch <- c(16, 16, 16, 15, 18, 17, 1, 2, 2)
  fig_open("fig14_performance.pdf", height = PLOT_SETTINGS$height1 + 0.45)
  layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.17))
  for (q in c("U0", "Z0")) {
    e <- abs(m[[q]]) * 100
    lo <- abs(m[[paste0(q, "_lo")]]) * 100; hi <- abs(m[[paste0(q, "_hi")]]) * 100
    yl <- range(c(e, lo, hi), na.rm = TRUE) * c(0.5, 2)
    panel(range(m$sec) * c(0.5, 2), yl, "seconds per value", "|relative error| (%)", log = "xy",
          main = if (q == "U0") expression("(a) value"~XVA[0]^"+") else expression("(b) control"~"Z\u0306"[0]))
    i <- m$group == "quant"; lines(m$sec[i], e[i], col = PAL[1], lwd = 1.2)
    for (j in seq_len(nrow(m))) {
      if (!is.na(lo[j])) arrows(m$sec[j], min(lo[j], hi[j]), m$sec[j], max(lo[j], hi[j]), angle = 90, code = 3,
                                length = 0.02, col = PAL[slot[j]])
      points(m$sec[j], e[j], col = PAL[slot[j]], pch = pch[j], cex = 1)
    }
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", legend = c("quantisation, (25,50) to (100,200)", "Richardson", expression("PDE"~801%*%1500),
                              "deep, plain", "deep, tied \u03b6 = 1000", "deep, untied \u03b6 = 1000", "deep, tied \u03b6 = 10"),
         col = PAL[c(1, 3, 5, 2, 4, 6, 4)], pch = c(16, 15, 18, 17, 1, 2, 2), bty = "n", cex = 0.8, ncol = 4)
  fig_close()
  record_output("fig14", "fig14_performance.pdf", c("verification", "deep"), "5.8",
                "accuracy against run time by method, seller, dynamic lambda_U = 0.30, alpha = 0; deep: mean and range over 4 seeds")
}

# ---- Figure 15 (App. B): bias against lambda_U with seed dispersion -----------------
fig15 <- function() {
  r <- suite_data("deep"); ref <- r$reference
  m <- r$runs[r$runs$design == "main", ]
  sol <- c("plain", "tied", "untied"); off <- c(-0.06, 0, 0.06)
  fig_open("fig15_deep_bias.pdf", height = PLOT_SETTINGS$height1 + 0.35)
  layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), heights = c(1, 0.12))
  for (q in c("U0", "Z0")) {
    col <- if (q == "U0") "err_U0_rel" else "err_Z0_rel"
    panel(c(0.15, 1.85), range(100 * m[[col]]) * 1.05, expression(lambda[U]~"(years)"),
          "error against quantisation (%)",
          main = if (q == "U0") expression("(a) value"~XVA[0]^"+") else expression("(b) control"~"Z\u0306"[0]))
    disc <- if (q == "U0") ref$disc_U0 / ref$quant_U0 - 1 else ref$disc_Z0 / ref$quant_Z0 - 1
    lines(ref$lambda_U, 100 * disc, col = MUTED, lty = 3, lwd = 1.4)
    for (j in 1:3) {
      x <- m[m$solver == sol[j], ]
      points(x$lambda_U + off[j], 100 * x[[col]], col = PAL[j], pch = PCH[j], cex = 0.55)
      a <- aggregate(x[[col]] * 100, list(l = x$lambda_U), function(v) c(mean(v), sd(v)))
      errbars(a$l + off[j] + 0.025, a$x[, 1] - a$x[, 2], a$x[, 1] + a$x[, 2], j)
    }
  }
  par(mar = c(0, 0, 0, 0)); plot.new()
  legend("center", c("plain", "jump-tied", "jump-untied", "time step alone (n = 20)"), col = c(PAL[1:3], MUTED),
         pch = c(PCH[1:3], NA), lty = c(NA, NA, NA, 3), bty = "n", cex = 0.8, ncol = 4)
  fig_close()
  record_output("fig15", "fig15_deep_bias.pdf", "deep", "App. B",
                "deep-solver errors against lambda_U: seeds (points), mean +/- sd (bars), n = 20 time step alone (dotted); penalty weight zeta = 1000")
}

# ---- Figure 16 (App. B): penalty sweep and final training loss ----------------------
fig16 <- function() {
  r <- suite_data("deep"); ref <- r$reference
  k <- r$runs[r$runs$lambda_U == 0.3 & r$runs$solver %in% c("tied", "untied"), ]
  fx <- c(tied = 0.92, untied = 1.08); sl <- c(tied = 2, untied = 3)
  fig_open("fig16_penalty.pdf")
  panel(c(0.6, 2e4), range(100 * k$err_U0_rel) * 1.05, expression("penalty weight"~zeta), "error in XVA (%)",
        log = "x", main = expression("(a) error,"~lambda[U]=="0.30"))
  abline(h = 100 * (ref$disc_U0 / ref$quant_U0 - 1)[ref$lambda_U == 0.3], col = MUTED, lty = 3, lwd = 1.4)
  for (s in names(fx)) { x <- k[k$solver == s, ]; points(x$K * fx[[s]], 100 * x$err_U0_rel, col = PAL[sl[[s]]], pch = PCH[sl[[s]]], cex = 0.55)
    a <- aggregate(100 * err_U0_rel ~ K, x, mean); series(a$K * fx[[s]], a[, 2], sl[[s]]) }
  panel(c(0.6, 2e4), range(k$final_mse) * c(0.8, 1.25), expression("penalty weight"~zeta),
        "final training loss (mean of last 100)", log = "xy", main = "(b) final training loss")
  for (s in names(fx)) { x <- k[k$solver == s, ]; points(x$K * fx[[s]], x$final_mse, col = PAL[sl[[s]]], pch = PCH[sl[[s]]], cex = 0.55)
    a <- aggregate(final_mse ~ K, x, mean); series(a$K * fx[[s]], a[, 2], sl[[s]]) }
  legend("right", c("jump-tied", "jump-untied", "time step alone (a)"), col = c(PAL[2:3], MUTED), lty = c(LTY[2:3], 3),
         pch = c(PCH[2:3], NA), lwd = 1.4, bty = "n", cex = 0.78, seg.len = 2.4)
  fig_close()
  record_output("fig16", "fig16_penalty.pdf", "deep", "App. B",
                "penalty sweep: error and final training loss against the penalty weight zeta for the jump solvers, lambda_U = 0.30, n = 20 time steps, four seeds each (loss histories were not stored)")
}

# ------------------------------------------------------------------------------

##  TABLES

# ------------------------------------------------------------------------------

# ---- Table 1 (5.1): parameters and solver settings ---------------------------------
tab01 <- function() {
  p <- base_params(); vs <- suite_data("verification")$meta$settings; ds <- suite_data("deep")$meta$settings
  rs <- suite_data("replication")$meta$settings
  rw <- function(a, b, c) data.frame(q = a, s = b, v = c, stringsAsFactors = FALSE)
  df <- rbind(
    rw("Initial stock price, strike", "$S_0$, $K$", sprintf("%g, %g", p$S0, p$K)),
    rw("Maturity (years)", "$T$", p$Tmat),
    rw("Volatility", "$\\sigma$", p$sigma),
    rw("Discount rate", "$r_D$", p$rD),
    rw("Repo rates", "$r_r^+$, $r_r^-$", sprintf("%g, %g", p$rr_p, p$rr_m)),
    rw("Collateral rates", "$r_c^+$, $r_c^-$", sprintf("%g, %g", p$rc_p, p$rc_m)),
    rw("Funding rates", "$r_f^+$, $r_f^-$", sprintf("%g; 0.05, 0.075, 0.10", p$rf_p)),
    rw("Bond return rates", "$\\mu_I$, $\\mu_C$", sprintf("%g, %g", p$muI, p$muC)),
    rw("Default intensities", "$h_I$, $h_C$", sprintf("%g, %g", p$hI, p$hC)),
    rw("Loss rates", "$L_I$, $L_C$", sprintf("%g, %g", p$LI, p$LC)),
    rw("Mean holding times (years)", "$\\lambda_U$, $\\lambda_V$", sprintf("%g, %g", p$mean_U, p$mean_V)),
    rw("Switching rates", "$a_U=1/\\lambda_U$, $a_V=1/\\lambda_V$", sprintf("%.4f, %.4f", p$aU, p$aV)),
    rw("Collateralisation level", "$\\alpha$", "0 to 1"),
    rw("MIDRULE", "", ""),
    rw("Quantisation: steps, grid size", "$(n,N)$", "(25,50), (50,100), (100,200)"),
    rw("PDE grid (reference; standard)", "$n_x\\times n_t$", sprintf("%d $\\times$ %d; %d $\\times$ %d",
       vs$pde_ref[1], vs$pde_ref[2], vs$pde_std[1], vs$pde_std[2])),
    rw("Closed form: quadrature time nodes", "$n_s$", paste(vs$cf_ns, collapse = ", ")),
    rw("Deep solvers: time steps; iterations", "$n$", sprintf("%d; %d", ds$N, ds$iters)),
    rw("Deep solvers: batch (plain; jump)", "$M$", sprintf("%d; %d", ds$M_plain, ds$M_jump)),
    rw("Deep solvers: seeds; penalty weight", "$\\zeta$", sprintf("%d seeds; %s", length(ds$seeds), paste(ds$K_sweep, collapse = ", "))),
    rw("Replication simulator: paths, fine steps", "", sprintf("%d; %d", rs$sim_M, rs$sim_nt)))
  tex_table("tab01_parameters.tex", df, "lll", c("Quantity", "Symbol", "Value"),
            note = "Model parameters (Pang and Sturm benchmark set) and solver settings, read from the saved suite metadata.")
  record_output("tab01", "tab01_parameters.tex", c("verification", "deep", "replication"), "5.1", "parameters and solver settings")
}

# ---- Table 2 (5.2): errors in value and control by regime and grid -----------------
tab02 <- function() {
  sm <- suite_data("verification")$summary; sm <- sm[sm$side == "+", ]
  key <- unique(sm[, c("regime", "mean_U", "alpha")])
  key <- key[order(match(ifelse(key$regime == "dynamic", sprintf("dynamic lamU=%.2f", key$mean_U), key$regime),
                         names(REGIME_SLOT)), key$alpha), ]
  df <- do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
    k <- key[i, ]; x <- sm[sm$regime == k$regime & (is.na(k$mean_U) | sm$mean_U == k$mean_U) & sm$alpha == k$alpha, ]
    u <- x[x$quantity == "U0", ]; z <- x[x$quantity == "Z0", ]
    data.frame(r = regime_tex(ifelse(k$regime == "dynamic", sprintf("dynamic lamU=%.2f", k$mean_U), k$regime)),
               a = fmt(k$alpha, 1),
               u1 = fmt(u$relerr_50_100, 3, TRUE, TRUE), u2 = fmt(u$relerr_100_200, 3, TRUE, TRUE),
               ur = fmt(u$relerr_rich_100_200, 4, TRUE, TRUE),
               z1 = fmt(z$relerr_50_100, 2, TRUE, TRUE), z2 = fmt(z$relerr_100_200, 2, TRUE, TRUE),
               zr = fmt(z$relerr_rich_100_200, 3, TRUE, TRUE), o = fmt(u$observed_order, 2), stringsAsFactors = FALSE)
  }))
  tex_table("tab02_verification.tex", df, "llrrrrrrr",
            c("Regime", "$\\alpha$", "(50,100)", "(100,200)", "Rich.", "(50,100)", "(100,200)", "Rich.", "Order"),
            header2 = "& & \\multicolumn{3}{c}{$\\mathrm{XVA}_0$ error (\\%)} & \\multicolumn{3}{c}{$\\breve Z_0$ error (\\%)} & \\\\ \\cmidrule(lr){3-5}\\cmidrule(lr){6-8}",
            note = "Seller, r_f^- = 0.10; errors against the closed form where its branch holds, otherwise the PDE on 1601 x 3000; Order: observed order of the value, p = log2[(XVA_0^(25) - XVA_0^(50)) / (XVA_0^(50) - XVA_0^(100))] over the levels (25,50), (50,100) and (100,200), the same three levels as Figure 1; it needs no reference value.")
  record_output("tab02", "tab02_verification.tex", "verification", "5.2", "errors in value and control by regime and grid")
}

# ---- Table 3 (5.3): seller's XVA by regime, alpha, r_f^- -----------------------------
tab03 <- function() {
  r <- suite_data("ps"); f <- rbind(r$figure4, r$figure5[r$figure5$convention == "mean", ])
  f <- f[f$alpha %in% c(0, 0.5, 1), ]
  g <- suite_data("bands")$grid
  bands_regime <- c("frozen0" = "frozen0", "frozen1" = "frozen1", "0.3" = "dynamic lamU=0.30", "1.39" = "dynamic lamU=1.39")
  rows <- list(list("frozen0", NA, "normal"), list("frozen1", NA, "crisis"),
               list("dynamic", 0.3, "dynamic, $\\lambda_U=0.30$"), list("dynamic", 1, "dynamic, $\\lambda_U=1.00$"),
               list("dynamic", 1.39, "dynamic, $\\lambda_U=1.39$ (calibrated)"), list("dynamic", 1.7, "dynamic, $\\lambda_U=1.70$"))
  ## The Pang-Sturm grid first; the bands suite (same solver) for lambda_U = 1.39
  ## and for the cells the grid lacks.
  df <- do.call(rbind, lapply(rows, function(rw) {
    vals <- unlist(lapply(RF_LEVELS, function(rf) sapply(c(0, 0.5, 1), function(a) {
      x <- f[f$regime == rw[[1]] & abs(f$rf_m - rf) < 1e-9 & f$alpha == a &
             (is.na(rw[[2]]) | (!is.na(f$lambda_U) & f$lambda_U == rw[[2]])), ]
      if (nrow(x)) return(x$quantisation_pct_Vhat0[1])
      bk <- bands_regime[if (is.na(rw[[2]])) rw[[1]] else as.character(rw[[2]])]
      y <- if (is.na(bk)) NULL else g[g$regime == bk & abs(g$rf_m - rf) < 1e-9 & abs(g$alpha - a) < 1e-9, ]
      if (NROW(y)) y$XVA_plus_pct_Vhat0[1] else NA })))
    data.frame(t(c(rw[[3]], fmt(vals, 2))), stringsAsFactors = FALSE)
  }))
  tex_table("tab03_xva.tex", df, "lrrrrrrrrr",
            c("$\\alpha$", rep(c("0", "0.5", "1"), 3)),
            header2 = paste0("& \\multicolumn{3}{c}{$r_f^-=0.05$} & \\multicolumn{3}{c}{$r_f^-=0.075$} & \\multicolumn{3}{c}{$r_f^-=0.10$} \\\\ ",
                             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}"),
            note = "Seller's XVA_0 in per cent of Vhat_0 (quantisation, Richardson); the values plotted in Figures 2 and 3; lambda_U = 1.39 and lambda_U = 0.30 at r_f^- = 0.075 from the bands suite; -- not computed.")
  record_output("tab03", "tab03_xva.tex", c("ps", "bands"), "5.3", "seller's XVA by regime, alpha and r_f^-")
}

# ---- Table 4 (5.3): model risk, frozen against dynamic ---------------------------------
tab04 <- function() {
  ## Exactly at the calibrated point (1.39, 0.99): the dynamic value from the
  ## regime suite's spot check there, the frozen values from the heat map
  ## (constant over it), and the crisis share pbar in closed form.
  r <- suite_data("regime"); p <- with_params(base_params(), alpha = 0.5, rf_m = 0.10)
  sc <- pct_Vhat0(); h <- r$heatmap[r$heatmap$side == "+" & r$heatmap$rf_m == 0.10, ]
  calm <- h$frozen_calm[1]; crisis <- h$frozen_crisis[1]
  df <- do.call(rbind, lapply(0:1, function(b) {
    dyn <- r$pde_check$quant[r$pde_check$side == "+" & r$pde_check$beta0 == b &
                             abs(r$pde_check$lambda_U - 1.39) < 1e-9 & abs(r$pde_check$lambda_V - 0.99) < 1e-9]
    pb  <- crisis_share(p, b); mix <- (1 - pb) * calm + pb * crisis
    init <- if (b == 0) calm else crisis
    data.frame(b0 = if (b == 0) "normal, $\\beta_0=0$" else "crisis, $\\beta_0=1$", dyn = fmt(sc * dyn, 2), pb = fmt(pb, 3),
               fi = fmt(sc * (init - dyn), 2, sign = TRUE), fc = fmt(sc * (crisis - dyn), 2, sign = TRUE),
               mx = fmt(sc * (mix - dyn), 3, sign = TRUE), stringsAsFactors = FALSE)
  }))
  tex_table("tab04_model_risk.tex", df, "lrrrrr",
            c("Initial regime", "Dynamic $\\mathrm{XVA}_0$", "$\\bar p$", "Frozen initial", "Frozen crisis", "Mixture"),
            header2 = "& & & \\multicolumn{3}{c}{Shortcut minus dynamic} \\\\ \\cmidrule(lr){4-6}",
            note = paste("Suggested caption: Model risk of three shortcuts at the calibrated holding times (lambda_U, lambda_V) = (1.39, 0.99),",
                         "seller, alpha = 0.5, r_f^- = 0.10, in per cent of Vhat_0. Frozen initial: the regime frozen at beta_0;",
                         "frozen crisis: the crisis frozen throughout; mixture: (1 - pbar) XVA_calm + pbar XVA_crisis with",
                         "pbar = (1/T) int_0^T P(beta_t = 1 | beta_0) dt. Entries are shortcut minus dynamic. For beta_0 = 1 the",
                         "two frozen shortcuts coincide. Ranges over the 15 x 15 map: x_model_risk_ranges.tex."))
  record_output("tab04", "tab04_model_risk.tex", "regime", "5.3", "model risk of three shortcuts at the calibrated point")
}

# The ranges over the heat map, formerly in Table 4 (appendix or a note).
x_model_risk_ranges <- function() {
  h <- suite_data("regime")$heatmap; m <- h[h$side == "+" & h$rf_m == 0.10, ]
  rng <- function(v, d = 2) sprintf("[%s, %s]", fmt(min(v), d, sign = TRUE), fmt(max(v), d, sign = TRUE))
  df <- do.call(rbind, lapply(0:1, function(b) { x <- m[m$beta0 == b, ]
    data.frame(b0 = if (b == 0) "normal, $\\beta_0=0$" else "crisis, $\\beta_0=1$",
               fi = rng(x$risk_frozen_initial_pct_Vhat0), fc = rng(x$risk_frozen_crisis_pct_Vhat0),
               mx = rng(x$risk_mixture_pct_Vhat0, 3), stringsAsFactors = FALSE) }))
  tex_table("x_model_risk_ranges.tex", df, "lrrr", c("Initial regime", "Frozen initial", "Frozen crisis", "Mixture"),
            note = "Shortcut minus dynamic XVA in per cent of Vhat_0, range over the 15 x 15 grid of (lambda_U, lambda_V) in [0.25, 3]^2; seller, alpha = 0.5, r_f^- = 0.10.")
  record_output("x_model_risk_ranges", "x_model_risk_ranges.tex", "regime", "not in the body", "model-risk ranges over the heat map")
}

# ---- Table 5 (5.4): seller, buyer and band width ----------------------------------
tab05 <- function() {
  g <- suite_data("bands")$grid; g <- g[g$alpha %in% c(0, 0.5, 1), ]
  df <- do.call(rbind, lapply(names(REGIME_SLOT), function(rg) do.call(rbind, lapply(RF_LEVELS, function(rf) {
    x <- g[g$regime == rg & abs(g$rf_m - rf) < 1e-9, ]; x <- x[order(x$alpha), ]
    cells <- as.vector(rbind(fmt(x$XVA_plus_pct_Vhat0, 2), fmt(x$XVA_minus_pct_Vhat0, 2), fmt(x$width_pct_Vhat0, 2, sign = TRUE)))
    data.frame(t(c(if (rf == RF_LEVELS[1]) regime_tex(rg) else "", fmt(rf, 3), cells)), stringsAsFactors = FALSE)
  }))))
  tex_table("tab05_band.tex", df, "llrrrrrrrrr",
            c("Regime", "$r_f^-$", rep(c("Seller", "Buyer", "Width"), 3)),
            header2 = paste0("& & \\multicolumn{3}{c}{$\\alpha=0$} & \\multicolumn{3}{c}{$\\alpha=0.5$} & \\multicolumn{3}{c}{$\\alpha=1$} \\\\ ",
                             "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}"),
            note = "XVA+ (seller), XVA- (buyer, admissible crisis convention) and the band width in per cent of Vhat_0; a negative width is an inverted band.")
  record_output("tab05", "tab05_band.tex", "bands", "5.4", "seller, buyer and band width")
}

# ---- Table 6 (5.5): six accounts at inception, with balance check -------------------
tab06 <- function() {
  t0 <- suite_data("replication")$t0; h <- t0[t0$alpha == 0.5 & t0$rf_m == 0.10, ]
  rows <- list(c("+", "frozen0", "Seller, normal"), c("+", "frozen1", "Seller, crisis"),
               c("+", "dynamic", "Seller, dynamic"), c("-", "frozen0", "Buyer, normal\\textsuperscript{b}"),
               c("-", "frozen1", "Buyer, crisis\\textsuperscript{b,c}"))
  df <- do.call(rbind, lapply(rows, function(rw) {
    x <- h[h$side == rw[1] & h$regime == rw[2], ]; v <- setNames(x$value, x$account)
    data.frame(case = rw[3], stock = fmt(v[["stock"]], 4), funding = fmt(v[["funding"]], 4),
               bondI = fmt(v[["investor_bond"]], 4), bondC = fmt(v[["counterparty_bond"]], 4),
               repo = fmt(v[["repo"]], 4), coll = fmt(v[["collateral"]], 4),
               xva = fmt(x$XVA0_pct_Vhat0[1], 2), bal = fmte(max(abs(x$residual)), 1), stringsAsFactors = FALSE)
  }))
  tex_table("tab06_accounts.tex", df, "lrrrrrrrr",
            header2 = "& \\multicolumn{6}{c}{Account value (per unit notional)} & & \\\\ \\cmidrule(lr){2-7}",
            c("", "Stock", "Funding", "$P^I$", "$P^C$", "Repo", "Collateral", "$\\mathrm{XVA}_0$ (\\%)", "Balance"),
            note = "Account values per unit notional (S_0 = K = 1, so the stock account also gives the shares held); XVA_0 in per cent of Vhat_0; Balance: |sum of the account values - wealth|. alpha = 0.5, r_f^- = 0.10. The buyer in the crisis holds no stock at any of the 4,901 crisis nodes of the (50,100) tree (0 inadmissible; likewise in the other eight (alpha, r_f^-) cases). (b) holder's position by the reflection (4.17): CVA loss, no DVA. (c) admissible crisis convention: no short stock in the crisis, residual crisis delta valued under Q.")
  record_output("tab06", "tab06_accounts.tex", "replication", "5.5", "six accounts at inception with balance check")
}

# ---- Table 7 (5.6): Shapley values and spread across orderings ------------------------
tab07 <- function() {
  v <- suite_data("shapley")$values; v <- v[v$version == "main" & v$alpha == 0.5 & v$rf_m == 0.10, ]
  V0 <- 100 / bs_value(0, base_params()$S0, base_params())
  rows <- list(c("+", "dynamic", "1.39", "Seller, dynamic, $\\lambda_U=1.39$"), c("+", "dynamic", "0.3", "Seller, dynamic, $\\lambda_U=0.30$"),
               c("+", "frozen1", NA, "Seller, crisis"), c("-", "dynamic", "1.39", "Buyer\\textsuperscript{b}"))
  feats <- c("funding_spread", "dva", "cva", "regime")
  df <- do.call(rbind, lapply(rows, function(rw) {
    x <- v[v$side == rw[1] & v$regime == rw[2] & (is.na(rw[3]) | abs(v$mean_U - as.numeric(rw[3])) < 1e-9), ]
    y <- x[match(feats, x$feature), ]
    rbind(data.frame(t(c(rw[4], fmt(x$pct_Vhat0[x$feature == "baseline"], 2), fmt(y$pct_Vhat0, 2, sign = TRUE),
                         fmt(x$XVA0_pct_Vhat0[1], 2))), stringsAsFactors = FALSE),
          data.frame(t(c("", "", sprintf("{\\scriptsize[%s, %s]}", fmt(y$min_marginal * V0, 2), fmt(y$max_marginal * V0, 2)), "")),
                     stringsAsFactors = FALSE))
  }))
  tex_table("tab07_shapley.tex", df, "lrrrrrr",
            c("", "Baseline", "Spread", "DVA", "CVA", "Regime", "$\\mathrm{XVA}_0$"),
            note = "Shapley values over 24 orderings in per cent of Vhat_0, with the range of the marginal contribution across orderings in brackets; alpha = 0.5, r_f^- = 0.10. (b) buyer under the admissible crisis convention (no short stock in the crisis).")
  record_output("tab07", "tab07_shapley.tex", "shapley", "5.6", "Shapley values and spread across orderings")
}

# ---- Table 8 (5.7): AAD sensitivities and agreement with bumps ------------------------
tab08 <- function() {
  r <- suite_data("aad"); c2 <- r$continuous; d1 <- r$discrete
  par_lab <- c(S0 = "$S_0$", sigma = "$\\sigma$", rf_m = "$r_f^-$", hI = "$h_I$", hC = "$h_C$",
               alpha = "$\\alpha$", LI = "$L_I$", LC = "$L_C$")
  sc <- c("seller frozen0", "seller frozen1", "seller dynamic", "buyer frozen0", "buyer frozen1")
  pars <- unique(c2$parameter)
  df <- data.frame(p = par_lab[pars], stringsAsFactors = FALSE)
  for (s0 in sc) { y <- c2[c2$scenario == s0, ]; y <- y[match(pars, y$parameter), ]
    df[[s0]] <- ifelse(abs(y$adjoint_rich) < 1e-12, "0", fmt(y$adjoint_rich, 5, sign = TRUE)) }
  df$disc <- vapply(pars, function(pp) fmte(max(abs(d1$rel_diff[d1$parameter == pp])), 1), "")
  df$cont <- vapply(pars, function(pp) { y <- c2[c2$parameter == pp & abs(c2$pde_bump) > 1e-9, ]
    fmte(max(abs(y$relerr_rich)), 1) }, "")
  tex_table("tab08_aad.tex", df, "lrrrrrrr",
            c("", "Normal", "Crisis", "Dynamic", "Normal\\textsuperscript{b}", "Crisis\\textsuperscript{b,c}", "Same scheme", "Continuous"),
            header2 = paste0("& \\multicolumn{3}{c}{Seller} & \\multicolumn{2}{c}{Buyer} & \\multicolumn{2}{c}{Largest difference to bumps} \\\\ ",
                             "\\cmidrule(lr){2-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}"),
            note = "Adjoint sensitivities (Richardson); Same scheme: largest relative difference to central bumps of the same scheme over three levels and six scenarios; Continuous: largest relative difference to PDE or closed-form bumps after Richardson. alpha = 0.5, r_f^- = 0.10. (b) holder's position by the reflection (4.17): CVA loss, no DVA. (c) admissible crisis convention: no short stock in the crisis, residual crisis delta valued under Q.")
  record_output("tab08", "tab08_aad.tex", "aad", "5.7", "AAD sensitivities and agreement with bumps")
}

# ---- Table 9 (5.8): method performance summary ------------------------------------------
tab09 <- function() {
  m <- method_performance()
  rng <- function(lo, hi) ifelse(is.na(lo), "", sprintf(" {\\scriptsize[%s, %s]}", fmt(lo, 2, TRUE, TRUE), fmt(hi, 2, TRUE, TRUE)))
  df <- data.frame(meth = m$method,
                   u = paste0(fmt(m$U0, 3, TRUE, TRUE), rng(m$U0_lo, m$U0_hi)),
                   z = paste0(fmt(m$Z0, 3, TRUE, TRUE), rng(m$Z0_lo, m$Z0_hi)),
                   s = ifelse(m$sec < 1, fmt(m$sec, 3), fmt(m$sec, 1)),
                   n = ifelse(is.na(m$seeds), "", m$seeds), stringsAsFactors = FALSE)
  df$meth <- sub("801 x 1500", "$801\\\\times1500$", df$meth)
  df$meth <- gsub("\u03b6", "$\\zeta$", df$meth, fixed = TRUE)
  tex_table("tab09_performance.tex", df, "lrrrr",
            c("Method", "$\\mathrm{XVA}_0$ error (\\%)", "$\\breve Z_0$ error (\\%)", "Seconds", "Seeds"),
            note = "Seller, dynamic regime, lambda_U = 0.30, alpha = 0, r_f^- = 0.10; errors against the PDE on 1601 x 3000 (quantisation, PDE) or Richardson quantisation (deep); deep: mean and [range] over the seeds; quantisation seconds include the tree build. The closed form does not apply to the dynamic regime.")
  record_output("tab09", "tab09_performance.tex", c("verification", "deep"), "5.8", "method performance summary")
}

# ---- Table 10 (App. C): dollar funding positions ------------------------------------------
tab10 <- function() {
  a <- suite_data("bands")$accounts_t0; a <- a[a$account == "funding" & a$alpha %in% c(0, 0.5, 1), ]
  df <- do.call(rbind, lapply(names(REGIME_SLOT), function(rg) do.call(rbind, lapply(c("+", "-"), function(sd) {
    vals <- unlist(lapply(c(0.05, 0.10), function(rf) sapply(c(0, 0.5, 1), function(al)
      a$value[a$regime == rg & a$side == sd & abs(a$rf_m - rf) < 1e-9 & a$alpha == al])))
    data.frame(t(c(if (sd == "+") regime_tex(rg) else "", if (sd == "+") "seller" else "buyer\\textsuperscript{b}",
                   fmt(vals, 4))), stringsAsFactors = FALSE)
  }))))
  tex_table("tab10_funding.tex", df, "llrrrrrr", c("Regime", "Side", paste0("$\\alpha=", c("0", "0.5", "1"), "$"), paste0("$\\alpha=", c("0", "0.5", "1"), "$")),
            header2 = paste0("& & \\multicolumn{3}{c}{$r_f^-=0.05$} & \\multicolumn{3}{c}{$r_f^-=0.10$} \\\\ ",
                             "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}"),
            note = "Funding account value at t = 0 per unit notional (S_0 = K = 1; positive: lending to the treasury); the buyer's hedge portfolio is the reflected one. (b) buyer under the admissible crisis convention (no short stock in the crisis).")
  record_output("tab10", "tab10_funding.tex", "bands", "App. C", "dollar funding positions, seller and buyer")
}

# ---- Table 11 (App. C): replication error against the rebalancing step ---------------
tab11 <- function() {
  r <- suite_data("replication"); d <- r$sim[r$sim$source == "pde", ]
  rk <- ifelse(d$regime == "dynamic", sprintf("dynamic lamU=%.2f", d$mean_U), d$regime)
  ms <- sort(unique(d$rebalances))
  df <- do.call(rbind, lapply(unique(rk), function(g) { x <- d[rk == g, ]; x <- x[order(x$rebalances), ]
    data.frame(t(c(regime_tex(g), fmt(x$rmse_pct_Vhat0, 2),
                   if (x$regime[1] == "dynamic") fmt(x$rmse_adj_pct_Vhat0[nrow(x)], 2) else "--",
                   fmt(r$decay$rmse_slope[r$decay$case == x$case[1]], 3))), stringsAsFactors = FALSE) }))
  tex_table("tab11_replication_error.tex", df, paste0("l", strrep("r", length(ms) + 2)),
            c("Regime", ms, "Adjusted", "Slope"),
            note = "RMSE of the seller's replication error in per cent of Vhat_0 by number of rebalancing dates; alpha = 0.5, r_f^- = 0.10; 20,000 paths on 1,600 fine steps, hedge ratios from the PDE (801 space nodes); Adjusted: regime-switch martingale added back, at the finest step; Slope: of log RMSE on log step.")
  record_output("tab11", "tab11_replication_error.tex", "replication", "App. C", "replication error against the rebalancing step")
}

# ---- Table 12 (App. B): deep solvers by seed and lambda_U ------------------------------
tab12 <- function() {
  ru <- suite_data("deep")$runs; ru <- ru[ru$design == "main", ]
  df <- do.call(rbind, lapply(sort(unique(ru$lambda_U)), function(l) do.call(rbind, lapply(c("plain", "tied", "untied"), function(s) {
    x <- ru[ru$lambda_U == l & ru$solver == s, ]; x <- x[order(x$seed), ]
    data.frame(t(c(if (s == "plain") lam_lab(l) else "", s, fmt(x$err_U0_rel, 2, TRUE, TRUE),
                   fmt(mean(x$err_U0_rel), 2, TRUE, TRUE), fmt(sd(x$err_U0_rel), 2, TRUE),
                   fmt(mean(x$err_Z0_rel), 1, TRUE, TRUE), fmt(mean(x$seconds), 0))), stringsAsFactors = FALSE)
  }))))
  tex_table("tab12_deep_seeds.tex", df, "llrrrrrrrr",
            c("$\\lambda_U$", "Solver", "Seed 1", "Seed 2", "Seed 3", "Seed 4", "Mean", "Sd", "$\\breve Z_0$ mean", "Seconds"),
            note = "Errors in per cent against Richardson quantisation; penalty weight zeta = 1000 for the jump solvers; n = 20 time steps; alpha = 0, r_f^- = 0.10.")
  record_output("tab12", "tab12_deep_seeds.tex", "deep", "App. B", "deep solvers by seed and lambda_U")
}

# ---- Table 13 (App. B): penalty-weight sweep -------------------------------------------------
tab13 <- function() {
  ru <- suite_data("deep")$runs; ru <- ru[ru$lambda_U == 0.3 & ru$solver %in% c("tied", "untied"), ]
  df <- do.call(rbind, lapply(sort(unique(ru$K)), function(K) {
    cells <- unlist(lapply(c("tied", "untied"), function(s) { x <- ru[ru$solver == s & ru$K == K, ]
      c(sprintf("%s (%s)", fmt(mean(x$err_U0_rel), 2, TRUE, TRUE), fmt(sd(x$err_U0_rel), 2, TRUE)),
        fmt(mean(x$err_U0_rel_disc), 2, TRUE, TRUE), fmte(mean(x$final_mse), 2)) }))
    data.frame(t(c(formatC(K, format = "d", big.mark = ","), cells)), stringsAsFactors = FALSE)
  }))
  tex_table("tab13_penalty.tex", df, "lrrrrrr", c("$\\zeta$", rep(c("Error (sd)", "Learning", "Final loss"), 2)),
            header2 = "& \\multicolumn{3}{c}{jump-tied} & \\multicolumn{3}{c}{jump-untied} \\\\ \\cmidrule(lr){2-4}\\cmidrule(lr){5-7}",
            note = "lambda_U = 0.30, four seeds; errors in per cent against Richardson quantisation; Learning: against the same scheme at n = 20 time steps; Final loss: mean squared terminal mismatch over the last 100 iterations.")
  record_output("tab13", "tab13_penalty.tex", "deep", "App. B", "penalty-weight sweep")
}

# ---- Table 14 (App. C): AAD against PDE bumps by regime; frozen against full ---------
tab14 <- function() {
  c2 <- suite_data("aad")$continuous
  par_lab <- c(S0 = "$S_0$", sigma = "$\\sigma$", rf_m = "$r_f^-$", hI = "$h_I$", hC = "$h_C$",
               alpha = "$\\alpha$", LI = "$L_I$", LC = "$L_C$")
  sc <- c("seller frozen0", "seller frozen1", "seller dynamic", "buyer frozen1")
  pars <- unique(c2$parameter)
  cellf <- function(y, col) ifelse(abs(y$pde_bump) < 1e-9, "--", fmt(y[[col]], 2, TRUE, TRUE))
  df <- data.frame(p = par_lab[pars], stringsAsFactors = FALSE)
  for (s0 in sc) { y <- c2[c2$scenario == s0, ]; y <- y[match(pars, y$parameter), ]
    df[[paste(s0, "l")]] <- cellf(y, "relerr_50_100"); df[[paste(s0, "r")]] <- cellf(y, "relerr_rich") }
  fr <- data.frame(p = "$\\sigma$, frozen grid", stringsAsFactors = FALSE)
  for (s0 in sc) { y <- c2[c2$scenario == s0 & c2$parameter == "sigma", ]
    ref <- ifelse(is.na(y$closed_form_bump), y$pde_bump, y$closed_form_bump)
    fr[[paste(s0, "l")]] <- fmt(y$adjoint_frozen_grid_50_100 / ref - 1, 1, TRUE, TRUE); fr[[paste(s0, "r")]] <- "" }
  names(fr) <- names(df)
  df <- rbind(df, data.frame(t(c("MIDRULE", rep("", ncol(df) - 1))), stringsAsFactors = FALSE) |> setNames(names(df)), fr)
  tex_table("tab14_aad_pde.tex", df, "lrrrrrrrr", c("", rep(c("(50,100)", "Rich."), 4)),
            header2 = paste0("& \\multicolumn{2}{c}{Seller, normal} & \\multicolumn{2}{c}{Seller, crisis} & ",
                             "\\multicolumn{2}{c}{Seller, dynamic} & \\multicolumn{2}{c}{Buyer, crisis\\textsuperscript{b,c}} \\\\ ",
                             "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}"),
            note = "Relative error in per cent of the adjoint against central bumps of the PDE (closed form for the seller in frozen regimes); -- the sensitivity vanishes on its branch; last row: vega with the grid frozen at (50,100). alpha = 0.5, r_f^- = 0.10. The buyer in the normal regime mirrors the seller's normal column (h_I <-> h_C, L_I <-> L_C) and is reported in Table 8; the buyer's crisis column tests the adjoint across the kink of the projection. (b) holder's position by the reflection (4.17): CVA loss, no DVA. (c) admissible crisis convention: no short stock in the crisis, residual crisis delta valued under Q.")
  record_output("tab14", "tab14_aad_pde.tex", "aad", "App. C", "AAD against PDE bumps by regime; frozen against full")
}

# ---- Table 15 (App. D): run-time record -------------------------------------------------
tab15 <- function() {
  mf <- read.csv(file.path(OUTPUT_DATA, "manifest.csv"), stringsAsFactors = FALSE)
  cmd <- c(verification = "verify", replication = "replicate", shapley = "shapley", bands = "bands",
           regime = "regime", deep = "deep", aad = "aad", ps = "valuation", buyerconv = "closeout")
  lab <- c(verification = "verification", replication = "replication", shapley = "Shapley", bands = "bands",
           regime = "regime", deep = "deep solvers", aad = "AAD", ps = "valuation", buyerconv = "close-out comparison")
  mod <- c(verification = 17, replication = 18, shapley = 19, bands = 20, regime = 21, deep = 22,
           aad = 23, ps = 24, buyerconv = 25)
  train <- c(deep = sum(suite_data("deep")$runs$seconds), ps = sum(suite_data("ps")$deep_runs$seconds))
  mf <- mf[mf$suite %in% names(cmd), ]
  mf <- mf[order(mod[mf$suite]), ]
  df <- data.frame(s = lab[mf$suite], m = mod[mf$suite], c = paste0("\\texttt{", cmd[mf$suite], "}"),
                   w = fmt(mf$seconds / 60, 1),
                   t = ifelse(mf$suite %in% names(train), fmt(train[mf$suite] / 60, 0), "--"),
                   k = mf$cores, sd = mf$seed, d = substr(mf$created, 1, 10), stringsAsFactors = FALSE)
  tex_table("tab15_runtime.tex", df, "lrlrrrrl",
            c("Suite", "Module", "Command", "Minutes", "Training core-min", "Cores", "Seed", "Run on"),
            note = "Wall time of the last run of each suite (R 4.3.3); deep solvers and valuation reused their per-run training caches on the last run, so their training time is listed separately as the sum of single-run times. valuation: the XVA grid over alpha, r_f^- and lambda_U behind Figures 2 and 3; close-out comparison: the buyer under both crisis conventions (Table 16).")
  record_output("tab15", "tab15_runtime.tex", c("verification", "deep", "ps", "buyerconv"), "App. D", "run-time record")
}

# ---- Table 16 (App. C): the buyer under the two crisis conventions ---------------
tab16 <- function() {
  r <- suite_data("buyerconv"); v <- r$values; v <- v[v$alpha == 0.5 & abs(v$rf_m - 0.10) < 1e-9, ]
  sc <- pct_Vhat0(); ad <- r$admissibility
  df <- do.call(rbind, lapply(names(REGIME_SLOT), function(rg) {
    x <- function(cv, col) v[v$convention == cv & v$regime == rg, col]
    n <- function(cv) ad$inadmissible_nodes[ad$convention == cv & ad$regime == rg]
    cn <- ad$crisis_nodes[ad$convention == "banned" & ad$regime == rg]
    data.frame(rg = regime_tex(rg), sel = fmt(sc * x("banned", "XVA_plus"), 2),
               be = fmt(x("extended", "XVA_minus_pct_Vhat0"), 2), bb = fmt(x("banned", "XVA_minus_pct_Vhat0"), 2),
               we = fmt(x("extended", "width_pct_Vhat0"), 2, sign = TRUE), wb = fmt(x("banned", "width_pct_Vhat0"), 2, sign = TRUE),
               ie = if (cn == 0) "--" else formatC(n("extended"), format = "d", big.mark = ","),
               ib = if (cn == 0) "--" else formatC(n("banned"), format = "d", big.mark = ","), stringsAsFactors = FALSE)
  }))
  ca <- r$carry[r$carry$convention == "banned", ]; inv <- tapply(r$values$inverted, r$values$convention, sum)
  tex_table("tab16_closeout.tex", df, "lrrrrrrr",
            c("Regime", "$\\mathrm{XVA}^+_0$", "Ext.", "Adm.", "Ext.", "Adm.", "Ext.", "Adm."),
            header2 = paste0("& Seller & \\multicolumn{2}{c}{Buyer $\\mathrm{XVA}^-_0$} & \\multicolumn{2}{c}{Band width} & ",
                             "\\multicolumn{2}{c}{Inadmissible crisis nodes} \\\\ \\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}"),
            note = sprintf(paste("Suggested caption: The buyer under the two crisis conventions, alpha = 0.5, r_f^- = 0.10, in per cent of Vhat_0.",
                                 "Ext. (extended): generator (4.15) as printed, a short stock position financed by reverse repo in the crisis;",
                                 "Adm. (admissible): the control projected onto Z >= 0 in the crisis (the convention of the thesis).",
                                 "The seller is the same under both. Inadmissible nodes: of %s crisis nodes of the (50,100) tree; --: no crisis nodes (normal regime).",
                                 "Under the admissible convention the crisis effect on XVA^- scales with the lost reverse-repo carry,",
                                 "from %.2f%% at r_r^+ = %.2f to %.2f%% at r_r^+ = %.2f; inverted band points over the %d points of the",
                                 "comparison grid: %d extended, %d admissible, all of the latter in the normal regime."),
                           formatC(max(ad$crisis_nodes), format = "d", big.mark = ","),
                           ca$effect_pct_Vhat0[1], ca$rr_p[1], ca$effect_pct_Vhat0[nrow(ca)], ca$rr_p[nrow(ca)],
                           nrow(r$values) / 2, inv[["extended"]], inv[["banned"]]))
  record_output("tab16", "tab16_closeout.tex", "buyerconv", "App. C", "the buyer under the two crisis conventions")
}

# ------------------------------------------------------------------------------

##  Extras: items no longer in the thesis, kept callable by name

#-------------------------------------------------------------------------------
x_published_ratios <- function() {
  pub <- suite_data("ps")$ratios
  df <- data.frame(f = pub$figure, rf = fmt(pub$rf_m, 2), a = fmt(pub$alpha, 1),
                   p = paste0("$\\approx", formatC(pub$published, format = "f", digits = 4), "$"),
                   q = fmt(pub$quantisation, 6), pde = fmt(pub$pde, 6),
                   ratio = fmt(pub$ratio_published_to_reproduced, 2), stringsAsFactors = FALSE)
  tex_table("x_published_ratios.tex", df, "lrrrrrr",
            c("Figure", "$r_f^-$", "$\\alpha$", "Published", "Quantisation", "PDE", "Ratio"),
            note = "Published values read off the article's figures by eye (V4 Table 6).")
  record_output("x_published_ratios", "x_published_ratios.tex", "ps", "not in the thesis", "published against reproduced")
}
x_lambda_convention <- function() {
  f5 <- suite_data("ps")$figure5; x <- f5[f5$rf_m == 0.10, ]
  fig_open("x_lambda_convention.pdf")
  for (cv in c("mean", "rate")) {
    panel(c(0, 1), range(x$quantisation), expression(alpha), expression(XVA[0]^"+"),
          main = if (cv == "mean") expression("(a) mean holding times"~lambda[U]~"(years)") else
                   expression("(b) read as rates"~a[U]~"(per year)"))
    y <- x[x$regime == "frozen0", ]; series(y$alpha, y$quantisation, 1)
    for (j in 1:3) { l <- c(0.3, 1, 1.7)[j]; y <- x[x$regime == "dynamic" & x$lambda_U == l & x$convention == cv, ]
      series(y$alpha, y$quantisation, j + 1) }
    leg("topleft", c(expression(beta[t]==0), lab_expr(paste0(if (cv == "mean") "lambda[U]=='" else "a[U]=='", lam_lab(c(0.3, 1, 1.7)), "'"))), 1:4)
  }
  fig_close()
  record_output("x_lambda_convention", "x_lambda_convention.pdf", "ps", "not in the thesis", "Figure 5 under the two readings")
}
x_buyer_convention <- function() {
  v <- suite_data("buyerconv")$values; v <- v[v$rf_m == 0.10, ]
  fig_open("x_buyer_convention.pdf")
  for (cv in c("extended", "banned")) {
    panel(c(0, 1), range(v$XVA_minus_pct_Vhat0), expression(alpha), expression(XVA[0]^"\u2212"~"(% of"~hat(V)[0]*")"),
          main = if (cv == "extended") "(a) extended equation" else "(b) short sale banned in the crisis")
    abline(h = 0, col = MUTED, lwd = 0.8)
    for (rg in names(REGIME_SLOT)) { x <- v[v$convention == cv & v$regime == rg, ]
      series(x$alpha, x$XVA_minus_pct_Vhat0, REGIME_SLOT[[rg]], type = "b") }
    if (cv == "banned") leg(0, y = 2.9, lab_expr(regime_name(names(REGIME_SLOT))), REGIME_SLOT, type = "b")
  }
  fig_close()
  record_output("x_buyer_convention", "x_buyer_convention.pdf", "buyerconv", "not in the thesis", "buyer XVA under the two crisis conventions")
}
x_funding_profile <- function() {
  pf <- suite_data("bands")$funding_profile; pf <- pf[pf$side == "+" & pf$alpha == 0.5, ]
  fig_open("x_funding_profile.pdf")
  for (col in c("E_borrowing", "E_lending")) {
    panel(c(0, 0.25), c(0, max(pf[[col]]) * 1.05), expression(t~"(years)"),
          if (col == "E_borrowing") expression(E*"["*A[t]^"\u2212"*"]") else expression(E*"["*A[t]^"+"*"]"),
          main = if (col == "E_borrowing") "(a) expected borrowing" else "(b) expected lending")
    for (rg in names(REGIME_SLOT)) { x <- pf[pf$regime == rg, ]; series(x$t, x[[col]], REGIME_SLOT[[rg]]) }
    if (col == "E_borrowing") leg(0, y = 0.5, lab_expr(regime_name(names(REGIME_SLOT))), REGIME_SLOT)
  }
  fig_close()
  record_output("x_funding_profile", "x_funding_profile.pdf", "bands", "not in the thesis", "expected funding profile")
}

# ------------------------------------------------------------------------------

##  The registry and the runner

# ------------------------------------------------------------------------------
PLOT_ITEMS <- c(
  fig01 = "5.2    Convergence of value and control",
  fig02 = "5.3    XVA against alpha: normal | crisis, lines by r_f^-",
  fig03 = "5.3    XVA against alpha: normal vs dynamic, lines by lambda_U",
  fig04 = "5.3    (lambda_U, lambda_V) heat map",
  fig05 = "App. C Regime effect against the spread",
  fig06 = "5.4    Seller and buyer against alpha by r_f^- (band)",
  fig07 = "App. C Repo-rate and credit-intensity panels",
  fig08 = "5.5    Stock units against S and time, normal vs crisis",
  fig09 = "5.5    Positions over alpha and r_f^-",
  fig10 = "5.5    Replication error and its decay",
  fig11 = "5.6    Stacked Shapley contributions",
  fig12 = "5.7    Frozen-grid against full-grid vega",
  fig13 = "5.7    AAD cost scaling and delta cross-check",
  fig14 = "5.8    Accuracy against run time, by method",
  fig15 = "App. B Bias against lambda_U with seed dispersion",
  fig16 = "App. B Penalty sweep and final training loss",
  tab01 = "5.1    Parameters and solver settings",
  tab02 = "5.2    Errors in value and control by regime and grid",
  tab03 = "5.3    Seller's XVA by regime, alpha, r_f^-",
  tab04 = "5.3    Model risk: frozen against dynamic",
  tab05 = "5.4    Seller, buyer and band width",
  tab06 = "5.5    Six accounts at inception, with balance check",
  tab07 = "5.6    Shapley values and spread across orderings",
  tab08 = "5.7    AAD sensitivities and agreement with bumps",
  tab09 = "5.8    Method performance summary",
  tab10 = "App. C Dollar funding positions, seller and buyer",
  tab11 = "App. C Replication error against rebalancing step",
  tab12 = "App. B Deep solvers by seed and lambda_U",
  tab13 = "App. B Penalty-weight sweep",
  tab14 = "App. C AAD against PDE bumps by regime; frozen vs full",
  tab15 = "App. D Run-time record",
  tab16 = "App. C Buyer under the two crisis conventions")
PLOT_EXTRAS <- c(x_model_risk_ranges = "model-risk ranges over the heat map (formerly in Table 4)",
                 x_published_ratios = "published against reproduced (Pang and Sturm)",
                 x_lambda_convention = "Figure 5 under the two readings of lambda",
                 x_buyer_convention  = "buyer XVA under the two crisis conventions",
                 x_funding_profile   = "expected funding profile")

## Locale for the plotmath symbols and a check of the working directory; called
## by run_plots(), and once by hand before drawing items one at a time.
plot_setup <- function() {
  if (.Platform$OS.type != "windows") try(Sys.setlocale("LC_CTYPE", "C.UTF-8"), silent = TRUE)
  if (!dir.exists(OUTPUT_DATA)) stop("run from the XVA-Code folder: ", OUTPUT_DATA, " not found")
  invisible(TRUE)
}

## run_plots(): every item of the one-pager, or those named ("extras" for the
## items no longer in the thesis).  The figure map is merged, so drawing one
## item keeps the others' rows.
run_plots <- function(items = character(0), ...) {
  plot_setup()
  if (length(items) == 0) items <- names(PLOT_ITEMS)
  if ("extras" %in% items) items <- c(setdiff(items, "extras"), names(PLOT_EXTRAS))
  bad <- setdiff(items, c(names(PLOT_ITEMS), names(PLOT_EXTRAS)))
  if (length(bad)) stop("unknown item(s): ", paste(bad, collapse = ", "),
                        "; known: ", paste(names(PLOT_ITEMS), collapse = " "), " extras")
  t0 <- proc.time()[["elapsed"]]; failed <- character(0)
  for (it in items) {
    res <- try(get(it)(), silent = TRUE)
    if (inherits(res, "try-error")) { while (dev.cur() > 1) dev.off(); failed <- c(failed, it)
      message("plots: ", it, ": ", conditionMessage(attr(res, "condition"))) }
  }
  new <- if (length(ls(.fig_map))) do.call(rbind, mget(ls(.fig_map), envir = .fig_map)) else NULL
  f <- file.path(PLOT_SETTINGS$fig_dir, "figure_map.csv")
  old <- if (file.exists(f)) tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL) else NULL
  if (!is.null(old) && !"item" %in% names(old)) old <- NULL           # earlier format
  mp <- rbind(if (!is.null(old) && !is.null(new)) old[!old$item %in% new$item, names(new)] else old, new)
  if (!is.null(mp)) { mp <- mp[order(mp$item), ]; write.csv(mp, f, row.names = FALSE) }
  structure(list(items = items, failed = failed, map = mp, seconds = proc.time()[["elapsed"]] - t0),
            class = "plots_result")
}

report_plots <- function(res) {
  ok <- setdiff(res$items, res$failed)
  cat(sprintf("\n%d item(s) drawn in %.1f s; map in %s\n", length(ok), res$seconds,
              file.path(PLOT_SETTINGS$fig_dir, "figure_map.csv")))
  if (length(ok)) print(res$map[res$map$item %in% ok, c("item", "file", "where")], row.names = FALSE)
  if (length(res$failed)) cat("FAILED:", paste(res$failed, collapse = ", "), "\n")
}
