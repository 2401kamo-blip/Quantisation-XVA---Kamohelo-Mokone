# Replication of Pang and Sturm, *XVA Valuation under Market Illiquidity*

Base-R implementation of the numerical study of arXiv:2011.03543v1, solved with
the deep BSDE method of E, Han and Jentzen (2017), the jump-augmented deep solver
of Andersson, Gnoatto, Patacca and Picarelli (2025) and a quantisation solver
(Callegaro, Gnoatto and Grasselli 2023), together with two independent
benchmarks that the article does not provide. Sensitivities of the XVA are
computed by adjoint algorithmic differentiation of the quantisation solver.

No packages beyond those shipped with R are required.

## Running

```
Rscript main.R          # assumption report, validation, Figure 4, Figure 5
Rscript main.R quick    # validation only, reduced settings
Rscript main.R regress  # regression gate: V4 tables + full-precision snapshot
Rscript main.R verify   # suite (a): verification against closed form and PDE (add 'force' to recompute)
Rscript main.R replicate # suite (b): replicating portfolio, admissibility, self-financing simulator
Rscript main.R shapley  # suite (c): Shapley decomposition (24 and 120 orderings)
Rscript main.R bands    # suite (d): price band, funding positions, repo and credit sensitivity
Rscript main.R regime   # suite (e): (lambda_U, lambda_V) heat maps, spread effect, model risk
Rscript main.R deep     # suite (f): deep solvers against quantisation, four seeds, K sweep (about an hour)
Rscript main.R aad      # suite (g): AAD sensitivities against bumps by level, frozen vega, cost, delta check
Rscript main.R valuation # suite (h): XVA grid over alpha, r_f^- and lambda_U with all solvers ("pangsturm" also works)
Rscript main.R closeout  # close-out comparison: the buyer under both crisis conventions ("buyerconv" also works; docs/buyer_closeout.md)
Rscript main.R plots    # figures (figures/, 6.30 in) and booktabs tables (tables/) from output/data
```

Output is written to `./output` as `figure4.pdf`, `figure5.pdf` and the
corresponding `.rds` result objects.

Indicative cost on one core: the validation suite takes a few minutes; each
deep BSDE training run at the default settings (`N = 20`, `M = 128`,
1500 iterations) takes roughly half a minute, and the two figure sweeps together
call it 27 times. The two benchmark solvers are much cheaper, and the
quantisation solver costs about one second per value at `n = 50`, `N = 100`, so
it runs on the fine alpha grid alongside the PDE. Reduce
`SETTINGS$dbsde_iters` in `R/09_experiments.R` for a fast pass.

### Figures from the committed results
### important for bypassing computational time needed for entire pipeline
The suite results (`output/data/*.rds`, about 5 MB) are committed, so the figures and
tables can be redrawn without rerunning any experiment:

```
cd XVA-Code
Rscript main.R plots              # all 16 figures and 16 tables of Chapter 5, in a few seconds
Rscript main.R plots fig02 tab03  # only the items named
Rscript main.R plots extras       # items not included in write up but are put of the methodology output
```

The items (figures and tables) follow this nubmering convention: `figures/fig01_convergence.pdf`
to `fig16_penalty.pdf` and `tables/tab01_parameters.tex` to `tab16_closeout.tex`. In
RStudio, run the lines of `main.R` that load the modules, then any line of the block
"Figures and tables, one at a time" (for example `fig02()`): each redraws one item.
`figures/figure_map.csv` lists each item with its section, data file and purpose.


## File map

| File | Content |
|---|---|
| `R/00_parameters.R` | benchmark coefficients of Section 7.2; checks on Assumption 12, Proposition 14 and Assumption 17 |
| `R/01_status_process.R` | the alternating renewal process, its two-state Markov representation and its compensator |
| `R/02_black_scholes.R` | the reference pair `(Vhat, Zhat)` of Equation (12) |
| `R/03_paths.R` | exact forward simulation of `(S, beta)` under Q |
| `R/04_driver.R` | the driver in the paper's three layers, with analytic partial derivatives |
| `R/05_network.R` | batch-normalised network, hand-coded backpropagation, Adam |
| `R/06_deep_bsde.R` | the deep BSDE solver, with and without the compensated-jump control |
| `R/07_pde_benchmark.R` | coupled regime-switching semilinear PDE benchmark |
| `R/08_linear_benchmark.R` | closed-form representation via a Girsanov shift |
| `R/09_experiments.R` | the sweeps behind Figures 4 and 5 |
| `R/10_validation.R` | component-by-component validation suite |
| `R/11_deep_bsde_jumps.R` | deep solver for BSDEs with jumps (Andersson, Gnoatto, Patacca & Picarelli 2025), plus `compare_deep_solvers()` and `compare_jump_schemes()` |
| `R/12_quantisation_grid.R` | recursive marginal quantisation of log S: Newton-Raphson on a tridiagonal Hessian, Lloyd fallback, transition probabilities |
| `R/13_quantisation_solver.R` | the quantisation solver: scheme (3.1) of Callegaro, Gnoatto and Grasselli on the product state (S, beta) |
| `R/14_quantisation_adjoint.R` | hand-coded adjoint (AAD) of the quantisation solver: delta, vega, r_f^-, h_I, h_C for both sides; grid differentiated by the implicit function theorem |
| `R/16_regression.R` | regression gate: the deterministic tables of the V4 record (half a unit in the last printed digit) and a full-precision snapshot (`regression/snapshot.rds`, relative 1e-9) |
| `R/17_verification.R` | suite (a): quantisation against the closed form and the PDE, frozen and dynamic, joint refinement (25,50)-(100,200), errors in value and Z0, Richardson extrapolation; writes `output/data/verification_*` |
| `R/18_replication.R` | suite (b): the six accounts in units and currency from Z, Z^I, Z^C with a balance check; inadmissible crisis nodes; self-financing simulator (seller) and the decay of the replication error with the rebalancing step |
| `R/19_shapley.R` | suite (c): exact Shapley decomposition of XVA_0 over funding spread, DVA, CVA, regime (24 orderings) and with collateral (120), with the definitions of "off" |
| `R/20_bands.R` | suite (d): seller and buyer XVA and the price band over alpha, r_f^- and regime (inversions flagged); funding positions at t = 0 and expected funding profiles; repo-rate and credit-intensity sensitivity |
| `R/21_regime.R` | suite (e): 15 x 15 (lambda_U, lambda_V) heat maps for both sides and initial regimes; regime effect against r_f^- - r_r^-; frozen-versus-dynamic model risk; PDE spot checks |
| `R/22_deep.R` | suite (f): plain, jump-tied, jump-untied deep solvers (and the exact-compensator variant) against quantisation and PDE over lambda_U, four seeds, penalty sweep; per-run cache |
| `R/23_aad.R` | suite (g): adjoint sensitivities against same-scheme bumps at every level and against PDE and extrapolated closed-form bumps, with Richardson; frozen-grid versus full vega across alpha; cost by level; delta against Zbreve_0/(sigma S0) |
| `R/24_pang_sturm.R` | suite (h): Pang and Sturm Figures 4 and 5 with closed form, PDE, quantisation and deep solvers (four seeds); both lambda conventions; published-versus-reproduced ratios |
| `R/25_buyer_convention.R` | the buyer's crisis convention: extended equation (default) against the short-sale ban on the control; values, band, carry, PDE check, admissibility, sensitivities |
| `R/26_plots.R` | reads output/data only; one function per item of the Chapter 5 one-pager (`fig01()` to `fig16()`, `tab01()` to `tab16()`), vector PDF figures at the thesis text width (Okabe-Ito palette, line-type and marker redundancy, plotmath symbols), booktabs .tex fragments, and figures/figure_map.csv |
| `R/15_sensitivity_validation.R` | the adjoint against bump-and-revalue on the same scheme (V8-adjoint sensitivities against bump-and-revalue) and on the PDE and closed form (V9-adjoint sensitivites against PDE and closed form); frozen-grid vega; cost |

