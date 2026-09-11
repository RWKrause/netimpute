# netimpute 1.1.0

Tie imputation was over-imputing: `netmice()` 1.0.0 filled missing cells at
about **1.5x the true tie rate**, and on the sparsest networks at 3-5x. The
true rate among the missing cells equalled the observed density, so there was
nothing for the imputation to correct for - the surplus was model bias. The
changes below address it. **Imputed tie values change for every network
target.**

## Breaking changes

* **The `PCA` budget is now the *total* predictor count, protected columns
  included.** `.clean_predictor_matrix()` previously applied the cap only to
  the columns it was about to collapse and then appended the protected block
  (a network target's own endogenous terms plus one `_tie`/`_recip` pair per
  other network) on top, so a tie model carried `budget + length(keep_raw)`
  predictors. With eight networks that is 16 extra, and the
  events-per-variable rule `PCA$ratio` documents (Peduzzi et al. 1996;
  Harrell 2015) was missed by that margin - worst on sparse networks, where
  the budget is smallest and the overshoot proportionally largest. On an
  eight-network dataset at 30% missing, the realised events per variable rose
  from 4.8-8.9 to 9.1-10.2 against an intended 10.

  Protected columns still keep their own named coefficients; they now spend
  the budget rather than riding on top of it, and the collapsible remainder
  is reduced accordingly. At least one component always survives, so a
  protected block that fills the budget cannot silently discard every
  attribute and cross-network predictor. When it exceeds the budget outright
  `netmice()` raises a `netimpute_budget_overflow` warning naming the target.

* **`models` terms are counted against the budget too.** They were appended
  with no cap at all, a second route past the same rule. They are still never
  collapsed into components - that would defeat the congeniality guarantee
  `models` exists to provide - but they are now built first and spend the
  budget, leaving the auto-generated block to collapse into what remains.

* **`net_sweeps` defaults to 5**, so a network's missing cells are drawn five
  times per visit rather than once. Imputed ties differ from 1.0.0 as a
  result. Set `net_sweeps = 1` for the old behaviour.

* **`reciprocity` is no longer emitted for an undirected target.** There
  `y_ji` is `y_ij`, so the column was a perfect predictor and the working
  model separated on it - visible only as a suppressed fit warning and a
  tolerated `chol()` failure. Undirected results change accordingly, and
  `dyad_regression()` no longer returns that column for symmetric input.

* **A failed Cholesky no longer disables the coefficient draw** - see Bug
  fixes. This changes imputed values wherever it used to trigger.

## New features

* **`net_sweeps`** (default `5`) makes K Gibbs passes over a network's missing
  cells per visit instead of one, each in a fresh random order, conditioning
  on every tie drawn so far. The tie updater is the full-conditional Gibbs
  sampler of an ERGM estimated by pseudo-likelihood; `net_sweeps` controls how
  far the network mixes toward the distribution its coefficients imply, while
  `maxit` controls the outer MICE chain and is where those coefficients are
  re-estimated. Coefficients are drawn once per visit and reused across the K
  passes, so raising it costs sweep time but no extra model fits (~0.14 s per
  pass per network at n = 100).

  Drawing each cell once per visit left the per-iteration tie density still
  drifting at `maxit = 20`. At K = 5 the imputed out-degree dispersion moves
  closer to truth and the correlation between true and imputed out-degree rose
  from 0.24 to 0.34 on the development data. **Imputed values change.**

* **`net_endo_terms` and `net_gw_decay`** expose the tie model's endogenous
  statistics: `reciprocity`, the bounded `twopath` indicator, and the
  geometrically weighted change statistics `gwesp`, `gwodegree`, `gwidegree`.
  These are true ERGM change statistics - the movement in the statistic if the
  cell went 0 to 1, evaluated on the network with that cell removed - which
  makes the degree terms leave-one-out by construction, so the tie being
  imputed is never a summand of its own predictor.

  The default is `c("reciprocity", "twopath")`, i.e. unchanged from 1.0.0.
  **The `gw*` terms are implemented, tested and documented but not on by
  default**: on the development data they reproduced triadic closure better
  than `twopath` did, yet bought it by inventing ties, worst where ties are
  rarest (at density 0.008 the imputed tie rate went from 2.1x truth to 4.3x,
  and mean Brier skill across eight networks from 0.12 to -0.07). `twopath` is
  bounded by 1; `gwesp`'s change statistic is bounded by
  `exp(decay) + 2(n - 2)`, so a positive coefficient on it is a strong tie
  generator. They may well pay off on denser networks or at a lower decay -
  that is what the argument is for.

  Which terms a target actually receives depends on its type: an undirected
  target has no `reciprocity` (there `y_ji` *is* `y_ij`, which makes the
  column a perfect predictor and separates the fit), a weighted target keeps
  `twopath` rather than the binary-ERGM `gw*` statistics, and an undirected
  binary target gets a single summed `gwdegree`. The set is filtered rather
  than rejected, because `networks` is routinely a mixed list; a `models`
  formula naming a term its target does not have is a clear error.

* **The linear predictor is clamped** to +/-30 before `plogis()` in the
  sequential updater, so no endogenous statistic can pin every probability at
  exactly 0 or 1 and turn the sweep deterministic.

* **`net_ridge`** adds an optional ridge penalty to the tie working model,
  **off by default** (`lambda = 0`). It is off because, on the data it was
  developed against, the total-budget fix alone already brought the imputed
  tie rate to within 2% of truth, and adding `lambda = 0.01` on top made the
  level bias slightly worse while cutting node-level discrimination sharply
  (mean Brier skill 0.17 to 0.04, mean correlation between true and imputed
  out-degree 0.26 to 0.02). A ridge shrinks the worst-identified directions
  hardest, and on a dyad design that is where the between-node signal lives.
  The argument is provided, tested and documented so the penalty can be
  measured properly rather than assumed. The fit is
  routinely done at a handful of events per predictor, and an unpenalised
  logistic fit there produces a linear predictor whose spread is far too
  large: fitted probabilities pile up near 0 and 1. Because `plogis()` is
  convex below 0.5, an over-dispersed linear predictor inflates the *mean*
  imputed tie probability wherever ties are rare - which is exactly the
  observed over-imputation.

  Predictors are standardised on the observed dyads before the penalty is
  applied and mapped back afterwards, so one `lambda` is meaningful across a
  design mixing 0/1 dyad indicators, principal components and attribute
  differences. Every non-intercept coefficient is penalised, including the
  endogenous and cross-network terms protected from the PCA collapse - those
  are the columns most likely to be numerous relative to the events
  available. The intercept is never penalised, so the baseline density stays
  free to match the observed rate. `scale = "epv"` makes `lambda` track
  predictors-per-event so a thinly identified target is penalised harder.

  The Bayesian coefficient draw uses the penalised covariance and is
  therefore also shrunk. That is intended, for the same convexity reason, but
  it does reduce between-imputation variance: a large `lambda` will
  eventually make the imputations improperly narrow.

* **`PCA = list(n = 0)`** is now accepted and means "no components": nothing
  is collapsed and only the protected predictors are kept (a network target's
  endogenous and cross-network dyad terms; an attribute model's isolate flags
  and `models` terms, or an intercept-only model if it has none). It works for
  `PCA_attributes` and `PCA_networks` too. Previously `n = 0` was rejected, and
  every budget floored at one component; a budget derived from `ratio` still
  does, and a fractional `n` below 1 still rounds up to 1 rather than silently
  becoming a zero request.

* **`PCA_networks`** budgets the tie models independently of `PCA`, mirroring
  `PCA_attributes`. Until now only the attribute side could be overridden,
  which had it backwards: an attribute model is budgeted against observed
  rows and a tie model against observed *events*, and it is the tie models
  that run closest to their budget, since the protected dyad block grows by
  two columns for every additional network and cannot be collapsed. `"none"`
  imposes no budget on tie models.

## Bug fixes

* **A failed Cholesky factorisation of the coefficient covariance no longer
  silently disables the Bayesian draw.** `.impute_ties_gibbs()` wrapped
  `chol(vcov(fit))` in `tryCatch()` and, on failure, left *every* coefficient
  at its point estimate - so between-imputation variance collapsed to zero
  with no warning, exactly when the fit was least trustworthy. The diagonal
  is now nudged and the draw proceeds.

## Also new in 1.1.0

These were listed as unreleased development changes; the version was never bumped for them, so they ship in 1.1.0 too.

### Breaking changes

* **`netmice()`'s `net_list` argument is now called `networks`**, and the same
  rename applies to `dyad_regression()`, `net_predictors()` and
  `netquickpred()`. The fitted object's `$net_list` element and
  `complete_netmice()`'s `$net_list` are now `$networks`. There is no
  deprecated alias: existing scripts that pass `net_list =` by name, or read
  `fit$net_list`, must be updated.

* **`n_components` is replaced by `PCA`**, a list with `n` (a fixed maximum
  number of predictors/components) and/or `ratio` (observations per
  predictor); the smaller budget wins. This is now the *single* dimensionality
  safeguard - it also replaces the internal `max(5, n/3)` cap that previously
  limited every imputation model.

  The default is `PCA = list(n = NULL, ratio = 10)`. For a **binary** target -
  any tie model, or a binary attribute - `ratio` counts *events*, the rarer of
  the two outcomes among the observed values, rather than rows. A logistic
  model is limited by its rarer class: a 30-node network at density 0.15 has
  ~870 dyad rows but only ~130 ties, and budgeting against the rows would
  permit enough predictors to separate the data perfectly. `ratio = 10` is the
  conventional events-per-variable floor (Peduzzi et al. 1996; Harrell 2015).

  **Imputed values change** as a result: models that previously carried up to
  `n/3` predictors now typically carry fewer, and mixed-model fits that were
  borderline may now fall back to standard PMM.

### New features

* **`netmice()` gains `PCA_attributes`**, an optional separate predictor
  budget for the attribute imputation models. `NULL` (the default) inherits
  `PCA` and reproduces earlier behaviour exactly; `"none"` imposes no budget,
  so attribute predictors keep their own coefficients instead of being
  collapsed to principal components; a `list(n=, ratio=)` sets an
  attribute-only budget. Tie models always use `PCA`.

  The two sides are not comparable: an attribute model is budgeted against
  observed *rows*, a tie model against observed *events*, which at typical
  network sizes differ by an order of magnitude, so one `ratio` need not suit
  both. `"none"` is intended for a short, deliberately chosen predictor set
  whose individual coefficients are the point; it removes the safeguard, so a
  wide predictor set against few observed rows can leave the univariate
  models near-singular.

* **Networks may be supplied as an edgelist.** `networks` accepts a data.frame
  instead of a list of matrices, described by the new `edgelist_options`
  argument (`edgelist_names`, `edgelist_format`, `edgelist_split`, `nodelist`,
  `missing`, `directed`). One edgelist can be split into many networks -
  `edgelist_format = "long"` splits on the interaction of the `edgelist_split`
  columns, `"wide"` makes one network per named column. `structural` accepts
  an edgelist too (its `missing` entry is ignored, since structural
  constraints are always known). See `vignette("netimpute")`.

* **New `id` argument** naming the column of `data` that holds node
  identifiers. With matrices it matches and reorders their row/column names;
  with an edgelist it is required, since node names cannot otherwise be
  resolved to rows. When it is `NULL`, positional alignment is assumed and
  said so out loud. The column is dropped from `data` after alignment, so it
  is never treated as an attribute.

* **R-hat convergence diagnostic.** `netmice()` now computes rank-normalized
  split-R-hat for every tracked quantity (attribute means and variances, the
  network diagnostics, and the imputed-tie traces), after discarding the first
  half of the iterations. It is stored on the fitted object as `$rhat`,
  reported by `print()`, and any value above 1.05 raises a warning naming the
  offending quantities and suggesting a larger `maxit`. Implemented inline -
  no new dependency - and verified to match `rstan::Rhat()` to ten decimal
  places. Constant traces (e.g. an isolate count that never moves) give `NA`
  rather than a spurious number.

* **The `isolate` flag is now always retained** whenever alter-based features
  are in play. Those features are `NA` for a node with no alters and are
  mean-filled downstream, which silently hands an isolate the *average
  neighbourhood of the connected nodes*; the flag is the missing-data
  indicator that lets a model offset that fill, so it now bypasses
  `netquickpred()`'s screening and is never absorbed into a principal
  component. A constant flag (no isolates, or nothing but isolates) is still
  dropped, since it offsets nothing.

  With several networks the flags are kept from multiplying. A network's flag
  is force-kept only when at least one alter-based feature from *that* network
  entered the model - offsetting that network's mean-fill is the flag's whole
  purpose - and near-duplicate flags prune each other, since the same people
  are typically isolated in every network. (Five identical flags are rank 1
  and would spend four parameters on nothing; they remain exempt from being
  pruned by ordinary predictors.) Exactly-duplicated flags are dropped even
  when `collin_method = "none"`.

* **`netmice()` now reports when a `predictor_selection` is still larger than
  the `PCA` budget**, once per call, naming the affected targets. The budget
  applies on top of the selection and at ordinary network sizes usually still
  binds, so a user who selected predictors expecting named coefficients would
  otherwise silently get principal components instead.

* **`netmice()` now recommends supplying `models`** when it is `NULL`, once
  per session, citing the congeniality requirement (Meng 1994) and, for the
  network case, Krause et al. (2020). `printFlag = FALSE` silences it.

* `netmice()` now also tracks the mean and variance of the **imputed ties
  only**, per network per iteration, in `netImpMean`/`netImpVar` — the tie-side
  counterpart of the `chainMean`/`chainVar` it already recorded for attributes.
  `plot()` draws them on a third page, laid out like the attributes page.

  The existing `netChain` diagnostics are unchanged, but they are computed on
  the whole completed network and so are diluted by the missingness fraction:
  with 7% of dyads missing, a 0.20 swing in the density of the imputed cells
  appears as a movement of only ~0.014 in overall density. Flat `netChain`
  traces were therefore never good evidence of convergence at low missingness.
  Prefer the imputed-tie traces when judging whether to raise `maxit`.

### Documentation

* New vignette, `vignette("netimpute")`: a worked introduction covering a
  first imputation, convergence checks, pooling results with `mice::pool()`,
  the predictor naming conventions, steering models with `models`, leaner
  predictor sets via `netquickpred()`, `structural`/`net_dependence`
  constraints, and directed vs. undirected networks.

* Every exported topic now has a `\seealso` section cross-linking the
  pipeline, so the help pages are navigable from any entry point.

### Bug fixes

* `plot()` no longer leaves the caller's `par()` settings modified. The
  networks page saved `par` *after* the attributes page had already applied its
  compact layout, and `on.exit()` handlers run in registration order, so the
  compact settings were restored last — every later plot in the session came out
  with `mfrow = c(2, 2)` and shrunken margins.

* `plot()` no longer plots chains along the iteration axis when `maxit = 1`.
  `chainMean[v, , ]` drops to a plain vector when `maxit` or `m` is 1, and
  `matplot()` read the resulting length-`m` vector as one chain of `m`
  iterations. Traces are now always coerced to an explicit `maxit` x `m` matrix.

* `plot()` no longer errors on a fit whose attribute has exactly **one** missing
  value, or whose network has exactly one imputed dyad. `var()` of a single
  value is `NA`, so the whole variance trace was `NA` and `matplot()` failed with
  "need finite 'ylim' values". Such panels are now drawn empty and labelled.

* Arguments passed through `...` to `plot()` now override the method's own
  `matplot()` defaults, as documented, instead of failing with "formal argument
  matched by multiple actual arguments" — `plot(fit, col = "red")` works.

* `net_diagnostics()`'s `avg_inv_geodesic` now follows tie direction on directed
  networks. It used `igraph::distances(mode = "all")` for every network, which
  computes distances on the underlying *undirected* graph, so the statistic
  silently reported an undirected quantity while `reciprocity` beside it was
  direction-aware — and the documentation promised ordered pairs. On a 3-cycle
  plus two isolates this is 0.225 (correct) rather than 0.3. Imputation itself is
  unaffected; the statistic is diagnostic-only, and values in `netChain` from
  earlier versions are not comparable for directed networks.

* Fixed silent data corruption when a network had **exactly one** missing
  tie, or an attribute exactly one observed value. `sample(x)` treats a
  length-1 numeric `x` as `1:x`, so the tie-wise sweep iterated over `1:k`
  (with `k` the missing cell's index) instead of that single cell, drawing
  into — and overwriting — observed ties. The same applied to
  `net_init = "sample"` and to the initial fill of an attribute with one
  observed value. All shuffles and resamples now go through length-safe
  helpers.

* Directedness is now derived **once**, from the data passed to `netmice()`,
  before any cell is filled, and is carried to everything that needs it: the
  tie updater, the initialisation, the igraph views the attribute predictors
  are built from, the per-iteration diagnostics, and `netquickpred()`'s
  screening. Previously each of these re-inferred it with `isSymmetric()`
  from the current, already-filled matrix. A directed network whose
  asymmetry lies entirely in its missing cells becomes symmetric once those
  are zeroed, and was then rebuilt as undirected — which returns
  `reciprocity_ratio` as `NA` for every node, dropping a real predictor from
  that visit. `net_diagnostics()` gains an optional `directed` argument and
  `netquickpred()` an optional `net_directed` argument; both default to the
  previous inference for standalone callers.

* Undirected networks now stay undirected after imputation. Previously the
  mirror cells of an undirected network were imputed independently, so the
  result was generally asymmetric — which silently reclassified the network
  as directed for every derived measure, both in the remaining sweeps and in
  the returned object. A network is treated as undirected when its input
  matrix is symmetric in its values and in its `NA` pattern. Under
  `net_update = "gibbs"` each unordered pair is now visited once, averaging
  the two directions' linear predictors (on the probability scale for binary
  ties, the linear predictor scale for weighted ties) and writing the single
  draw to both cells. Under `net_update = "simultaneous"` disagreeing pairs
  are reconciled afterwards, at random for binary ties and by averaging for
  weighted ties. `net_init = "sample"` also mirrors its initial fills.
  Observed cells are never overwritten by any of these steps.

# netimpute 1.0.0

* First release.

* `netmice()` jointly imputes missing nodal attributes and missing network
  ties in a single chained-equations loop: within each iteration every
  attribute and every network with missing data is visited once, in a
  randomised interleaved order, and the network-derived predictors of an
  attribute and the attribute-derived predictors of a network are rebuilt
  at every visit.

* Missing ties are updated tie-wise by default (`net_update = "gibbs"`):
  each tie is redrawn from its full conditional given the ties imputed so
  far, with reciprocity and the shared-partner indicator refreshed after
  every draw via change statistics. `net_update = "simultaneous"` imputes
  all missing cells of a network at once.

* `net_measures()` computes any subset of 28 node-level structural measures
  plus the per-attribute homophily block; `net_measures_core()` and
  `net_measures_full()` are wrappers. Binary and non-negative weighted
  networks are supported.

* `dyad_regression()` fits a dyad-level (cell-level) tie model, optionally
  with social relations model random intercepts via `lme4`.

* `netquickpred()` selects predictors per imputation target, generalising
  `mice::quickpred()` to network measures and dyad-level terms.

* `net_predictors()` builds node-level predictor sets across several
  networks, optionally reduced to principal components.

* Additional `netmice()` features: structural zeros (`structural`), logical
  constraints between networks (`net_dependence`), user-specified
  imputation models with interactions (`models`), per-target univariate
  methods (`method`), and parallel chains (`ncores`).
