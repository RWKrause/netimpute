# netimpute (development version)

## Breaking changes

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

## New features

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

## Documentation

* New vignette, `vignette("netimpute")`: a worked introduction covering a
  first imputation, convergence checks, pooling results with `mice::pool()`,
  the predictor naming conventions, steering models with `models`, leaner
  predictor sets via `netquickpred()`, `structural`/`net_dependence`
  constraints, and directed vs. undirected networks.

* Every exported topic now has a `\seealso` section cross-linking the
  pipeline, so the help pages are navigable from any entry point.

## Bug fixes

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
