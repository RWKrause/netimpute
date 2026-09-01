# netimpute (development version)

## New features

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
