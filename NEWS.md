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
