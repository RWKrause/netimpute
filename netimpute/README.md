# netimpute

Joint multiple imputation of network ties and node attributes in R.

Missing data in network studies usually affects both sides at once: the
respondent who skipped the attribute survey is often the one whose
nominations are missing too. `netimpute` imputes both together, in one
chained-equations loop, so that each side informs the other.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("RWKrause/netimpute")
```

## Usage

Supply a data frame of node attributes and a list of adjacency matrices,
with `NA` marking both missing attribute values and unknown ties:

``` r
library(netimpute)

set.seed(1)
n <- 30
friends <- matrix(rbinom(n * n, 1, 0.15), n, n)
diag(friends) <- 0
friends[sample(which(row(friends) != col(friends)), 40)] <- NA

attrs <- data.frame(age = rnorm(n, 35, 8),
                    performance = rnorm(n),
                    gender = sample(c("F", "M"), n, TRUE))
attrs$performance[sample(n, 4)] <- NA

fit <- netmice(attrs, list(friends = friends), m = 5, maxit = 20)
completed <- complete_netmice(fit, 1)
plot(fit)   # convergence diagnostics
```

Each attribute is imputed from the other attributes and from node-level
network measures (degree, centrality, brokerage, closure, and homophily
measures computed over a node's alters). Each network's ties are imputed
from dyad-level predictors built from the attributes and the other
networks. Both predictor sets are rebuilt at every visit.

Missing ties are redrawn **tie-wise** by default: one tie at a time,
conditional on the ties imputed so far, with reciprocity and the
shared-partner indicator refreshed after every single draw.

## Further options

- `models` — protect a hypothesised relationship, including interactions
  between any internally created terms (e.g. `"friends ~ age_ego:reciprocity"`).
- `predictor_selection = "quickpred"` / `netquickpred()` — lean,
  correlation-screened predictor sets per target instead of everything.
- `structural`, `net_dependence` — structural zeros and logical
  constraints between networks.
- `net_random_intercepts` — social relations model random intercepts in
  the tie model.
- `method` — any of `mice`'s numeric-response univariate methods, globally
  or per target.

See `?netmice` for the full argument reference, including the complete
lists of dyad-level term names and implemented network measures.
