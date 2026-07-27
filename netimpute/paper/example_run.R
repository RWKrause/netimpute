# Generates the real, verified output embedded in the JSS paper's
# illustrative example (paper/netimpute-jss.tex). Run from the package root;
# stdout of this script is what the paper's Routput blocks show.
options(width = 80)
suppressMessages(library(netimpute))
set.seed(20260701)

n <- 40
friends <- matrix(rbinom(n * n, 1, 0.08), n, n)
diag(friends) <- 0
advice <- matrix(rpois(n * n, 0.35), n, n)
diag(advice) <- 0

attrs <- data.frame(
  age         = round(rnorm(n, 34, 7)),
  gender      = sample(c("F", "M"), n, replace = TRUE),
  department  = sample(c("sales", "eng", "hr"), n, replace = TRUE, prob = c(0.4, 0.4, 0.2)),
  performance = round(rnorm(n, 50, 10), 1)
)

cat("---- net_measures_core() on the friendship network ----\n")
core_out <- suppressMessages(net_measures_core(friends, attrs))
print(round(head(core_out[c("node_id", "outdegree", "indegree", "reciprocity_ratio",
                             "betweenness", "constraint")]), 3))

cat("\n---- net_measures(): a bespoke measure selection ----\n")
print(names(net_measures(friends, measure_set = c("indegree", "betweenness", "pagerank"))))

cat("\n---- dyad_regression(): friendship ties on attributes + advice network ----\n")
dr <- suppressMessages(
  dyad_regression(list(friends = friends, advice = advice), attrs, target = "friends"))
print(summary(dr$model)$coefficients[c("age_absdiff", "department_same", "advice_tie"), ])

cat("\n---- Introduce missingness in attributes and in the friendship network ----\n")
attrs_miss <- attrs
set.seed(1)
attrs_miss$age[sample(n, 4)] <- NA
attrs_miss$performance[sample(n, 5)] <- NA
attrs_miss$department[sample(n, 3)] <- NA

friends_miss <- friends
off <- which(row(friends_miss) != col(friends_miss))
friends_miss[sample(off, round(length(off) * 0.05))] <- NA

cat("Missing values per attribute:\n")
print(colSums(is.na(attrs_miss)))
cat("Missing ties in the friendship network:",
    sum(is.na(friends_miss[row(friends_miss) != col(friends_miss)])), "\n")

cat("\n---- netmice(): joint imputation (tie-wise Gibbs updating by default) ----\n")
# maxit = 5 here (below the default of 20) keeps this illustration's
# convergence tables compact; see ?netmice for why 20 is recommended.
fit <- suppressMessages(netmice(
  attrs_miss, list(friends = friends_miss, advice = advice),
  m = 5, maxit = 5, donors = 5, seed = 20260701, printFlag = FALSE,
  models = list("performance ~ friends_indegree + age * department",
                "friends ~ age_ego:reciprocity")
))
print(fit)

cat("\n---- One completed dataset ----\n")
completed <- complete_netmice(fit, 1)
print(head(completed$data))
cat("Any remaining NA in attributes:", anyNA(completed$data), "\n")
cat("Any remaining NA in friendship ties:",
    anyNA(completed$net_list$friends[row(friends) != col(friends)]), "\n")

cat("\n---- Convergence diagnostics: chain means for 'age' across iterations ----\n")
print(round(fit$chainMean["age", , ], 2))

cat("\n---- Network-level diagnostics across iterations (chain 1, friendship net) ----\n")
print(round(t(fit$netChain["friends", , , 1]), 3))

cat("\n---- netquickpred(): per-target predictor selection ----\n")
# mincor well above the 0.1 default: with only n = 40 nodes, chance
# correlations of network measures easily exceed 0.1, so a stricter screen
# keeps the illustration honest (and compact)
qp <- suppressMessages(netquickpred(
  attrs_miss, list(friends = friends_miss, advice = advice),
  targets = "performance", mincor = 0.4
))
print(qp)
