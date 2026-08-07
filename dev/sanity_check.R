# Manual smoke test - run this after devtools::load_all() to sanity-check
# the package. Not part of the built package (see .Rbuildignore).
#
# I could not execute R in the environment where this file was written, so
# please run this yourself and report any failures.

devtools::load_all()

set.seed(1)

## ---- single network / net_measures_core() -------------------------------
g1 <- igraph::sample_gnp(60, p = 0.06, directed = TRUE)
attrs1 <- data.frame(
  age    = round(rnorm(60, 35, 8)),
  gender = sample(c("F", "M"), 60, replace = TRUE),
  dept   = sample(c("sales", "eng", "hr"), 60, replace = TRUE)
)

core1 <- net_measures_core(g1, attrs1)
stopifnot(nrow(core1) == 60)
stopifnot(all(c("outdegree", "indegree", "reciprocity_ratio", "bonacich_power",
                "betweenness", "isolate", "constraint", "harmonic_closeness",
                "local_clustering") %in% names(core1)))
stopifnot(all(c("age_diff_out", "age_diff_in", "age_alter_mean", "age_alter_min",
                "age_alter_max") %in% names(core1)))
stopifnot(all(c("gender_ei_index", "gender_blau_index", "gender_alter_mode") %in% names(core1)))
stopifnot(all(c("dept_ei_index", "dept_blau_index", "dept_alter_mode") %in% names(core1)))
print(head(core1))
cat("core measures: OK -", ncol(core1), "columns\n\n")

## ---- single network / net_measures_full() --------------------------------
full1 <- net_measures_full(g1, attrs1, use_sna = TRUE)
stopifnot(nrow(full1) == 60)
print(names(full1))
cat("full measures: OK -", ncol(full1), "columns (sna installed:",
    requireNamespace("sna", quietly = TRUE), ")\n\n")

## ---- undirected network edge case ----------------------------------------
g2 <- igraph::sample_gnp(30, p = 0.1, directed = FALSE)
attrs2 <- data.frame(score = rnorm(30))
core2 <- net_measures_core(g2, attrs2)
stopifnot(all(is.na(core2$reciprocity_ratio)))  # concept N/A for undirected
cat("undirected edge case: OK\n\n")

## ---- isolate edge case ----------------------------------------------------
g3 <- igraph::make_empty_graph(n = 10, directed = TRUE)
g3 <- igraph::add_edges(g3, c(1, 2, 2, 3))
attrs3 <- data.frame(x = rnorm(10), grp = sample(c("a", "b"), 10, replace = TRUE))
core3 <- net_measures_core(g3, attrs3)
stopifnot(core3$isolate[4:10] |> all())
cat("isolate edge case: OK\n\n")

## ---- net_predictors() across a list of networks, three output modes ------
nets <- list(
  school_a = igraph::sample_gnp(40, 0.08, directed = TRUE),
  school_b = igraph::sample_gnp(55, 0.06, directed = TRUE)
)
attr_list <- lapply(c(40, 55), function(n) {
  data.frame(
    age    = round(rnorm(n, 35, 8)),
    gender = sample(c("F", "M"), n, replace = TRUE)
  )
})

meas_only <- net_predictors(nets, attr_list, measure_set = "full", output = "measures")
stopifnot(nrow(meas_only) == 95)
cat("net_predictors(output='measures'): OK -", nrow(meas_only), "rows\n\n")

pca_only <- net_predictors(nets, attr_list, measure_set = "full", output = "pca",
                            n_components = 5)
stopifnot(all(c("network_id", "node_id", paste0("PC", 1:5)) %in% names(pca_only$predictors)))
print(summary(pca_only$pca_model)$importance[, 1:5])
cat("net_predictors(output='pca'): OK\n\n")

both <- net_predictors(nets, attr_list, measure_set = "full", output = "both",
                        n_components = 5, keep_vars = c("total_degree", "betweenness"),
                        residualize_kept = TRUE)
stopifnot(all(c("total_degree", "betweenness", paste0("PC", 1:5)) %in% names(both$predictors)))
print(head(both$predictors))
cat("net_predictors(output='both', residualize_kept=TRUE): OK\n\n")

cat("ALL SANITY CHECKS PASSED\n")
