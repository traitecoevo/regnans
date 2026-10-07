# What does a future plan buy a classification on the plant SCM?
#
# Classifies the lma singular strategy (one resident: a Jacobian of two
# equilibrium solves) and a two-resident lma community (a Jacobian of four,
# plus two leave-one-out solves for protected coexistence), and takes that
# community's Jacobian alone, sequentially and under multicore plans of two
# and four workers. Reports the wall time, the equilibrium solves and how far
# the Jacobian moves from the sequential one. The second resident is near
# extinction, so one leave-one-out solve starts from almost nothing and
# dominates its step: the classification is limited by that solve, the
# Jacobian alone by its slowest stencil point.
#
#   Rscript scripts/parallel-benchmark.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

model_support <- list(p = plant_default_assembly_pars(max_patch_lifetime = 30),
                      plant_control = plant_default_assembly_control())
start <- function() {
  community_start(bounds(lma = c(0.02, 0.6)), model_support = model_support,
                  derivative_control = list(d_second = 1e-2, eps_second = 1e-2))
}
one <- community_solve_singularity(start(), x0 = 0.08, tol = 0.3)
pair <- start() |>
  community_add(trait_matrix(c(0.06, 0.3), "lma"), birth_rate = c(200, 200)) |>
  community_demography()

## the pair is not a singular coalition; the classifier says so, and the
## cost of classifying it is what is measured here
classify <- function(community) {
  cl <- suppressWarnings(community_classify_singularity(community))
  list(jacobian = cl$jacobian, solves = cl$evaluations, verdict = cl$classification)
}
jacobian <- function(community) {
  J <- community_selection_gradient_jacobian(community)
  list(jacobian = J, solves = attr(J, "evaluations"), verdict = "")
}
timed <- function(f, community, workers) {
  if (workers > 1L) {
    old <- future::plan(future::multicore, workers = workers)
    on.exit(future::plan(old), add = TRUE)
  }
  seconds <- system.time(out <- f(community))[["elapsed"]]
  c(out, seconds = seconds)
}

cases <- list(list("classify one", classify, one), list("classify pair", classify, pair),
              list("jacobian of pair", jacobian, pair))
rows <- list()
for (case in cases) {
  sequential <- NULL
  for (workers in c(1L, 2L, 4L)) {
    run <- timed(case[[2]], case[[3]], workers)
    if (is.null(sequential)) sequential <- run
    rows[[length(rows) + 1L]] <- data.frame(
      case = case[[1]], workers = workers, solves = run$solves,
      seconds = round(run$seconds, 1), speedup = round(sequential$seconds / run$seconds, 2),
      jacobian_change = signif(max(abs(run$jacobian - sequential$jacobian)) /
                                 max(abs(sequential$jacobian)), 2),
      classification = run$verdict)
  }
}
print(do.call(rbind, rows), row.names = FALSE)
