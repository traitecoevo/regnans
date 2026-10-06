# How should a pairwise invasibility plot spend its resident solves?
#
# Compares the seeding strategies of community_pip() on the Geritz et al. 1999
# seed-size model with its equilibrium found by fixed-point iteration of the
# model's one-generation demography rather than from the closed form, so that
# every resident is a genuine solve whose cost depends on where it starts. Reports demography evaluations and
# wall time, sequentially and under a multicore plan.
#
#   Rscript scripts/pip-benchmark.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

comm <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log",
                        harness = harness_gm99(alpha = 7, beta = 15),
                        demography_control = demographic_step_control(
                          list(equilibrium_solver_name = "equilibrium_iteration",
                               equilibrium_eps = 1e-8, equilibrium_nsteps = 1000)))

run <- function(seed, refine, workers) {
  if (workers > 1L) {
    old <- future::plan(future::multicore, workers = workers)
    on.exit(future::plan(old), add = TRUE)
  }
  pip <- community_pip(comm, control = pip_control(list(
    n_resident = 41, n_mutant = 201, n_coarse = 9, seed = seed, refine = refine)))
  data.frame(seed = seed, refine = refine, workers = workers,
             residents = nrow(pip$residents),
             evaluations = sum(pip$residents$n_evals),
             evals_per_resident = round(mean(pip$residents$n_evals), 1),
             contour_points = nrow(pip$contours),
             seconds = round(pip$elapsed, 2))
}

grid <- expand.grid(seed = c("cold", "neighbour", "interpolate"),
                    refine = c(0L, 2L), workers = c(1L, 4L),
                    stringsAsFactors = FALSE)
results <- do.call(rbind, Map(run, grid$seed, grid$refine, grid$workers))
rownames(results) <- NULL
print(results, row.names = FALSE)
