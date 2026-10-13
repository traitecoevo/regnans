# What does a warm start buy an equilibrium solve on the plant SCM?
#
# The gate for #60's continuation stages: a corrector step re-solves a
# community at parameters a step away from one already solved. This measures
# that solve started cold (birth_rate_initial, the default node schedule)
# against one started from the solved community's equilibrium birth rates,
# and from its birth rates and node schedule together, at the same hmat and
# at hmat moved by 1%. Reported in SCM runs (the demography runner's calls,
# what the solve costs) and seconds. The ratio warm / cold at the moved hmat
# is the number the roadmap records: near 0.1 continuation on plant pays,
# near 0.5 parallel ray-shooting is the better investment.
#
#   Rscript scripts/warm-start-benchmark.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

by_hmat <- function(hmat) {
  community_start(bounds(lma = c(0.02, 0.6)),
                  model_support = list(p = plant_default_assembly_pars(hmat = hmat, max_patch_lifetime = 30),
                                       plant_control = plant_default_assembly_control()),
                  derivative_control = list(d_second = 1e-2, eps_second = 1e-2))
}

## the lma singular strategy at hmat = 10, the community a continuation
## would be stepping from
root <- community_solve_singularity(by_hmat(10), x0 = 0.08, tol = 0.3)
x <- root$traits
n <- as.numeric(root$birth_rate)
schedule <- root$model_support$node_schedule_times

solve <- function(hmat, birth_rate, schedule = NULL) {
  comm <- by_hmat(hmat) |> community_add(x, birth_rate = birth_rate)
  if (!is.null(schedule)) comm$model_support$node_schedule_times <- schedule
  seconds <- system.time(out <- community_demography(comm))[["elapsed"]]
  data.frame(runs = nrow(attr(out, "progress")), seconds = round(seconds, 1),
             birth_rate = signif(as.numeric(out$birth_rate), 6),
             converged = isTRUE(attr(out, "converged")))
}

cases <- list(
  `same hmat, cold`                 = list(hmat = 10, birth_rate = comm_start <- by_hmat(10)$birth_rate_initial),
  `same hmat, birth rates`          = list(hmat = 10, birth_rate = n),
  `same hmat, birth rates+schedule` = list(hmat = 10, birth_rate = n, schedule = schedule),
  `hmat +1%, cold`                  = list(hmat = 10.1, birth_rate = by_hmat(10)$birth_rate_initial),
  `hmat +1%, birth rates`           = list(hmat = 10.1, birth_rate = n),
  `hmat +1%, birth rates+schedule`  = list(hmat = 10.1, birth_rate = n, schedule = schedule)
)
rows <- lapply(cases, function(cs) do.call(solve, cs))
out <- do.call(rbind, rows)
out$ratio_runs <- round(out$runs / out$runs[grepl("cold", rownames(out))][c(1, 1, 1, 2, 2, 2)], 2)
out$ratio_seconds <- round(out$seconds / out$seconds[grepl("cold", rownames(out))][c(1, 1, 1, 2, 2, 2)], 2)
print(out)
cat(sprintf("\nlma singular strategy at hmat = 10: %s (equilibrium birth rate %s)\n", signif(x, 5), signif(n, 6)))
