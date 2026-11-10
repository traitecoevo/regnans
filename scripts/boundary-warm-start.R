# The warm-start gate at a coexistence boundary point on the plant SCM.
#
# scripts/warm-start-benchmark.R measures it at a singular strategy. This
# measures it where #60's Newton solve evaluates its residual: the FF16 lma
# boundary along hmat (test-plant-smoke-boundary.R), solved first from its
# bracketed seed. A point one step away, in the resident's lma or in hmat, is
# then solved as boundary_residual_fn() solves it, from the centre's
# equilibrium birth rates and the default node schedule, against the same
# solve started cold (birth_rate_initial). Reported in SCM runs and seconds.
#
#   Rscript scripts/boundary-warm-start.R

suppressMessages(devtools::load_all(".", quiet = TRUE))

by_hmat <- function(hmat) {
  community_start(bounds(lma = c(0.02, 0.6)),
                  model_support = list(p = plant_default_assembly_pars(hmat = hmat, max_patch_lifetime = 30),
                                       plant_control = plant_default_assembly_control()),
                  harness = harness_plant(fitness_gradient = "model"),
                  derivative_control = list(d_second = 1e-2, eps_second = 1e-2,
                                            d_parameter = 1e-2, eps_parameter = 1e-2))
}

## the boundary point, from the bracket's ends as the smoke test seeds it
seed <- by_hmat(15.6875) |> community_add(trait_matrix(0.0775, "lma"), birth_rate = 1.14)
at <- community_solve_boundary(seed, by_hmat, y = 0.23, p = c(hmat = 15.6875), tol = 1e-3)
b <- attr(at, "boundary")
hmat <- b$p[["hmat"]]
x <- as.numeric(at$traits)
n <- as.numeric(at$birth_rate)

solve <- function(hmat, lma, birth_rate) {
  comm <- by_hmat(hmat) |> community_add(trait_matrix(lma, "lma"), birth_rate = birth_rate)
  seconds <- system.time(out <- community_demography(comm))[["elapsed"]]
  data.frame(runs = nrow(attr(out, "progress")), seconds = round(seconds, 2),
             birth_rate = signif(as.numeric(out$birth_rate), 6),
             converged = isTRUE(attr(out, "converged")))
}
step <- 1e-2
cold <- by_hmat(hmat)$birth_rate_initial
cases <- list(
  `lma +1%, cold`  = list(hmat = hmat, lma = x * (1 + step), birth_rate = cold),
  `lma +1%, warm`  = list(hmat = hmat, lma = x * (1 + step), birth_rate = n),
  `hmat +1%, cold` = list(hmat = hmat * (1 + step), lma = x, birth_rate = cold),
  `hmat +1%, warm` = list(hmat = hmat * (1 + step), lma = x, birth_rate = n)
)
out <- do.call(rbind, lapply(cases, function(cs) do.call(solve, cs)))
is_cold <- grepl("cold", rownames(out))
out$ratio_runs <- round(out$runs / rep(out$runs[is_cold], each = 2), 2)
out$ratio_seconds <- round(out$seconds / rep(out$seconds[is_cold], each = 2), 2)
print(out)
cat(sprintf("\nboundary at hmat = %s: resident lma %s (birth rate %s), y* = %s; over both points warm / cold is %.2f in runs, %.2f in seconds\n",
            signif(hmat, 6), signif(x, 5), signif(n, 5), signif(b$y, 5),
            sum(out$runs[!is_cold]) / sum(out$runs[is_cold]), sum(out$seconds[!is_cold]) / sum(out$seconds[is_cold])))
