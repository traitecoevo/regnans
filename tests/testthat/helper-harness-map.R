# A minimal test-only harness whose demography runner is an arbitrary map
# `n -> map(n)`.
#
# The shipped explicit harnesses (DD99, GK98, ...) hand back their equilibrium
# analytically from the demography runner, so the equilibrium solvers never
# actually have to *iterate* on them: the DD99 tests in test-demography-solvers.R
# verify dispatch, util_nlsolve integration and cleanup, but not convergence.
# This harness closes that gap. With a Ricker map the fixed point is known
# (n* = K) and the iteration's convergence rate is tunable via r, so we can set
# up the case the root-finding solvers exist for -- one where
# equilibrium_iteration has not converged in the allotted steps.
#
# It also lets a species be driven to zero, which is what exercises the
# extinct-species re-check inside equilibrium_hybrid.

harness_map <- function(map, trait_names = "x", label = "map") {
  h <- list(
    map = map,
    trait_names = trait_names,
    label = label,
    fns = list(
      parameters = function(community) {
        list(traits = community$traits, birth_rate = community$birth_rate)
      },
      make_demography_runner = function(community) {
        f <- community$harness$map
        n_calls <- 0L
        last_offspring_production <- NULL
        history <- list()
        function(birth_rates) {
          n_calls <<- n_calls + 1L
          out <- f(birth_rates)
          last_offspring_production <<- out
          history[[length(history) + 1L]] <<- list(`in` = birth_rates, out = out)
          out
        }
      },
      demography_runner_cleanup = function(community, runner, converged = TRUE) {
        e <- environment(runner)
        if (is.function(e$runner_full)) {
          runner <- e$runner_full
          e <- environment(runner)
        }
        community$birth_rate <- e$last_offspring_production
        attr(community, "converged") <- converged
        attr(community, "n_calls") <- e$n_calls
        attr(community, "progress") <- e$history
        community
      },
      viable_bounds = function(community) community,
      check_for_inviable_strategies = function(community) community$birth_rate,
      update_fitness_function = function(community) {
        # not used by the equilibrium tests; keep the pipeline happy
        community$resident_fitness <- rep(0, nrow(community$traits))
        community$fitness_function <- function(x) rep(0, NROW(x))
        community
      }
    )
  )
  class(h) <- c("harness_map", "harness")
  h
}

## Ricker: n_{t+1} = n exp(r (1 - n / K)). Fixed point n* = K, approached with
## multiplier (1 - r), so r near 2 converges slowly and oscillates.
ricker_map <- function(r, K) {
  function(n) n * exp(r * (1 - n / K))
}
