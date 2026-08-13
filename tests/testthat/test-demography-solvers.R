# Demographic equilibrium solvers (issue #27).
#
# community_demography() dispatches on demography_control$equilibrium_solver_name
# to one of five backends. Previously only the default (equilibrium_iteration)
# was tested, and only via the plant SCM (seconds per solve); the
# equilibrium_solve_* paths -- which are the only callers of util_nlsolve -- were
# untested against the community object (see #27 "blocked/deferred").
#
# Here we drive all five through the DD99 toy harness, whose multi-resident
# equilibrium has an analytic answer (the competition linear-solve), so every
# solver can be checked against a known target in milliseconds.
#
# NOTE: the explicit (toy) harness returns its equilibrium directly from the
# demography runner, so for DD99 the solvers do not have to *iterate* to a fixed
# point the way the plant SCM does. What these tests verify is the dispatch, the
# util_nlsolve integration, and the cleanup write-back -- all against the
# analytic oracle. The genuinely-iterative convergence of the default solver is
# covered on the real SCM in test-plant-smoke.R.

dd99_pars <- list(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4)

solve_with <- function(solver, x = c(-0.5, 0.5)) {
  community_start(
    bounds(x = c(-2, 2)),
    harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4),
    demography_control = demographic_step_control(
      list(equilibrium_solver_name = solver))) |>
    community_add(trait_matrix(x, "x"), birth_rate = rep(100, length(x))) |>
    community_demography()
}

test_that("every equilibrium solver recovers the analytic DD99 equilibrium", {
  target <- dd99_equilibrium(c(-0.5, 0.5), dd99_pars)
  solvers <- c("equilibrium_iteration", "single_step", "equilibrium_hybrid",
               "equilibrium_solve_nleqslv", "equilibrium_solve_dfsane")
  for (s in solvers) {
    comm <- solve_with(s)
    expect_true(attr(comm, "converged"), info = s)
    expect_equal(as.numeric(comm$birth_rate), target, tolerance = 1e-5, info = s)
    # at the equilibrium each resident's invasion fitness is ~0
    expect_equal(comm$resident_fitness, c(0, 0), tolerance = 1e-6, info = s)
  }
})

test_that("the nleqslv and dfsane solvers agree with the default iteration", {
  ref <- as.numeric(solve_with("equilibrium_iteration")$birth_rate)
  expect_equal(as.numeric(solve_with("equilibrium_solve_nleqslv")$birth_rate),
               ref, tolerance = 1e-5)
  expect_equal(as.numeric(solve_with("equilibrium_solve_dfsane")$birth_rate),
               ref, tolerance = 1e-5)
})

test_that("community_demography rejects an unknown solver", {
  comm <- community_start(
    bounds(x = c(-2, 2)), harness = harness_dd99(),
    demography_control = demographic_step_control(
      list(equilibrium_solver_name = "not_a_solver"))) |>
    community_add(trait_matrix(0, "x"), birth_rate = 100)
  expect_error(community_demography(comm), "Unknown solver")
})

test_that("a three-resident community also solves to the analytic equilibrium", {
  x <- c(-1, 0, 1)
  target <- dd99_equilibrium(x, dd99_pars)
  comm <- solve_with("equilibrium_solve_nleqslv", x = x)
  expect_equal(as.numeric(comm$birth_rate), target, tolerance = 1e-5)
  expect_equal(comm$resident_fitness, rep(0, 3), tolerance = 1e-6)
})

# ---- genuine fixed-point solving -------------------------------------------
#
# The DD99 tests above check dispatch and write-back, but the explicit harness
# hands its equilibrium straight back, so no solver ever has to iterate. These
# tests use helper-harness-map.R, whose demography runner is a real map with a
# known fixed point, to check that the root finders actually solve -- including
# the case they exist for, where the fixed-point iteration is still far from
# convergence when its step budget runs out.

map_community <- function(map, solver, n0, nsteps = 30) {
  community_start(bounds(x = c(-2, 2)), harness = harness_map(map),
                  trait_scale = "linear",
                  demography_control = demographic_step_control(list(
                    equilibrium_solver_name = solver,
                    equilibrium_nsteps = nsteps,
                    equilibrium_eps = 1e-5))) |>
    community_add(trait_matrix(rep(0, length(n0)), "x"), birth_rate = n0) |>
    community_demography()
}

test_that("the root finders solve a fixed point the iteration has not reached", {
  # Ricker with r = 1.9 approaches n* = K with multiplier (1 - r) = -0.9, so
  # after 30 steps it is still ~0.9^30 = 4% of the way out.
  map <- ricker_map(r = 1.9, K = 10)

  iter <- map_community(map, "equilibrium_iteration", n0 = 4)
  expect_false(attr(iter, "converged"))
  expect_gt(abs(as.numeric(iter$birth_rate) - 10), 1e-3)

  for (solver in c("equilibrium_solve_nleqslv", "equilibrium_solve_dfsane")) {
    sol <- map_community(map, solver, n0 = 4)
    expect_true(attr(sol, "converged"), info = solver)
    expect_equal(as.numeric(sol$birth_rate), 10, tolerance = 1e-5, info = solver)
  }
})

test_that("the hybrid solver reaches the fixed point the iteration missed", {
  sol <- map_community(ricker_map(r = 1.9, K = 10), "equilibrium_hybrid", n0 = 4)
  expect_true(attr(sol, "converged"))
  expect_equal(as.numeric(sol$birth_rate), 10, tolerance = 1e-5)
})

test_that("the root finders solve a two-species fixed point", {
  map <- function(n) c(n[1] * exp(1.5 * (1 - n[1] / 10)),
                       n[2] * exp(1.5 * (1 - n[2] / 4)))
  for (solver in c("equilibrium_solve_nleqslv", "equilibrium_solve_dfsane")) {
    sol <- map_community(map, solver, n0 = c(3, 8))
    expect_true(attr(sol, "converged"), info = solver)
    expect_equal(as.numeric(sol$birth_rate), c(10, 4), tolerance = 1e-5,
                 info = solver)
  }
})

# ---- equilibrium_hybrid's extinct-species re-check --------------------------
#
# This is the block that was broken: it read `eq_solution$strategies` and called
# run_scm() on the result, treating a `community` as a plant `Parameters`. It
# now works off community$birth_rate and the harness connectors, so it runs on
# any model -- and, crucially, actually runs at all.

## Species 1 is a Ricker (fixed point 10); species 2 declines by the factor
## `growth_when_rare` when rare and 0.4x when common, so a root finder sends it
## to zero. If it grows when re-introduced, the solver killed a viable species.
two_species_map <- function(growth_when_rare) {
  function(n) {
    c(n[1] * exp(1.5 * (1 - n[1] / 10)),
      if (n[2] <= 1e-2) growth_when_rare * n[2] else 0.4 * n[2])
  }
}

hybrid_community <- function(growth_when_rare, n0 = c(3, 5), nattempts = 2) {
  community_start(bounds(x = c(-2, 2)),
                  harness = harness_map(two_species_map(growth_when_rare)),
                  trait_scale = "linear",
                  demography_control = demographic_step_control(list(
                    equilibrium_solver_name = "equilibrium_hybrid",
                    equilibrium_nsteps = 30,
                    equilibrium_nattempts = nattempts,
                    equilibrium_solver_logN = FALSE,
                    equilibrium_extinct_birth_rate = 1e-3))) |>
    community_add(trait_matrix(c(0, 1), "x"), birth_rate = n0) |>
    community_demography()
}

test_that("equilibrium_hybrid accepts a solution whose extinction is genuine", {
  # re-introduced at 1e-3 it shrinks to 5e-4, so the extinction is real
  sol <- hybrid_community(growth_when_rare = 0.5)
  expect_true(attr(sol, "converged"))
  expect_equal(as.numeric(sol$birth_rate[1]), 10, tolerance = 1e-5)
  expect_lt(as.numeric(sol$birth_rate[2]), 1e-3)
})

test_that("equilibrium_hybrid rejects a solution that killed a viable species", {
  # re-introduced at 1e-3 it doubles, so the solver was wrong to zero it; the
  # hybrid rejects every attempt and falls back on the iteration result, which
  # keeps species 2 alive. (The rejected attempts warn from the solvers -- that
  # is the speculative-solve path doing its job.)
  sol <- suppressWarnings(hybrid_community(growth_when_rare = 2))
  expect_gt(as.numeric(sol$birth_rate[2]), 1e-3)
})
