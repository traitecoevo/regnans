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
# The toy harness's demography runner is its one-generation map, so every
# solver here genuinely iterates (or root-finds) to the fixed point, as on the
# plant SCM; the "model" solver returns the model's own equilibrium and is the
# reference. A single step is just that: one generation.

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
  solvers <- c("model", "equilibrium_iteration", "equilibrium_hybrid",
               "equilibrium_solve_newton", "equilibrium_solve_nleqslv",
               "equilibrium_solve_dfsane")
  for (s in solvers) {
    comm <- solve_with(s)
    expect_true(attr(comm, "converged"), info = s)
    expect_equal(as.numeric(comm$birth_rate), target, tolerance = 1e-4, info = s)
    # at the equilibrium each resident's invasion fitness is ~0
    expect_equal(comm$resident_fitness, c(0, 0), tolerance = 1e-4, info = s)
    if (s != "model") expect_gt(NROW(attr(comm, "progress")), 1L, label = s)
  }
})

test_that("a single step is one generation, not the equilibrium", {
  target <- dd99_equilibrium(c(-0.5, 0.5), dd99_pars)
  one <- solve_with("single_step")
  expect_equal(NROW(attr(one, "progress")), 1L)
  expect_false(isTRUE(all.equal(as.numeric(one$birth_rate), target, tolerance = 1e-3)))
  # but it moves from the start (100, 100) towards it
  expect_lt(sum(abs(as.numeric(one$birth_rate) - target)), sum(abs(c(100, 100) - target)))
})

test_that("the nleqslv and dfsane solvers agree with the default iteration", {
  ref <- as.numeric(solve_with("equilibrium_iteration")$birth_rate)
  expect_equal(as.numeric(solve_with("equilibrium_solve_nleqslv")$birth_rate),
               ref, tolerance = 1e-5)
  expect_equal(as.numeric(solve_with("equilibrium_solve_dfsane")$birth_rate),
               ref, tolerance = 1e-5)
})

test_that("an unknown solver is rejected at the control, and again if it reaches dispatch", {
  expect_error(demographic_step_control(list(equilibrium_solver_name = "not_a_solver")),
               "should be one of")
  comm <- community_start(bounds(x = c(-2, 2)), harness = harness_dd99()) |>
    community_add(trait_matrix(0, "x"), birth_rate = 100)
  comm$demography_control$equilibrium_solver_name <- "not_a_solver"
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

  for (solver in c("equilibrium_solve_newton", "equilibrium_solve_nleqslv", "equilibrium_solve_dfsane")) {
    sol <- map_community(map, solver, n0 = 4)
    expect_true(attr(sol, "converged"), info = solver)
    expect_equal(as.numeric(sol$birth_rate), 10, tolerance = 1e-5, info = solver)
  }
  # Newton gets there in a handful of map evaluations where the iteration
  # needed more than its 30-step budget
  newton <- map_community(map, "equilibrium_solve_newton", n0 = 4)
  expect_lt(attr(newton, "n_calls"), 15L)
  expect_equal(dim(newton$demography_state$jacobian), c(1L, 1L))
})

test_that("a carried-over Jacobian saves the finite-difference pass on the next solve", {
  map <- function(n) c(n[1] * exp(1.5 * (1 - n[1] / 10)),
                       n[2] * exp(1.5 * (1 - n[2] / 4)))
  first <- map_community(map, "equilibrium_solve_newton", n0 = c(3, 8))
  expect_equal(dim(first$demography_state$jacobian), c(2L, 2L))
  again <- first
  again$birth_rate <- c(3, 8)
  again <- community_demography(again)
  expect_equal(as.numeric(again$birth_rate), c(10, 4), tolerance = 1e-5)
  # two finite-difference evaluations saved
  expect_lte(attr(again, "n_calls"), attr(first, "n_calls") - 2L)
  # a hint for a different number of residents is discarded, not misused
  one <- map_community(ricker_map(r = 1.5, K = 10), "equilibrium_solve_newton", n0 = 4)
  one$demography_state <- first$demography_state
  one$birth_rate <- 4
  one <- community_demography(one)
  expect_equal(as.numeric(one$birth_rate), 10, tolerance = 1e-5)
})

test_that("the carried-over Jacobian is keyed on the residual it describes", {
  map <- function(n) c(n[1] * exp(1.5 * (1 - n[1] / 10)),
                       n[2] * exp(1.5 * (1 - n[2] / 4)))
  first <- map_community(map, "equilibrium_solve_newton", n0 = c(3, 8))
  state <- first$demography_state
  expect_named(state, c("jacobian", "key"))
  expect_equal(state$key$i_keep, 1:2)
  expect_type(state$key$keep, "logical")
  expect_true(state$key$logN)

  # a kept species' residual is relative and a free one's absolute, so a hint
  # from solves whose keep flags differ is for another system: discarded, and
  # the solve costs the same as a cold one
  flipped <- first
  flipped$birth_rate <- c(3, 8)
  flipped$demography_state$key$keep <- !state$key$keep
  flipped <- community_demography(flipped)
  expect_equal(as.numeric(flipped$birth_rate), c(10, 4), tolerance = 1e-5)
  expect_equal(attr(flipped, "n_calls"), attr(first, "n_calls"))

  # changing the residents clears the hint; canonical_rhs() re-attaches it
  # deliberately after community_add() when it wants the warm start
  expect_null(community_add(first, trait_matrix(0, "x"), birth_rate = 2)$demography_state)
  expect_null(community_reset(first)$demography_state)
})

test_that("a failed Newton solve leaves no Jacobian behind", {
  # finite only at the starting point, so the finite-difference Jacobian is
  # NaN and so is every refresh: the solve fails rather than converging
  map <- function(n) if (all(abs(n - 4) < 1e-12)) 5 else NaN
  expect_warning(
    sol <- map_community(map, "equilibrium_solve_newton", n0 = 4),
    "did not converge")
  expect_false(attr(sol, "converged"))
  expect_null(sol$demography_state)
})

test_that("the hybrid solver reaches the fixed point the iteration missed", {
  sol <- map_community(ricker_map(r = 1.9, K = 10), "equilibrium_hybrid", n0 = 4)
  expect_true(attr(sol, "converged"))
  expect_equal(as.numeric(sol$birth_rate), 10, tolerance = 1e-5)
})

test_that("the root finders solve a two-species fixed point", {
  map <- function(n) c(n[1] * exp(1.5 * (1 - n[1] / 10)),
                       n[2] * exp(1.5 * (1 - n[2] / 4)))
  for (solver in c("equilibrium_solve_newton", "equilibrium_solve_nleqslv", "equilibrium_solve_dfsane")) {
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
