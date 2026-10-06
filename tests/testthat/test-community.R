# Tests for the community interface (community_start/add/drop, demography,
# selection gradient). These target regnans's own API, unlike the older
# test files which call plant's (now removed) internal functions directly.

# `assembly_model_support()` is defined in helper-assembly.R.

# ---- fast unit tests (no SCM runs) ----------------------------------------

test_that("trait_matrix builds a named single-trait matrix", {
  m <- trait_matrix(c(0.1, 0.2), "lma")
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 1L))
  expect_equal(colnames(m), "lma")
  expect_equal(as.numeric(m), c(0.1, 0.2))
})

test_that("bounds() produces lower/upper rows named by trait", {
  b <- bounds(lma = c(0.01, 2))
  expect_equal(rownames(b), "lma")
  expect_equal(unname(b["lma", ]), c(0.01, 2))
})

test_that("demographic_step_control has expected defaults and rejects unknowns", {
  ctrl <- demographic_step_control()
  expect_null(ctrl$equilibrium_solver_name)      # resolved per harness by community_start()
  for (nm in c("equilibrium_eps", "equilibrium_nsteps",
               "equilibrium_large_birth_rate_change",
               "equilibrium_extinct_birth_rate",
               "equilibrium_min_offspring_arriving")) {
    expect_true(nm %in% names(ctrl), info = nm)
  }
  expect_equal(demographic_step_control(list(equilibrium_eps = 1e-3))$equilibrium_eps,
               1e-3)
  expect_error(demographic_step_control(list(not_a_param = 1)),
               "Unknown control parameters")
})

test_that("community_start creates an empty community", {
  comm <- community_start(bounds(lma = c(0.01, 2)))
  expect_s3_class(comm, "community")
  expect_equal(comm$trait_names, "lma")
  expect_equal(length(comm), 0L)
  expect_equal(nrow(comm$traits), 0L)
})

test_that("community_add appends species and validates input", {
  comm <- community_start(bounds(lma = c(0.01, 2)))
  comm <- community_add(comm, trait_matrix(c(0.05, 0.1), "lma"), birth_rate = 5)
  expect_equal(length(comm), 2L)
  expect_equal(comm$birth_rate, c(5, 5))                 # recycled scalar
  expect_error(community_add(comm, c(0.05, 0.1)), "must be a matrix")
  expect_error(community_add(comm, trait_matrix(0.05, "lma"), birth_rate = c(1, 2)),
               "Incompatible length")
  expect_error(community_add(comm, trait_matrix(c(0.05, 0.1), c("lma", "x"))),
               "Incorrect size")
})

test_that("community_drop validates indices (no solve needed)", {
  # Index validation happens before any demography re-solve, so no model_support.
  comm <- community_start(bounds(lma = c(0.01, 2)))
  comm <- community_add(comm, trait_matrix(c(0.05, 0.1, 0.2), "lma"),
                        birth_rate = c(1, 2, 3))
  expect_error(community_drop(comm, 99), "Invalid indicies")
  expect_error(community_drop(comm, c(TRUE, FALSE)), "Invalid length")
})

test_that("community_drop keeps the correct species (regression: numeric index)", {
  # community_drop resets and re-solves the survivors, so it needs a model.
  # Run it on the fast DD99 harness (issue #27) -- the drop/re-solve logic is
  # model-agnostic, and DD99 traits are stable under the re-solve (only the
  # equilibrium density changes), so we assert on which traits survive.
  comm <- community_start(bounds(x = c(-2, 2)), harness = harness_dd99()) |>
    community_add(trait_matrix(c(-0.5, 0.5), "x"), birth_rate = c(50, 50))

  d <- community_drop(comm, 2)            # drop the 2nd species
  expect_equal(length(d), 1L)
  expect_equal(as.numeric(d$traits[, "x"]), -0.5)

  d2 <- community_drop(comm, c(FALSE, TRUE))  # logical mask, drop the 2nd
  expect_equal(as.numeric(d2$traits[, "x"]), -0.5)
})

# ---- regression: schedule defaults track max_patch_lifetime ----------------

test_that("plant_default_assembly_pars regenerates the node schedule default", {
  # Lowering max_patch_lifetime must shrink node_schedule_times_default to match,
  # otherwise the SCM errors with "time_max must be greater than current time"
  # when the equilibrium runner resets the schedule.
  p <- plant_default_assembly_pars(max_patch_lifetime = 30)
  expect_equal(p$max_patch_lifetime, 30)
  expect_lte(max(p$node_schedule_times_default), 30)
})

# ---- pipeline on the fast DD99 harness (issue #27) -------------------------
# The demography / selection-gradient pipeline is model-agnostic, so it is
# exercised here against DD99 (instant, with analytic oracles) rather than the
# SCM. The genuine plant-SCM equivalents (equilibrium birth rate, fundamental
# niche, attractor) live in test-plant-smoke.R.

test_that("community_demography on an empty community gives no resident fitness", {
  comm <- community_start(bounds(x = c(-2, 2)), harness = harness_dd99()) |>
    community_demography()
  expect_length(comm$resident_fitness, 0L)
  expect_true(is.function(comm$fitness_function))
})

test_that("community_demography solves a single resident to its analytic equilibrium", {
  comm <- community_start(bounds(x = c(-2, 2)),
                          harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4)) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 200) |>
    community_demography()

  expect_true(attr(comm, "converged"))
  expect_length(comm$birth_rate, 1L)
  # DD99 single-resident equilibrium is the carrying capacity K(x)
  Kx <- 500 * exp(-(0.5)^2 / (2 * 1^2))
  expect_equal(as.numeric(comm$birth_rate), Kx, tolerance = 1e-8)
  expect_equal(comm$resident_fitness, 0, tolerance = 1e-8)  # ~0 at equilibrium
})

test_that("community_selection_gradient recovers the analytic DD99 gradient", {
  # DD99 selection gradient is -(x - x0) / sigma_K^2; at x = 0.5 (x0 = 0,
  # sigma_K = 1) that is exactly -0.5, and it vanishes at the singular point x0.
  grad_at <- function(x) {
    community_start(bounds(x = c(-2, 2)),
                    harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4)) |>
      community_add(trait_matrix(x, "x"), birth_rate = 200) |>
      community_demography() |>
      community_selection_gradient()
  }
  g <- grad_at(0.5)
  expect_length(g$selection_gradient, 1L)
  expect_equal(g$selection_gradient, -0.5, tolerance = 1e-4)
  expect_equal(grad_at(0)$selection_gradient, 0, tolerance = 1e-6)
})
