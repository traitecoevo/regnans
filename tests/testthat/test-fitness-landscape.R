# community_fitness_landscape, grid method (issue #27). The landscape machinery
# is model-agnostic, so it is exercised on the fast DD99 harness rather than the
# SCM. (The bayesopt/surrogate method is not covered here -- it pulls in
# mlr3mbo and is exercised separately.)

dd99_resident <- function(...) {
  community_start(bounds(x = c(-2, 2)),
                  harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4),
                  trait_scale = "linear",
                  fitness_control = list(method = "grid", ...))
}

test_that("community_fitness_landscape (grid) samples fitness over the bounds", {
  comm <- dd99_resident(n_evals = 12) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography() |>
    community_fitness_landscape()

  pts <- comm$fitness_points
  expect_s3_class(pts, "tbl_df")
  # first column is named by the trait, plus fitness + resident flag
  expect_equal(names(pts), c("x", "fitness", "resident"))

  # grid spans the bounds (augmented by the resident +/- offset points)
  expect_gte(nrow(pts), 12L)
  expect_gte(min(pts$x), -2)
  expect_lte(max(pts$x), 2)
  expect_true(all(is.finite(pts$fitness)))

  # the resident is flagged and (at equilibrium) has fitness ~0
  expect_true(any(pts$resident))
  expect_equal(pts$x[pts$resident], 0.5, tolerance = 1e-8)
  expect_lt(abs(pts$fitness[pts$resident]), 1e-6)
})

test_that("community_fitness_landscape solves demography if needed", {
  # No community_demography() call: the landscape function should solve it.
  comm <- community_fitness_landscape(dd99_resident(n_evals = 8))
  expect_true(is.function(comm$fitness_function))
  expect_s3_class(comm$fitness_points, "tbl_df")
})

test_that("community_fitness_landscape rejects an unknown method", {
  comm <- dd99_resident(n_evals = 8) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography()
  expect_error(community_fitness_landscape(comm, method = "nope"),
               "Unknown fitness landscape method")
})

# ---- fitness_control defaults ----------------------------------------------
#
# community_start() used to leave fitness_control NULL, so
# community_fitness_landscape() failed on `fitness_control$method` unless the
# caller happened to know to supply one.

test_that("fitness_landscape_control fills in defaults and rejects unknowns", {
  ctrl <- fitness_landscape_control()
  expect_equal(ctrl$method, "grid")
  expect_true(is.numeric(ctrl$n_evals) && ctrl$n_evals > 1)
  expect_true(is.numeric(ctrl$n_init))

  expect_equal(fitness_landscape_control(list(n_evals = 7))$n_evals, 7)
  # overriding one option leaves the others at their defaults
  expect_equal(fitness_landscape_control(list(n_evals = 7))$method, "grid")
  expect_error(fitness_landscape_control(list(nope = 1)),
               "Unknown fitness control parameters")
})

test_that("community_start supplies a working fitness_control by default", {
  comm <- community_start(bounds(x = c(-2, 2)),
                          harness = harness_dd99(x0 = 0, sigma_K = 1,
                                                 sigma_C = 0.4),
                          trait_scale = "linear")
  expect_equal(comm$fitness_control$method, "grid")

  # the whole point: this used to fail with a NULL `method`
  out <- comm |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography() |>
    community_fitness_landscape()
  expect_s3_class(out$fitness_points, "tbl_df")
  expect_gte(nrow(out$fitness_points), fitness_landscape_control()$n_evals)
})

test_that("a community built without going through community_start still works", {
  comm <- dd99_resident(n_evals = 8) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography()
  comm$fitness_control <- NULL
  out <- community_fitness_landscape(comm)
  expect_equal(out$fitness_control$method, "grid")
  expect_s3_class(out$fitness_points, "tbl_df")
})
