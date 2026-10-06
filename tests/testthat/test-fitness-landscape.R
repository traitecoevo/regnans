# community_fitness_landscape, grid and bayesopt methods (issue #27). The
# landscape machinery is model-agnostic, so it is exercised on the fast DD99
# harness rather than the SCM.

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

# ---- bayesopt --------------------------------------------------------------
#
# The Gaussian-process landscape samples the true fitness function where the
# expected improvement is highest and fits a surrogate through the samples. The
# oracle for the samples is the fitness function itself; the surrogate (a GP
# without a nugget) must interpolate them.

skip_bayesopt <- function() {
  for (pkg in c("mlr3mbo", "DiceKriging", "nloptr", "lgr", "withr")) {
    skip_if_not_installed(pkg)
  }
}

quiet_bbotk <- function(env = parent.frame()) {
  logger <- lgr::get_logger("mlr3/bbotk")
  old <- logger$threshold
  logger$set_threshold("warn")
  withr::defer(logger$set_threshold(old), envir = env)
}

dd99_bayesopt <- function(trait_bounds = c(-2, 2), x0 = 0, scale = "linear",
                          resident = -0.5) {
  community_start(bounds(x = trait_bounds), trait_scale = scale,
                  harness = harness_dd99(x0 = x0, sigma_K = 1, sigma_C = 0.4),
                  fitness_control = list(method = "bayesopt", n_evals = 12,
                                         n_init = 6)) |>
    community_add(trait_matrix(resident, "x"), birth_rate = 100) |>
    community_demography()
}

test_that("bayesopt samples true fitness on a linear trait scale", {
  skip_bayesopt(); quiet_bbotk()
  comm <- dd99_bayesopt()
  set.seed(1)
  expect_no_warning(out <- community_fitness_landscape(comm))

  pts <- out$fitness_points
  expect_equal(names(pts), c("x", "fitness", "batch_nr", "resident"))
  expect_equal(nrow(pts), 12L)
  expect_false(is.unsorted(pts$x))
  # negative trait values: only reachable now the search follows trait_scale
  expect_true(all(pts$x >= -2 & pts$x <= 2))
  expect_lt(min(pts$x), 0)
  expect_equal(pts$fitness, out$fitness_function(pts$x), tolerance = 1e-12)
  # searched beyond the initial design
  expect_gt(max(pts$batch_nr), 1L)

  expect_equal(sum(pts$resident), 1L)
  expect_equal(pts$x[pts$resident], -0.5, tolerance = 1e-12)
  expect_lt(abs(pts$fitness[pts$resident]), 1e-6)

  expect_equal(out$fitness_surrogate_function(pts$x), pts$fitness,
               tolerance = 1e-6)
})

test_that("bayesopt searches on a log trait scale and flags the resident", {
  skip_bayesopt(); quiet_bbotk()
  comm <- dd99_bayesopt(trait_bounds = c(0.1, 3), x0 = 1, scale = "log",
                        resident = 0.8)
  set.seed(1)
  pts <- community_fitness_landscape(comm)$fitness_points
  # exp(log(3)) may land an ulp past the bound
  expect_true(all(pts$x >= 0.1 * (1 - 1e-12) & pts$x <= 3 * (1 + 1e-12)))
  expect_equal(pts$fitness, comm$fitness_function(pts$x), tolerance = 1e-12)
  # flagged on the search scale, so the exp(log(x)) round trip cannot lose it
  expect_equal(sum(pts$resident), 1L)
  expect_equal(pts$x[pts$resident], 0.8, tolerance = 1e-12)
})

test_that("bayesopt honours the bounds argument", {
  skip_bayesopt(); quiet_bbotk()
  set.seed(1)
  out <- community_fitness_landscape(dd99_bayesopt(),
                                     bounds = bounds(x = c(-1, 0)))
  pts <- out$fitness_points
  expect_equal(range(pts$x), c(-1, 0))
})

test_that("bayesopt is reproducible from a seed set by the caller", {
  skip_bayesopt(); quiet_bbotk()
  comm <- dd99_bayesopt()
  set.seed(7)
  a <- community_fitness_landscape(comm)$fitness_points
  set.seed(7)
  b <- community_fitness_landscape(comm)$fitness_points
  expect_identical(a, b)
})

test_that("a landscape refuses bounds that are not finite on the trait scale", {
  comm <- dd99_resident(n_evals = 8) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography()
  expect_error(community_fitness_landscape(comm,
                                           bounds = bounds(x = c(-Inf, 2))),
               "finite bounds")
})
