test_that("positive_1d", {
  f <- function(x) -x^2 + 1
  tol <- 1e-8
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  tol <- 1e-3
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  expect_error(positive_1d(f, -2, 0.1, tol=tol), "no positive values")
})

test_that("bounds", {

  # First test object with infinite bounds
  expect_silent({
    bounds0 <- bounds_infinite("lma")
  })
  expect_true(is.matrix(bounds0))
  expect_equal(bounds0,
    matrix(c(-Inf, Inf), nrow=1, dimnames =list("lma",c("lower", "upper")))
  )

  # Manually created bounds

  expect_silent(
    bounds_1d <- bounds(lma=c(0.01, 10))
  )
  expect_silent(
    bounds_2d <- bounds(lma=c(0.01, 10), rho=c(1, 1000))
  )

  expect_true(is.matrix(bounds_1d))
  expect_equal(
    bounds_1d,
    matrix(c(0.01, 10), nrow = 1, dimnames = list("lma", c("lower", "upper")))
  )

  expect_true(is.matrix(bounds_2d))
  expect_equal(
    bounds_2d,
    matrix(c(0.01, 10, 1, 1000), byrow=TRUE, nrow = 2, dimnames = list(c("lma", "rho"), c("lower", "upper")))
  )

  expect_silent(
    check_bounds(bounds_1d)
  )
  expect_silent(
    check_bounds(bounds_2d)
  )

  expect_error(
    check_bounds(c(0,1))
  )
  expect_error(
    check_bounds(matrix(c(0.01, 10, 1, 1000)))
  )
  expect_error(
    check_bounds(bounds0, finite=TRUE)
  )
  expect_silent(
    check_bounds(bounds0)
  )


  # Points lie within bounds
  expect_silent(
    check_point(0.02, bounds_1d)
  )
  expect_error(
    check_point(0.001, bounds_1d)
  )
})

# NOTE (issue #27): the SCM-backed tests that used to live here --
# community_viable_fitness_1D, max_growth_rate, max_fitness -- are log-scale /
# fundamental-niche behaviours specific to the plant path, so they now live in
# test-plant-smoke.R (the consolidated, minimal set of genuine SCM tests). The
# model-agnostic pipeline is covered fast on DD99 elsewhere.

# ---- max_fitness -----------------------------------------------------------
#
# Oracle: with sigma_C > sigma_K, DD99's singular strategy x0 is an ESS, so a
# resident at x0 at equilibrium is the fitness maximum, where fitness is 0. The
# bounds are asymmetric so the search does not start on the answer.

test_that("max_fitness finds the DD99 ESS in one trait", {
  comm <- community_start(bounds(x = c(-1, 3)), trait_scale = "linear",
                          harness = harness_dd99(x0 = 0, sigma_K = 1,
                                                 sigma_C = 1.5)) |>
    community_add(trait_matrix(0, "x"), birth_rate = 500) |>
    community_demography()
  mx <- max_fitness(comm, log_scale = FALSE, tol = 1e-8)
  expect_named(mx, "x")
  expect_equal(as.numeric(mx), 0, tolerance = 1e-6)
  expect_equal(attr(mx, "fitness"), 0, tolerance = 1e-10)
})

test_that("max_fitness finds the DD99 ESS in two traits", {
  h <- harness_dd99_nd(x0 = c(0, 0), sigma_K = c(1, 1), sigma_C = c(1.5, 1.5))
  comm <- community_start(bounds(x1 = c(-1, 3), x2 = c(-1, 3)),
                          trait_scale = "linear", harness = h) |>
    community_add(trait_matrix(c(0, 0), c("x1", "x2")), birth_rate = 500) |>
    community_demography()
  mx <- max_fitness(comm, log_scale = FALSE)
  expect_named(mx, c("x1", "x2"))
  expect_equal(as.numeric(mx), c(0, 0), tolerance = 1e-6)
  expect_equal(attr(mx, "fitness"), 0, tolerance = 1e-10)
})
