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
# community_viable_fitness, max_growth_rate, max_fitness -- are log-scale /
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
  mx <- max_fitness(comm, tol = 1e-8)
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
  mx <- max_fitness(comm)
  expect_named(mx, c("x1", "x2"))
  expect_equal(as.numeric(mx), c(0, 0), tolerance = 1e-6)
  expect_equal(attr(mx, "fitness"), 0, tolerance = 1e-10)
})

# ---- community_viable_fitness ----------------------------------------------
#
# Oracle: a model whose fitness is a quadratic bowl, s(y) = 1 - (y - c)' M (y - c),
# is positive inside an ellipse. The smallest box containing the ellipse has
# half-widths sqrt(diag(M^-1)); a slice through the centre along each axis would
# give 1 / sqrt(diag(M)) instead, so a tilted ellipse tells the two apart.

bowl_community <- function(M, centre, trait_bounds, scale = "linear",
                           log_traits = FALSE) {
  k <- length(centre)
  names <- if (k == 1L) "x" else paste0("x", seq_len(k))
  fitness <- function(x_mut, x_res, n_res, pars) {
    y <- matrix(x_mut, ncol = k)
    if (pars$log_traits) y <- log(y)
    d <- sweep(y, 2, pars$centre)
    1 - rowSums((d %*% pars$M) * d)
  }
  h <- harness_explicit(fitness, function(x_res, pars) numeric(0),
                        pars = list(M = M, centre = centre,
                                    log_traits = log_traits),
                        trait_names = names, label = "bowl")
  b <- matrix(rep(trait_bounds, each = k), nrow = k,
              dimnames = list(names, c("lower", "upper")))
  community_start(b, harness = h, trait_scale = scale)
}

test_that("community_viable_fitness finds the 1-D viable interval", {
  comm <- bowl_community(matrix(1 / 0.8^2), centre = 0.5,
                         trait_bounds = c(-3, 3))
  vb <- community_viable_fitness(comm)
  expect_equal(dimnames(vb), list("x", c("lower", "upper")))
  expect_equal(as.numeric(vb), c(0.5 - 0.8, 0.5 + 0.8), tolerance = 1e-3)
})

test_that("community_viable_fitness searches on a log trait scale", {
  # positive where |log(y) - log(0.2)| < 1.5
  comm <- bowl_community(matrix(1 / 1.5^2), centre = log(0.2),
                         trait_bounds = c(1e-3, 50), scale = "log",
                         log_traits = TRUE)
  vb <- community_viable_fitness(comm)
  expect_equal(as.numeric(vb), 0.2 * exp(c(-1.5, 1.5)), tolerance = 1e-3)
})

test_that("community_viable_fitness bounds a tilted 2-D region", {
  M <- solve(matrix(c(1, 0.8, 0.8, 1), 2) * 0.5^2)   # correlated, half-axes ~0.5
  centre <- c(0.3, -0.2)
  vb <- community_viable_fitness(bowl_community(M, centre, c(-3, 3)))
  half <- sqrt(diag(solve(M)))
  expect_equal(dimnames(vb), list(c("x1", "x2"), c("lower", "upper")))
  expect_equal(vb[, "lower"], centre - half, tolerance = 1e-3,
               ignore_attr = TRUE)
  expect_equal(vb[, "upper"], centre + half, tolerance = 1e-3,
               ignore_attr = TRUE)
  # strictly wider than the axis slices through the centre
  expect_true(all(vb[, "upper"] - centre > 1 / sqrt(diag(M)) + 0.05))
})

test_that("community_viable_fitness stops at the community's bounds", {
  comm <- bowl_community(matrix(1 / 2^2), centre = 0.5, trait_bounds = c(0, 1))
  expect_equal(as.numeric(community_viable_fitness(comm)), c(0, 1))
})

test_that("community_viable_fitness starts from the maximum when x is inviable", {
  comm <- bowl_community(matrix(1 / 0.8^2), centre = 0.5,
                         trait_bounds = c(-3, 3))
  vb <- community_viable_fitness(comm, x = 2.5)
  expect_equal(as.numeric(vb), c(-0.3, 1.3), tolerance = 1e-3)
})

test_that("community_viable_fitness returns NULL where nothing is viable", {
  comm <- bowl_community(matrix(1 / 0.5^2), centre = 5, trait_bounds = c(-1, 1))
  expect_null(community_viable_fitness(comm))
})
