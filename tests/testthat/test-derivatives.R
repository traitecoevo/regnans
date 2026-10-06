# The derivative layer: every derivative in the package comes from here.
#
# DD99 oracles, with a single resident x at equilibrium density N = K(x):
#   s(y)   = r (1 - K(x) C(y - x) / K(y))
#   ds/dy  = -r (K(x)/K(y)) C(y - x) [ (y - x0)/sigma_K^2 - (y - x)/sigma_C^2 ]
#   G(x)   = ds/dy |_{y=x} = -r (x - x0) / sigma_K^2
#   dG/dx  = -r / sigma_K^2
#   d2s/dy2 |_{x = x0} = r (1/sigma_C^2 - 1/sigma_K^2)

dd99_resident <- function(x, x0 = 0, sigma_K = 1, sigma_C = 0.4, r = 1, ...) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  harness = harness_dd99(r = r, x0 = x0, sigma_K = sigma_K,
                                         sigma_C = sigma_C), ...) |>
    community_add(trait_matrix(x, "x"), birth_rate = 100) |>
    community_demography()
}

dd99_dsdy <- function(y, x, x0 = 0, sigma_K = 1, sigma_C = 0.4, r = 1) {
  K <- function(z) exp(-(z - x0)^2 / (2 * sigma_K^2))
  C <- exp(-(y - x)^2 / (2 * sigma_C^2))
  -r * (K(x) / K(y)) * C * ((y - x0) / sigma_K^2 - (y - x) / sigma_C^2)
}

# The multi-trait model has product kernels, so the same expression holds per
# dimension with the shared prefactor K(x) C(y - x) / K(y) taken over all of them.
dd99_nd_dsdy <- function(y, x, x0, sigma_K = c(1, 1), sigma_C = c(0.4, 0.4), r = 1) {
  t(apply(y, 1, function(yi) {
    K_ratio <- exp(-sum((x - x0)^2 / (2 * sigma_K^2)) + sum((yi - x0)^2 / (2 * sigma_K^2)))
    C <- exp(-sum((yi - x)^2 / (2 * sigma_C^2)))
    -r * K_ratio * C * ((yi - x0) / sigma_K^2 - (yi - x) / sigma_C^2)
  }))
}

test_that("derivative_control supplies defaults and validates overrides", {
  ctrl <- derivative_control()
  expect_equal(ctrl$d_gradient, 1e-4)
  expect_equal(ctrl$r_hessian, 2L)
  expect_equal(ctrl$r_jacobian, 1L)

  ctrl <- derivative_control(list(d_second = 1e-2, r_hessian = 3))
  expect_equal(ctrl$d_second, 1e-2)
  expect_identical(ctrl$r_hessian, 3L)
  expect_equal(ctrl$d_gradient, 1e-4)

  expect_error(derivative_control(list(dx = 1)), "Unknown control parameters dx")
  expect_error(derivative_control(list(d_gradient = -1)), "single positive number")
  expect_error(derivative_control(list(r_jacobian = 1.5)), "positive whole number")
  expect_error(derivative_control(list(r_hessian = 0)), "positive whole number")
})

test_that("community_start carries derivative settings", {
  comm <- community_start(bounds(x = c(-2, 2)), harness = harness_dd99())
  expect_equal(comm$derivative_control, derivative_control())

  comm <- community_start(bounds(x = c(-2, 2)), harness = harness_dd99(),
                          derivative_control = list(r_gradient = 2))
  expect_identical(comm$derivative_control$r_gradient, 2L)
})

test_that("community_fitness_gradient matches the DD99 slope at and away from the resident", {
  x <- 0.5
  comm <- dd99_resident(x)

  g <- community_fitness_gradient(comm)
  expect_equal(dim(g), c(1L, 1L))
  expect_equal(colnames(g), "x")
  expect_equal(as.numeric(g), -0.5, tolerance = 1e-6)
  expect_equal(as.numeric(attr(g, "value")), 0, tolerance = 1e-8)

  y <- c(-0.3, 0.2, 0.8, 1.4)
  g <- community_fitness_gradient(comm, y)
  expect_equal(dim(g), c(4L, 1L))
  expect_equal(as.numeric(g), dd99_dsdy(y, x), tolerance = 1e-6)
  expect_equal(as.numeric(attr(g, "value")), comm$fitness_function(y))
})

test_that("community_fitness_gradient handles several traits and several points", {
  comm <- community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)),
                          trait_scale = "linear",
                          harness = harness_dd99_nd(x0 = c(0.3, -0.5))) |>
    community_add(trait_matrix(c(0.5, 0.1), c("x1", "x2")), birth_rate = 100) |>
    community_demography()

  g <- community_fitness_gradient(comm)
  expect_equal(dim(g), c(1L, 2L))
  expect_equal(as.numeric(g), -(c(0.5, 0.1) - c(0.3, -0.5)), tolerance = 1e-6)

  g <- community_fitness_gradient(comm, c(0.5, 0.1))
  expect_equal(dim(g), c(1L, 2L))

  pts <- rbind(c(0.5, 0.1), c(0.3, -0.5), c(-1, 1))
  g <- community_fitness_gradient(comm, pts)
  expect_equal(dim(g), c(3L, 2L))
  expect_equal(unname(g[, ]), dd99_nd_dsdy(pts, c(0.5, 0.1), c(0.3, -0.5)),
               tolerance = 1e-6)
})

test_that("community_fitness_gradient refuses malformed input", {
  comm <- dd99_resident(0.5)
  expect_error(community_fitness_gradient(comm, matrix(1, 1, 2)),
               "Expected points with 1 trait column")
  expect_error(community_fitness_gradient(comm, numeric(0)),
               "No trait points")
  bare <- community_start(bounds(x = c(-2, 2)), harness = harness_dd99()) |>
    community_add(trait_matrix(0, "x"))
  expect_error(community_fitness_gradient(bare),
               "run community_demography\\(\\) first")
})

test_that("community_fitness_hessian matches the DD99 curvature", {
  H <- community_fitness_hessian(dd99_resident(0, sigma_C = 0.4))
  expect_equal(dim(H), c(1L, 1L))
  expect_equal(dimnames(H), list("x", "x"))
  expect_equal(as.numeric(H), 1 / 0.4^2 - 1, tolerance = 1e-4)

  expect_error(community_fitness_hessian(dd99_resident(0), c(0, 1)),
               "single trait point")
})

test_that("community_selection_gradient_jacobian matches the DD99 resident derivative", {
  comm <- dd99_resident(0.4, sigma_K = 1.5)
  J <- community_selection_gradient_jacobian(comm)
  expect_equal(dim(J), c(1L, 1L))
  expect_equal(as.numeric(J), -1 / 1.5^2, tolerance = 1e-4)
  g0 <- attr(J, "selection_gradient")
  expect_equal(names(g0), "x")
  expect_equal(as.numeric(g0), -0.4 / 1.5^2, tolerance = 1e-6)

  two <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                         harness = harness_dd99()) |>
    community_add(trait_matrix(c(-1, 1), "x"), birth_rate = c(100, 100)) |>
    community_demography()
  expect_error(community_selection_gradient_jacobian(two),
               "needs exactly one resident")
})

test_that("community_selection_gradient returns one row per resident", {
  one <- dd99_resident(0.5) |> community_selection_gradient()
  expect_equal(one$selection_gradient, -0.5, tolerance = 1e-6)
  expect_null(dim(one$selection_gradient))
  expect_length(one$resident_fitness, 1L)

  two <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                         harness = harness_dd99(sigma_C = 1.5)) |>
    community_add(trait_matrix(c(-0.6, 0.6), "x"), birth_rate = c(100, 100)) |>
    community_demography() |>
    community_selection_gradient()
  expect_equal(dim(two$selection_gradient), c(2L, 1L))
  expect_length(two$resident_fitness, 2L)
  expect_equal(as.numeric(two$resident_fitness), c(0, 0), tolerance = 1e-8)
  # a symmetric pair about x0 = 0 feels mirror-image selection
  expect_equal(two$selection_gradient[1, 1], -two$selection_gradient[2, 1],
               tolerance = 1e-6)
})

test_that("derivative settings on the community change the stencil", {
  coarse <- dd99_resident(0.5, derivative_control = list(d_gradient = 0.2,
                                                         eps_gradient = 0.2))
  fine <- dd99_resident(0.5)
  g_coarse <- as.numeric(community_fitness_gradient(coarse))
  g_fine <- as.numeric(community_fitness_gradient(fine))
  expect_gt(abs(g_coarse + 0.5), abs(g_fine + 0.5))
  expect_equal(g_coarse, -0.5, tolerance = 0.05)
})
