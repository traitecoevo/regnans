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

test_that("derivative settings on the community change the finite-difference stencil", {
  # GK98 supplies no derivatives, so its gradient is a finite difference
  gk98_at <- function(...) {
    community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                    harness = harness_gk98(d = 1.5), ...) |>
      community_add(trait_matrix(0.7, "x"), birth_rate = 1) |>
      community_demography()
  }
  fine <- gk98_at()
  coarse <- gk98_at(derivative_control = list(d_gradient = 0.3, eps_gradient = 0.3))
  reference <- as.numeric(regnans:::fd_fitness_gradient(
    fine$fitness_function, matrix(0.7), derivative_control(list(r_gradient = 4))))
  g_fine <- as.numeric(community_fitness_gradient(fine))
  g_coarse <- as.numeric(community_fitness_gradient(coarse))
  expect_gt(abs(g_coarse - reference), abs(g_fine - reference))
  expect_equal(g_fine, reference, tolerance = 1e-6)
  expect_equal(g_coarse, reference, tolerance = 0.1)
})

# ---- model-supplied derivatives and the finite-difference fill --------------

test_that("harnesses advertise what they can differentiate", {
  expect_equal(harness_provides(harness_dd99()),
               c("fitness_gradient", "fitness_hessian"))
  expect_true(harness_provides(harness_dd99_nd(), "fitness_gradient"))
  expect_equal(harness_provides(harness_gk98()), character(0))
  expect_false(harness_provides(harness_jj12(), "fitness_hessian"))
  expect_error(harness_provides(harness_dd99(), "jacobian"), "should be one of")
})

test_that("community_start records the source of every derivative", {
  dd <- dd99_resident(0.5)
  expect_equal(dd$fitness_derivatives$source,
               c(fitness_gradient = "model", fitness_hessian = "model"))

  gk <- community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                        harness = harness_gk98()) |>
    community_add(trait_matrix(0.2, "x"), birth_rate = 1) |>
    community_demography()
  expect_equal(gk$fitness_derivatives$source,
               c(fitness_gradient = "finite difference",
                 fitness_hessian = "finite difference"))
  expect_true(is.function(gk$fitness_derivatives$gradient))
  expect_true(is.function(gk$fitness_derivatives$hessian))
})

test_that("harness_fd is idempotent and required", {
  h <- harness_fd(harness_gk98())
  expect_true(h$fd)
  expect_identical(harness_fd(h), h)
  expect_error(harness_fd(list()), "needs a harness object")

  bare <- list(bounds = bounds(x = c(-4, 4)), trait_names = "x",
               traits = trait_matrix(0.2, "x"), birth_rate = 1,
               demography_control = demographic_step_control(),
               harness = harness_gk98(), trait_scale = "linear")
  class(bare) <- "community"
  bare <- community_demography(bare)
  expect_error(community_fitness_gradient(bare), "harness_fd")
})

test_that("DD99 model derivatives agree with finite differences and are exact at the resident", {
  comm <- dd99_resident(0.5)
  y <- c(-0.7, 0.1, 0.5, 1.2)
  g_model <- community_fitness_gradient(comm, y)
  g_fd <- regnans:::fd_fitness_gradient(comm$fitness_function,
                                        matrix(y, ncol = 1), derivative_control())
  expect_equal(unname(g_model[, ]), unname(g_fd[, ]), tolerance = 1e-6)
  expect_equal(as.numeric(g_model), dd99_dsdy(y, 0.5), tolerance = 1e-12)

  sg <- community_selection_gradient(comm)
  expect_equal(sg$selection_gradient, -0.5, tolerance = 1e-12)

  H <- community_fitness_hessian(dd99_resident(0, sigma_C = 0.4))
  expect_equal(as.numeric(H), 1 / 0.4^2 - 1, tolerance = 1e-12)
})

test_that("DD99 nD derivatives match the analytic slope and curvature", {
  x <- c(0.5, 0.1); x0 <- c(0.3, -0.5)
  comm <- community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)),
                          trait_scale = "linear",
                          harness = harness_dd99_nd(x0 = x0)) |>
    community_add(trait_matrix(x, c("x1", "x2")), birth_rate = 100) |>
    community_demography()
  pts <- rbind(c(0.5, 0.1), c(0.3, -0.5), c(-1, 1))
  expect_equal(unname(community_fitness_gradient(comm, pts)[, ]),
               dd99_nd_dsdy(pts, x, x0), tolerance = 1e-12)

  at <- community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)),
                        trait_scale = "linear",
                        harness = harness_dd99_nd(sigma_C = c(0.4, 1.5))) |>
    community_add(trait_matrix(c(0, 0), c("x1", "x2")), birth_rate = 100) |>
    community_demography()
  expect_equal(unname(community_fitness_hessian(at)),
               diag(c(1 / 0.4^2 - 1, 1 / 1.5^2 - 1)), tolerance = 1e-12)
})

test_that("a provider returning the wrong shape is refused", {
  h <- harness_explicit(
    fitness = dd99_fitness,
    equilibrium = function(x_res, pars) dd99_equilibrium(x_res, pars),
    pars = list(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4),
    trait_names = "x", label = "bad",
    fitness_gradient = function(x_mut, x_res, n_res, pars) cbind(x_mut, x_mut))
  comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                          harness = h) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography()
  expect_error(community_fitness_gradient(comm), "must return a 1 by 1 matrix")
})

test_that("print.harness reports derivatives", {
  out <- paste(utils::capture.output(print(harness_dd99())), collapse = "\n")
  expect_match(out, "derivatives: fitness_gradient, fitness_hessian")
  out <- paste(utils::capture.output(print(harness_fd(harness_gk98()))), collapse = "\n")
  expect_match(out, "derivatives: none \\(finite differences for the rest\\)")
})
