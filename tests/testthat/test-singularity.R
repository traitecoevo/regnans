# N-dimensional singular-strategy solving and classification.
#
# Both functions are model-agnostic, so they are developed and validated here
# against the toy harnesses, whose singular strategies and second-order
# conditions are known analytically (the plant SCM anchor is in
# test-plant-smoke-singularity.R). The oracles used below:
#
# DD99, single resident at the singular strategy x* = x0, N = K(x0):
#   s(y) = r (1 - exp(-(1/sigma_C^2 - 1/sigma_K^2) y^2 / 2))   (measuring y from x0)
#   => d2s/dy2 = r (1/sigma_C^2 - 1/sigma_K^2)                  evolutionary stability
#      G(x)    = r dlogK/dx = -r (x - x0)/sigma_K^2
#   => dG/dx   = -r / sigma_K^2                                 convergence stability
# so DD99 is always convergence stable, and is a branching point exactly when
# sigma_C < sigma_K -- the documented oracle.
#
# GK98, symmetric three-patch (mu = (-d, 0, d), equal capacities, x* = 0):
#   S(y) = -y^2/(2 sigma^2) + log((1 + 2 cosh(y d / sigma^2))/3)
#   => d2S/dy2 = -1/sigma^2 + 2 d^2 / (3 sigma^4)
# which crosses zero at d/sigma = sqrt(3/2), matching the branching threshold.
#
# JJ12: a continuously stable strategy with the closed form x* = x_opt - a sigma^2.

dd99_1d <- function(sigma_C = 0.4, sigma_K = 1, x0 = 0, r = 1) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  harness = harness_dd99(r = r, x0 = x0,
                                         sigma_K = sigma_K, sigma_C = sigma_C))
}

dd99_2d <- function(sigma_C = c(0.4, 1.5), sigma_K = c(1, 1), x0 = c(0, 0)) {
  community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)), trait_scale = "linear",
                  harness = harness_dd99_nd(x0 = x0, sigma_K = sigma_K,
                                            sigma_C = sigma_C))
}

# ---- community_solve_singularity --------------------------------------------

test_that("community_solve_singularity recovers the 1D DD99 singular strategy", {
  out <- community_solve_singularity(dd99_1d(x0 = 0.6))
  expect_true(attr(out, "converged"))
  expect_equal(unname(attr(out, "singularity")), 0.6, tolerance = 1e-6)
  expect_equal(names(attr(out, "singularity")), "x")
  # the returned community is the one *at* the root, with a vanishing gradient
  expect_equal(as.numeric(out$traits), 0.6, tolerance = 1e-6)
  expect_lt(abs(out$selection_gradient), 1e-6)
})

test_that("community_solve_singularity agrees with the 1D uniroot solver", {
  nd <- community_solve_singularity(dd99_1d(x0 = -0.35))
  one <- community_solve_singularity_1D(dd99_1d(x0 = -0.35), tol = 1e-8)
  expect_equal(as.numeric(nd$traits), as.numeric(one$traits), tolerance = 1e-5)
})

test_that("community_solve_singularity recovers the 2-trait DD99 singularity", {
  x0 <- c(0.3, -0.5)
  for (solver in c("nleqslv", "dfsane")) {
    out <- community_solve_singularity(dd99_2d(x0 = x0), solver = solver)
    expect_true(attr(out, "converged"), info = solver)
    expect_equal(unname(attr(out, "singularity")), x0, tolerance = 1e-5,
                 info = solver)
    expect_equal(names(attr(out, "singularity")), c("x1", "x2"), info = solver)
    expect_equal(attr(out, "solver"), solver)
    expect_lt(max(abs(out$selection_gradient)), 1e-5)
  }
})

test_that("community_solve_singularity works on the JJ12 closed form", {
  a <- 0.1; sigma <- 1.4; x_opt <- 0.5
  out <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                         harness = harness_jj12(a = a, x_opt = x_opt,
                                                sigma = sigma)) |>
    community_solve_singularity()
  expect_equal(unname(attr(out, "singularity")), x_opt - a * sigma^2,
               tolerance = 1e-6)
})

test_that("community_solve_singularity honours the trait scale", {
  # a log-scale community searches in log(x); the root is the same
  h <- harness_dd99(x0 = 0.5, sigma_K = 1, sigma_C = 0.4, trait_name = "lma")
  out <- community_start(bounds(lma = c(0.05, 2)), trait_scale = "log",
                         harness = h) |>
    community_solve_singularity()
  expect_equal(unname(attr(out, "singularity")), 0.5, tolerance = 1e-5)
})

test_that("community_solve_singularity warns when the bounds exclude the root", {
  expect_warning(
    out <- community_solve_singularity(dd99_1d(x0 = 0), bounds = c(0.5, 1.5)),
    "Bounds do not include a singularity")
  expect_equal(as.numeric(out$traits), 0.5, tolerance = 1e-8)

  expect_error(
    suppressWarnings(
      community_solve_singularity(dd99_1d(x0 = 0), bounds = c(0.5, 1.5),
                                  edge_ok = FALSE)),
    "Bounds do not include a singularity")
})

test_that("community_solve_singularity validates x0 and bounds", {
  expect_error(community_solve_singularity(dd99_2d(), x0 = 1),
               "x0 must have one value per trait")
  expect_error(community_solve_singularity(dd99_2d(), bounds = c(-1, 1)),
               "bounds must be a 2 x 2 matrix")
})

# ---- community_classify_singularity: 1-D ------------------------------------

test_that("classification of DD99 matches the analytic second-order conditions", {
  # sigma_C < sigma_K: fitness minimum -> branching point
  cl <- community_solve_singularity(dd99_1d(sigma_C = 0.4, sigma_K = 1)) |>
    community_classify_singularity()

  expect_s3_class(cl, "singularity_classification")
  expect_equal(as.numeric(cl$hessian), 1 / 0.4^2 - 1 / 1^2, tolerance = 1e-4)
  expect_equal(as.numeric(cl$jacobian), -1 / 1^2, tolerance = 1e-4)
  expect_false(cl$evolutionarily_stable)
  expect_true(cl$convergence_stable)
  expect_true(cl$strongly_convergence_stable)
  expect_equal(cl$classification, "branching point")
  expect_equal(unname(cl$branching_direction), 1)
  expect_equal(dim(cl$hessian), c(1L, 1L))
  expect_equal(dimnames(cl$hessian), list("x", "x"))
})

test_that("DD99 with a wide competition kernel classifies as a CSS", {
  cl <- community_solve_singularity(dd99_1d(sigma_C = 1.5, sigma_K = 1)) |>
    community_classify_singularity()
  expect_equal(as.numeric(cl$hessian), 1 / 1.5^2 - 1, tolerance = 1e-4)
  expect_true(cl$evolutionarily_stable)
  expect_true(cl$convergence_stable)
  expect_equal(cl$classification, "CSS")
  expect_null(cl$branching_direction)
})

test_that("GK98 classification tracks the d/sigma = sqrt(3/2) threshold", {
  gk98_H <- function(d, sigma = 1) -1 / sigma^2 + 2 * d^2 / (3 * sigma^4)

  css <- community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                         harness = harness_gk98(d = 1.0, sigma = 1)) |>
    community_solve_singularity() |>
    community_classify_singularity()
  expect_equal(as.numeric(css$hessian), gk98_H(1.0), tolerance = 1e-4)
  expect_true(css$evolutionarily_stable)
  expect_equal(css$classification, "CSS")

  branch <- community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                            harness = harness_gk98(d = 1.5, sigma = 1)) |>
    community_solve_singularity() |>
    community_classify_singularity()
  expect_equal(as.numeric(branch$hessian), gk98_H(1.5), tolerance = 1e-4)
  expect_false(branch$evolutionarily_stable)
  expect_true(branch$convergence_stable)
  expect_equal(branch$classification, "branching point")
})

test_that("JJ12 classifies as a CSS", {
  cl <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                        harness = harness_jj12(a = 0.1, x_opt = 0,
                                               sigma = 1)) |>
    community_solve_singularity() |>
    community_classify_singularity()
  expect_equal(unname(cl$traits), -0.1, tolerance = 1e-6)
  expect_lt(as.numeric(cl$hessian), 0)
  expect_lt(as.numeric(cl$jacobian), 0)
  expect_equal(cl$classification, "CSS")
})

# ---- community_classify_singularity: N-D ------------------------------------

test_that("2-trait DD99 classification recovers the branching direction", {
  # dimension 1 branches (sigma_C < sigma_K), dimension 2 is an ESS direction
  cl <- community_solve_singularity(
    dd99_2d(sigma_C = c(0.4, 1.5), sigma_K = c(1, 1))) |>
    community_classify_singularity()

  expect_equal(dim(cl$hessian), c(2L, 2L))
  expect_equal(dimnames(cl$hessian), list(c("x1", "x2"), c("x1", "x2")))
  # product-Gaussian kernels -> diagonal Hessian, one entry per dimension
  expect_equal(unname(diag(cl$hessian)),
               c(1 / 0.4^2 - 1, 1 / 1.5^2 - 1), tolerance = 1e-4)
  expect_lt(max(abs(cl$hessian[upper.tri(cl$hessian)])), 1e-6)
  expect_equal(unname(diag(cl$jacobian)), c(-1, -1), tolerance = 1e-4)

  expect_false(cl$evolutionarily_stable)
  expect_true(cl$convergence_stable)
  expect_true(cl$strongly_convergence_stable)
  expect_equal(cl$classification, "branching point")

  # the eigen-decomposition, not just the verdict: the unstable direction is
  # trait 1 alone
  expect_equal(sort(cl$hessian_eigen$values),
               sort(c(1 / 0.4^2 - 1, 1 / 1.5^2 - 1)), tolerance = 1e-4)
  expect_equal(unname(abs(cl$branching_direction)), c(1, 0), tolerance = 1e-6)
  expect_equal(names(cl$branching_direction), c("x1", "x2"))
})

test_that("2-trait DD99 with both kernels wide is a CSS in every direction", {
  cl <- community_solve_singularity(
    dd99_2d(sigma_C = c(1.5, 2), sigma_K = c(1, 1))) |>
    community_classify_singularity()
  expect_true(all(cl$hessian_eigen$values < 0))
  expect_true(cl$evolutionarily_stable)
  expect_equal(cl$classification, "CSS")
  expect_null(cl$branching_direction)
})

test_that("2-trait DD99 branches in both directions when both kernels are narrow", {
  cl <- community_solve_singularity(
    dd99_2d(sigma_C = c(0.4, 0.5), sigma_K = c(1, 1))) |>
    community_classify_singularity()
  expect_true(all(cl$hessian_eigen$values > 0))
  expect_false(cl$evolutionarily_stable)
  expect_equal(cl$classification, "branching point")
  # steepest disruptive direction is the narrower kernel (trait 1)
  expect_equal(unname(abs(cl$branching_direction)), c(1, 0), tolerance = 1e-6)
})

test_that("community_classify_singularity needs exactly one resident", {
  comm <- dd99_1d() |>
    community_add(trait_matrix(c(-1, 1), "x"), birth_rate = c(100, 100)) |>
    community_demography()
  expect_error(community_classify_singularity(comm),
               "needs exactly one resident")
})

test_that("community_classify_singularity solves demography if needed", {
  comm <- dd99_1d() |> community_add(trait_matrix(0, "x"), birth_rate = 100)
  expect_null(comm$fitness_function)
  cl <- community_classify_singularity(comm)
  expect_equal(cl$classification, "branching point")
})

test_that("the classification prints its verdict and eigenvalues", {
  cl <- community_solve_singularity(dd99_1d()) |>
    community_classify_singularity()
  out <- paste(utils::capture.output(print(cl)), collapse = "\n")
  expect_match(out, "branching point")
  expect_match(out, "Hessian eigenvalues")
  expect_match(out, "branching direction")
})
