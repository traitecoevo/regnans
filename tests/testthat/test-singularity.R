# Singular-strategy solving (any dimension; Newton or a 1-D bracket) and
# classification.
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

test_that("the bracket solver recovers the 1D DD99 singular strategy", {
  out <- community_solve_singularity(dd99_1d(x0 = -0.35), solver = "bracket",
                                     tol = 1e-8)
  expect_true(attr(out, "converged"))
  expect_equal(attr(out, "solver"), "bracket")
  expect_equal(unname(attr(out, "singularity")), -0.35, tolerance = 1e-6)
  expect_equal(as.numeric(out$traits), -0.35, tolerance = 1e-6)
  expect_lt(abs(out$selection_gradient), 1e-6)
  # Newton from the midpoint finds the same root
  nd <- community_solve_singularity(dd99_1d(x0 = -0.35))
  expect_equal(as.numeric(nd$traits), as.numeric(out$traits), tolerance = 1e-5)
})

test_that("the bracket solver searches on the trait scale", {
  # GM99 seed size is strictly positive; the root is the same on either scale
  h <- harness_gm99(alpha = 4.5, beta = 15)
  log_root <- community_start(bounds(x = c(0.06, 0.9)), harness = h,
                              trait_scale = "log") |>
    community_solve_singularity(solver = "bracket", tol = 1e-8)
  lin_root <- community_start(bounds(x = c(0.06, 0.9)), harness = h,
                              trait_scale = "linear") |>
    community_solve_singularity(solver = "bracket", tol = 1e-8)
  expect_equal(as.numeric(log_root$traits), as.numeric(lin_root$traits),
               tolerance = 1e-6)
})

test_that("the bracket solver returns the bound selection pushes towards", {
  # x* = 0 lies below [0.5, 1.5]: the gradient is negative at both ends
  expect_warning(
    edge <- community_solve_singularity(dd99_1d(), bounds = c(0.5, 1.5),
                                        solver = "bracket"),
    "Bounds do not include a singularity")
  expect_equal(as.numeric(edge$traits), 0.5)
  expect_false(attr(edge, "converged"))
  # and above it, the upper bound
  expect_warning(
    edge <- community_solve_singularity(dd99_1d(), bounds = c(-1.5, -0.5),
                                        solver = "bracket"),
    "Bounds do not include a singularity")
  expect_equal(as.numeric(edge$traits), -0.5)

  expect_error(
    community_solve_singularity(dd99_1d(), bounds = c(0.5, 1.5),
                                solver = "bracket", edge_ok = FALSE),
    "Bounds do not include a singularity")
})

test_that("the bracket solver needs a single trait", {
  expect_error(community_solve_singularity(dd99_2d(), solver = "bracket"),
               "needs a single trait")
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

test_that("community_classify_singularity needs a resident", {
  expect_error(community_classify_singularity(dd99_1d()), "at least one resident")
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

# ---- singular coalitions ----------------------------------------------------
#
# Oracles in helper-coalition.R: the DD99 mirror-image pair in closed form, the
# DD99 symmetric triple and the GK98 pair as scalar roots of their closed-form
# gradients, and closed-form curvatures at each resident.

gk98_1d <- function(d = 1.5) {
  community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                  harness = harness_gk98(d = d, sigma = 1))
}

test_that("community_solve_singularity finds the DD99 dimorphic coalition in closed form", {
  a <- dd99_pair_root(sigma_C = 0.4, sigma_K = 1)
  for (solver in c("nleqslv", "newton")) {
    out <- community_solve_singularity(dd99_1d(), x0 = c(-0.3, 0.6), solver = solver, tol = 1e-10)
    expect_true(attr(out, "converged"), info = solver)
    expect_equal(as.numeric(out$traits), c(-a, a), tolerance = 1e-8, info = solver)
    expect_equal(out$birth_rate[1], out$birth_rate[2], tolerance = 1e-8, info = solver)
    expect_lt(max(abs(out$selection_gradient)), 1e-9)
    root <- attr(out, "singularity")
    expect_equal(dim(root), c(2L, 1L))
    expect_equal(colnames(root), "x")
    # a Newton solve: one Jacobian across two residents and a few steps
    expect_lte(attr(out, "evaluations"), 12L)
  }
  # derivative-free, and a looser residual tolerance
  df <- community_solve_singularity(dd99_1d(), x0 = c(-0.3, 0.6), solver = "dfsane")
  expect_equal(as.numeric(df$traits), c(-a, a), tolerance = 1e-5)
  # by default the community's own residents are the starting coalition, and
  # their densities the starting equilibrium
  start <- dd99_1d() |>
    community_add(trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(100, 100)) |>
    community_demography()
  expect_equal(as.numeric(community_solve_singularity(start)$traits), c(-a, a), tolerance = 1e-6)
})

test_that("community_solve_singularity finds the DD99 trimorphic coalition", {
  b <- dd99_triple_root(sigma_C = 0.3)
  out <- community_solve_singularity(dd99_1d(sigma_C = 0.3), x0 = c(-0.5, 0.05, 0.8), tol = 1e-10)
  expect_true(attr(out, "converged"))
  expect_equal(as.numeric(out$traits), c(-b, 0, b), tolerance = 1e-8)
  expect_equal(out$birth_rate, dd99_coalition_density(c(-b, 0, b), sigma_C = 0.3), tolerance = 1e-6)
})

test_that("community_solve_singularity finds the GK98 dimorphic coalition", {
  a <- gk98_pair_root(d = 1.5)
  out <- community_solve_singularity(gk98_1d(), x0 = c(-0.5, 1), tol = 1e-10)
  expect_true(attr(out, "converged"))
  expect_equal(as.numeric(out$traits), c(-a, a), tolerance = 1e-8)
})

test_that("community_solve_singularity searches a coalition on a log trait scale", {
  # GM99 at alpha R = 7 branches at m* ~ 0.645; its daughters settle at a
  # dimorphism no closed form gives, so the oracle is internal: every
  # resident's gradient vanishes, every resident has invasion fitness zero, and
  # a Newton solve with a correct trait-scale Jacobian gets there in few solves
  gm <- community_start(bounds(x = c(0.1, 0.95)), harness = harness_gm99(alpha = 7, beta = 15))
  out <- community_solve_singularity(gm, x0 = c(0.5, 0.8), tol = 1e-8)
  expect_true(attr(out, "converged"))
  x <- as.numeric(out$traits)
  expect_lt(max(abs(out$selection_gradient * x)), 1e-7)
  expect_lt(max(abs(out$resident_fitness)), 1e-8)
  expect_true(all(out$birth_rate > 0))
  expect_lte(attr(out, "evaluations"), 20L)
})

test_that("a coalition that cannot coexist is reported lost, not returned as converged", {
  # sigma_C > sigma_K: no dimorphism exists, and the pair merges at x0 (or one
  # resident is driven out on the way, depending on the solver's path)
  for (solver in c("nleqslv", "newton", "dfsane")) {
    expect_warning(
      out <- community_solve_singularity(dd99_1d(sigma_C = 1.5), x0 = c(-0.3, 0.6), solver = solver),
      "lost the coalition", info = solver)
    expect_false(attr(out, "converged"), info = solver)
    expect_equal(nrow(out$traits), 2L)
  }
  expect_error(community_solve_singularity(dd99_1d(), x0 = c(0.2, 0.2)), "distinct residents")
  expect_error(community_solve_singularity(dd99_1d(), x0 = c(-0.3, 0.6), solver = "bracket"),
               "needs a single resident")
  expect_error(community_solve_singularity(dd99_2d(), x0 = c(0.1, 0.2, 0.3)), "one value per trait")
})

test_that("the DD99 pair classifies as a convergence-stable branching coalition", {
  a <- dd99_pair_root()
  sol <- community_solve_singularity(dd99_1d(), x0 = c(-0.3, 0.6), tol = 1e-10)
  cl <- community_classify_singularity(sol)
  expect_s3_class(cl, "singularity_classification")
  expect_equal(cl$classification, "branching point")
  expect_false(cl$evolutionarily_stable)
  expect_equal(cl$resident_evolutionarily_stable, c(FALSE, FALSE))
  # each resident sits at a fitness minimum of the closed-form curvature
  curv <- dd99_coalition_curvature(c(-a, a))
  expect_true(all(curv > 0))
  expect_equal(vapply(cl$hessian, as.numeric, numeric(1)), curv, tolerance = 1e-5)
  # the coalition Jacobian against the closed form: both eigenvalues negative
  J <- oracle_jacobian(dd99_coalition_gradient, c(-a, a))
  expect_equal(sort(Re(cl$jacobian_eigen$values)), sort(eigen(J)$values), tolerance = 1e-5)
  expect_true(cl$convergence_stable)
  expect_true(cl$strongly_convergence_stable)
  expect_equal(dim(cl$branching_direction), c(2L, 1L))
  expect_equal(abs(as.numeric(cl$branching_direction)), c(1, 1))
  expect_equal(dim(cl$traits), c(2L, 1L))
  expect_equal(cl$evaluations, 1L + 2L * 2L)
  out <- paste(utils::capture.output(print(cl)), collapse = "\n")
  expect_match(out, "coalition of 2")
  expect_match(out, "branching direction of resident 2")
})

test_that("the GK98 pair classifies as an evolutionarily stable coalition", {
  a <- gk98_pair_root(d = 1.5)
  cl <- community_solve_singularity(gk98_1d(), x0 = c(-0.5, 1), tol = 1e-10) |>
    community_classify_singularity()
  expect_equal(cl$classification, "CSS")
  expect_equal(cl$resident_evolutionarily_stable, c(TRUE, TRUE))
  expect_equal(vapply(cl$hessian, as.numeric, numeric(1)),
               rep(gk98_pair_curvature(a), 2), tolerance = 1e-5)
  expect_null(cl$branching_direction)
})

test_that("the DD99 trimorphic coalition classifies as a branching coalition", {
  b <- dd99_triple_root(sigma_C = 0.3)
  sol <- community_solve_singularity(dd99_1d(sigma_C = 0.3), x0 = c(-0.5, 0.05, 0.8), tol = 1e-10)
  cl <- community_classify_singularity(sol)
  expect_equal(cl$classification, "branching point")
  expect_equal(vapply(cl$hessian, as.numeric, numeric(1)),
               dd99_coalition_curvature(c(-b, 0, b), sigma_C = 0.3), tolerance = 1e-5)
  J <- oracle_jacobian(function(x) dd99_coalition_gradient(x, sigma_C = 0.3), c(-b, 0, b))
  expect_equal(sort(Re(cl$jacobian_eigen$values)), sort(Re(eigen(J)$values)), tolerance = 1e-5)
  expect_true(cl$convergence_stable)

  # the middle resident is denser than the outer two, so it evolves faster:
  # convergence stability is decided on diag(speeds) J, speeds the densities
  n <- dd99_coalition_density(c(-b, 0, b), sigma_C = 0.3)
  w <- n / mean(n)
  expect_equal(cl$speeds, w, tolerance = 1e-6)
  expect_gt(max(abs(w - 1)), 0.1)
  expect_equal(sort(Re(cl$jacobian_weighted_eigen$values)), sort(Re(eigen(w * J)$values)),
               tolerance = 1e-5)
  expect_match(paste(utils::capture.output(print(cl)), collapse = "\n"), "speed-weighted")
  # equal speeds, or any given ones
  eq <- community_classify_singularity(sol, speeds = "equal")
  expect_equal(eq$speeds, rep(1, 3))
  expect_equal(eq$jacobian_weighted_eigen$values, eq$jacobian_eigen$values)
  given <- community_classify_singularity(sol, speeds = c(1, 4, 1))
  expect_equal(given$speeds, c(0.5, 2, 0.5))
  expect_equal(sort(Re(given$jacobian_weighted_eigen$values)),
               sort(Re(eigen(c(0.5, 2, 0.5) * given$jacobian)$values)), tolerance = 1e-10)
  expect_error(community_classify_singularity(sol, speeds = c(1, 1)), "3 positive numbers")
  expect_error(community_classify_singularity(sol, speeds = c(1, -1, 1)), "3 positive numbers")
})

test_that("a two-trait coalition is solved and classified with residents stacked trait by trait", {
  # sigma_C = (0.4, 1.5): the DD99 pair splits along x1 only, at the 1-D
  # closed-form root, with both residents at x2 = 0; across x2 each sits at a
  # fitness maximum of curvature r (1/sigma_C2^2 - 1/sigma_K2^2)
  a <- dd99_pair_root(sigma_C = 0.4)
  out <- community_solve_singularity(dd99_2d(), x0 = rbind(c(-0.3, 0.1), c(0.6, -0.1)), tol = 1e-10)
  expect_true(attr(out, "converged"))
  expect_equal(unname(out$traits), cbind(c(-a, a), c(0, 0)), tolerance = 1e-8)
  root <- attr(out, "singularity")
  expect_equal(dim(root), c(2L, 2L))
  expect_equal(colnames(root), c("x1", "x2"))

  cl <- community_classify_singularity(out)
  expect_equal(cl$classification, "branching point")
  expect_equal(dim(cl$traits), c(2L, 2L))
  expect_equal(rownames(cl$jacobian), c("x1[1]", "x1[2]", "x2[1]", "x2[2]"))
  # the x1 block is the 1-D coalition Jacobian; x1 and x2 decouple at x2 = 0
  J1 <- oracle_jacobian(dd99_coalition_gradient, c(-a, a))
  expect_equal(unname(cl$jacobian[1:2, 1:2]), J1, tolerance = 1e-5)
  expect_equal(unname(cl$jacobian[1:2, 3:4]), matrix(0, 2, 2), tolerance = 1e-6)
  expect_equal(unname(cl$jacobian[3:4, 1:2]), matrix(0, 2, 2), tolerance = 1e-6)
  curv <- dd99_coalition_curvature(c(-a, a))
  for (i in 1:2) {
    expect_equal(unname(diag(cl$hessian[[i]])), c(curv[i], 1 / 1.5^2 - 1), tolerance = 1e-5)
  }
  expect_equal(unname(abs(cl$branching_direction)), cbind(c(1, 1), c(0, 0)), tolerance = 1e-6)
})

test_that("the resident Jacobian refuses residents too close to tell apart", {
  # 8e-4 apart at x ~ 0.5, inside two finite-difference steps (1e-3 |x| each)
  close <- dd99_1d() |>
    community_add(trait_matrix(c(0.5, 0.5008), "x"), birth_rate = c(100, 100)) |>
    community_demography()
  expect_error(community_selection_gradient_jacobian(close),
               "residents 1 and 2 are within the finite-difference step")
  expect_error(community_classify_singularity(close), "within the finite-difference step")
})
