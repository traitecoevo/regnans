# Tests for the Geritz et al. 1998 soft-selection toy model (issue #33).
# Analytic oracles (symmetric three-patch, mu = (-d, 0, d), equal capacities):
# singular strategy x* = sum c_j mu_j = 0; convergence stable for all d/sigma;
# branching point iff d/sigma > sqrt(3/2) ~= 1.2247.

curvature <- function(comm, at = 0, d = 1e-4) {
  f <- comm$fitness_function
  (f(at + d) - 2 * f(at) + f(at - d)) / (d * d)
}

test_that("gk98 resident fitness is zero at the singular strategy", {
  pars <- list(sigma = 1, mu = c(-1.5, 0, 1.5), K = c(1, 1, 1))
  n <- gk98_equilibrium(0, pars)
  expect_equal(gk98_fitness(0, 0, n, pars), 0, tolerance = 1e-12)
})

test_that("gk98 singular strategy is the capacity-weighted mean optimum", {
  # symmetric -> x* = 0
  h <- harness_gk98(d = 1.5, sigma = 1)
  out <- community_start(bounds(x = c(-4, 4)), harness = h, trait_scale = "linear") |>
    community_solve_singularity()
  expect_equal(as.numeric(out$traits), 0, tolerance = 1e-4)

  # asymmetric capacities -> x* = sum(c_j mu_j)
  mu <- c(-1, 0, 2); K <- c(1, 2, 1)
  xstar <- sum(K / sum(K) * mu)
  h2 <- harness_gk98(sigma = 1, mu = mu, K = K)
  out2 <- community_start(bounds(x = c(-4, 4)), harness = h2, trait_scale = "linear") |>
    community_solve_singularity()
  expect_equal(as.numeric(out2$traits), xstar, tolerance = 1e-3)
})

test_that("gk98 branches iff d/sigma > sqrt(3/2)", {
  css <- community_start(bounds(x = c(-4, 4)),
                         harness = harness_gk98(d = 1.0, sigma = 1)) |>
    community_add(trait_matrix(0, "x")) |>
    community_demography()
  expect_lt(curvature(css), 0) # CSS (fitness maximum)

  branch <- community_start(bounds(x = c(-4, 4)),
                            harness = harness_gk98(d = 1.5, sigma = 1)) |>
    community_add(trait_matrix(0, "x")) |>
    community_demography()
  expect_gt(curvature(branch), 0) # branching (fitness minimum)
})

test_that("gk98 branching threshold matches sqrt(3/2) analytically", {
  curv0 <- function(d) {
    pars <- list(sigma = 1, mu = c(-d, 0, d), K = c(1, 1, 1))
    h <- 1e-4
    n <- gk98_equilibrium(0, pars)
    (gk98_fitness(h, 0, n, pars) - 2 * gk98_fitness(0, 0, n, pars) +
       gk98_fitness(-h, 0, n, pars)) / (h * h)
  }
  expect_lt(curv0(1.2), 0)
  expect_gt(curv0(1.25), 0)
  # curvature crosses zero at d = sqrt(3/2)
  root <- uniroot(curv0, c(1.0, 1.5))$root
  expect_equal(root, sqrt(3 / 2), tolerance = 1e-3)
})
