# Tests for the Dieckmann & Doebeli 1999 toy model (issue #33). Analytic
# oracles: singular strategy x* = x0; branching point iff sigma_C < sigma_K.

curvature <- function(comm, at = 0, d = 1e-4) {
  f <- comm$fitness_function
  (f(at + d) - 2 * f(at) + f(at - d)) / (d * d)
}

test_that("dd99 resident fitness is zero and equilibrium is K(x)", {
  pars <- list(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4)
  for (x in c(-1.5, 0, 0.7)) {
    Kx <- pars$K0 * exp(-(x - pars$x0)^2 / (2 * pars$sigma_K^2))
    expect_equal(dd99_equilibrium(x, pars), Kx, tolerance = 1e-10)
    expect_equal(dd99_fitness(x, x, dd99_equilibrium(x, pars), pars), 0,
                 tolerance = 1e-12)
  }
})

test_that("dd99 multi-resident equilibrium solves the competition system", {
  pars <- list(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4)
  x <- c(-1, 1)
  N <- dd99_equilibrium(x, pars)
  # at equilibrium every resident's invasion fitness is ~0
  expect_equal(dd99_fitness(x, x, N, pars), c(0, 0), tolerance = 1e-8)
})

test_that("dd99 singular strategy is x0", {
  h <- harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4)
  out <- community_start(bounds(x = c(-2, 2)), harness = h, trait_scale = "linear") |>
    community_solve_singularity()
  expect_equal(as.numeric(out$traits), 0, tolerance = 1e-4)

  h2 <- harness_dd99(x0 = 0.6, sigma_K = 1, sigma_C = 0.4)
  out2 <- community_start(bounds(x = c(-2, 2)), harness = h2, trait_scale = "linear") |>
    community_solve_singularity()
  expect_equal(as.numeric(out2$traits), 0.6, tolerance = 1e-4)
})

test_that("dd99 branches iff sigma_C < sigma_K", {
  branch <- community_start(bounds(x = c(-2, 2)),
                            harness = harness_dd99(sigma_K = 1, sigma_C = 0.4)) |>
    community_add(trait_matrix(0, "x")) |>
    community_demography()
  expect_gt(curvature(branch), 0) # fitness minimum -> branching

  ess <- community_start(bounds(x = c(-2, 2)),
                         harness = harness_dd99(sigma_K = 1, sigma_C = 1.5)) |>
    community_add(trait_matrix(0, "x")) |>
    community_demography()
  expect_lt(curvature(ess), 0)    # fitness maximum -> ESS
})
