# Tests for the toy-model harness abstraction and the JJ12 bird arrival-time model
# (issue #33). These run no SCM, so they are fast and validate the assembly /
# attractor machinery against the model's analytic answers.

# ---- C++ model oracles (direct calls) -------------------------------------

test_that("jj12 model invasion fitness is zero for resident at equilibrium", {
  pars <- list(a = 0.125, x_opt = 0, sigma = 10, R0 = 1, K = 5, p = 0.5)
  # w(x, x) == 1 for ANY x and any p in (0,1)  =>  log fitness == 0
  for (x in c(-20, -5, 0, 3.7, 15)) {
    n <- jj12_equilibrium(x, pars)
    expect_equal(jj12_fitness(x, x, n, pars), 0, tolerance = 1e-12)
  }
})

test_that("jj12 single-resident equilibrium matches K R(x)/(1-p)", {
  pars <- list(a = 0.1, x_opt = 2, sigma = 5, R0 = 3, K = 7, p = 0.4)
  for (x in c(-3, 2, 6)) {
    Rx <- pars$R0 * exp(-(x - pars$x_opt)^2 / (2 * pars$sigma^2))
    expect_equal(jj12_equilibrium(x, pars),
                 pars$K * Rx / (1 - pars$p), tolerance = 1e-12)
  }
})

test_that("jj12 selection gradient vanishes at x* = x_opt - a sigma^2", {
  pars <- list(a = 0.125, x_opt = 0, sigma = 10, R0 = 1, K = 5, p = 0.5)
  xstar <- pars$x_opt - pars$a * pars$sigma^2
  d <- 1e-5
  grad <- function(r) {
    n <- jj12_equilibrium(r, pars)
    (jj12_fitness(r + d, r, n, pars) -
       jj12_fitness(r - d, r, n, pars)) / (2 * d)
  }
  expect_equal(grad(xstar), 0, tolerance = 1e-6)
  # convergence stable: gradient points back toward x* on both sides
  expect_gt(grad(xstar - 10), 0)
  expect_lt(grad(xstar + 10), 0)
})

test_that("jj12 fitness is vectorised over mutants", {
  pars <- list(a = 0.1, x_opt = 0, sigma = 10, R0 = 1, K = 5, p = 0.5)
  x <- 0
  n <- jj12_equilibrium(x, pars)
  muts <- c(-5, 0, 5)
  w <- jj12_fitness(muts, x, n, pars)
  expect_length(w, 3L)
  expect_equal(w[2], 0, tolerance = 1e-12) # mutant == resident
})

# ---- harness object --------------------------------------------------------

test_that("harness_jj12 builds an explicit harness with the connectors and its own equilibrium", {
  h <- harness_jj12(a = 0.1, x_opt = 0, sigma = 10)
  expect_s3_class(h, "harness_jj12")
  expect_s3_class(h, "harness_explicit")
  expect_s3_class(h, "harness")
  expect_setequal(
    names(h$fns),
    c("parameters", "make_demography_runner", "demography_runner_cleanup",
      "viable_bounds", "check_for_inviable_strategies",
      "update_fitness_function", "equilibrium"))
})

test_that("community_start defaults to the plant harness", {
  comm <- community_start(bounds(lma = c(0.01, 2)))
  expect_s3_class(comm$harness, "harness_plant")
})

# ---- full pipeline against the analytic singular strategy ------------------

test_that("community_demography solves the jj12 resident to equilibrium", {
  h <- harness_jj12(a = 0.125, x_opt = 0, sigma = 10, R0 = 1, K = 5, p = 0.5)
  comm <- community_start(bounds(x = c(-40, 40)), harness = h, trait_scale = "linear") |>
    community_add(trait_matrix(0, "x")) |>
    community_demography()
  expect_true(attr(comm, "converged"))
  expect_equal(comm$birth_rate, 10, tolerance = 1e-8)       # K R(0)/(1-p)
  expect_equal(comm$resident_fitness, 0, tolerance = 1e-10) # at equilibrium
})

test_that("community_solve_singularity recovers the jj12 CSS", {
  for (a in c(0, 0.04, 0.125)) {
    h <- harness_jj12(a = a, x_opt = 0, sigma = 10, R0 = 1, K = 5, p = 0.5)
    out <- community_start(bounds(x = c(-40, 40)), harness = h, trait_scale = "linear") |>
      community_solve_singularity()
    expect_equal(as.numeric(out$traits), -a * 100, tolerance = 1e-3)
  }
})

test_that("jj12 singular strategy shifts with x_opt and sigma", {
  h <- harness_jj12(a = 0.05, x_opt = 8, sigma = 6, R0 = 1, K = 5, p = 0.5)
  out <- community_start(bounds(x = c(-40, 40)), harness = h, trait_scale = "linear") |>
    community_solve_singularity()
  expect_equal(as.numeric(out$traits), 8 - 0.05 * 36, tolerance = 1e-3)
})
