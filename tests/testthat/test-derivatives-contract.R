# The contract every harness that advertises a derivative signs: it agrees
# with a finite difference of its own fitness function. Run over every shipped
# harness that provides anything; a harness author elsewhere calls
# harness_check_derivatives() from their own tests the same way.

dd99_harnesses <- list(
  narrow   = harness_dd99(sigma_C = 0.4, sigma_K = 1),
  wide     = harness_dd99(sigma_C = 1.5, sigma_K = 1, x0 = 0.6, r = 2, K0 = 50),
  nd       = harness_dd99_nd(x0 = c(0.3, -0.5), sigma_C = c(0.4, 1.5))
)

community_for <- function(h, trait_scale = "linear") {
  k <- length(h$trait_names)
  b <- if (k == 1L) bounds(x = c(-2, 2)) else bounds(x1 = c(-2, 2), x2 = c(-2, 2))
  comm <- community_start(b, trait_scale = trait_scale, harness = h)
  resident <- if (k == 1L) 0.4 else c(0.4, -0.2)
  comm |>
    community_add(trait_matrix(resident, h$trait_names), birth_rate = 100) |>
    community_demography()
}

test_that("every shipped DD99 harness honours the derivative contract", {
  for (nm in names(dd99_harnesses)) {
    res <- harness_check_derivatives(community_for(dd99_harnesses[[nm]]),
                                     n_points = 6)
    expect_true(all(res$pass), info = nm)
    expect_setequal(unique(res$quantity), c("fitness_gradient", "fitness_hessian"))
    expect_lt(max(res$rel_err), 1e-6)
  }
})

test_that("the contract holds on a log trait scale and for a multi-resident community", {
  h <- harness_dd99(x0 = 0.5, sigma_K = 1, sigma_C = 0.4, trait_name = "lma")
  comm <- community_start(bounds(lma = c(0.05, 2)), trait_scale = "log",
                          harness = h) |>
    community_add(trait_matrix(c(0.2, 0.9), "lma"), birth_rate = c(100, 100)) |>
    community_demography()
  res <- harness_check_derivatives(comm, n_points = 4)
  expect_true(all(res$pass))
  expect_equal(nrow(res), 4 * (1 + 1))
})

test_that("an empty community is solved before checking", {
  comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                          harness = harness_dd99())
  res <- harness_check_derivatives(comm, n_points = 3)
  expect_true(all(res$pass))
})

test_that("a harness without model derivatives has nothing to check", {
  comm <- community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                          harness = harness_gk98())
  res <- harness_check_derivatives(comm)
  expect_equal(nrow(res), 0L)
  expect_named(res, c("quantity", "point", "entry", "model", "fd",
                      "abs_err", "rel_err", "pass"))
})

test_that("a wrong derivative fails the contract with the evidence attached", {
  pars <- list(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4)
  h <- harness_explicit(
    fitness = dd99_fitness,
    equilibrium = function(x_res, pars) dd99_equilibrium(x_res, pars),
    pars = pars, trait_names = "x", label = "skewed",
    fitness_gradient = function(x_mut, x_res, n_res, pars)
      1.1 * dd99_fitness_gradient(x_mut, x_res, n_res, pars))
  comm <- community_for(h)
  err <- tryCatch(harness_check_derivatives(comm), error = identity)
  expect_s3_class(err, "derivative_check_error")
  expect_match(conditionMessage(err), "disagree with finite differences")
  expect_match(conditionMessage(err), "fitness_gradient\\[x\\]")
  expect_true(is.data.frame(err$results))
  expect_false(all(err$results$pass))
  expect_equal(unique(err$results$quantity), "fitness_gradient")
})

test_that("the check refuses a community whose closures do not match its harness", {
  comm <- community_for(harness_dd99())
  comm$fitness_derivatives$source[] <- "finite difference"
  expect_error(harness_check_derivatives(comm), "no model-supplied closure")
})
