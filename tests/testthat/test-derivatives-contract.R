# The contract every harness that advertises a derivative signs: it agrees
# with a finite difference of its own fitness function. Run over every shipped
# harness that provides anything; a harness author elsewhere calls
# harness_check_derivatives() from their own tests the same way.

# One entry per shipped harness: the harness, the bounds to test inside, the
# trait scale, and the residents to solve before checking. GK98 and GM99 solve
# multi-resident equilibria numerically, so each also gets a pair.
shipped <- list(
  dd99_narrow = list(harness_dd99(sigma_C = 0.4, sigma_K = 1),
                     bounds(x = c(-2, 2)), "linear", c(0.4, -0.9)),
  dd99_wide   = list(harness_dd99(sigma_C = 1.5, sigma_K = 1, x0 = 0.6, r = 2, K0 = 50),
                     bounds(x = c(-2, 2)), "linear", 0.4),
  dd99_nd     = list(harness_dd99_nd(x0 = c(0.3, -0.5), sigma_C = c(0.4, 1.5)),
                     bounds(x1 = c(-2, 2), x2 = c(-2, 2)), "linear",
                     matrix(c(0.4, -0.2), 1)),
  gk98_css    = list(harness_gk98(d = 1.0, sigma = 1), bounds(x = c(-3, 3)), "linear", 0.3),
  gk98_branch = list(harness_gk98(d = 1.5, sigma = 1), bounds(x = c(-3, 3)), "linear",
                     c(-0.8, 0.8)),
  gk98_uneven = list(harness_gk98(mu = c(-1, 0.5, 2), K = c(1, 2, 0.5)),
                     bounds(x = c(-3, 3)), "linear", 0.3),
  jj12        = list(harness_jj12(a = 0.1, x_opt = 0.5, sigma = 1.4),
                     bounds(x = c(-3, 3)), "linear", 0.2),
  gm99_css    = list(harness_gm99(alpha = 4.5, beta = 15, R = 1),
                     bounds(x = c(0.08, 0.9)), "log", 0.3),
  gm99_branch = list(harness_gm99(alpha = 7, beta = 15, R = 1),
                     bounds(x = c(0.08, 0.9)), "log", c(0.2, 0.5))
)

community_for <- function(spec) {
  h <- spec[[1]]
  comm <- community_start(spec[[2]], trait_scale = spec[[3]], harness = h)
  residents <- spec[[4]]
  if (!is.matrix(residents)) residents <- trait_matrix(residents, h$trait_names)
  colnames(residents) <- h$trait_names
  comm |>
    community_add(residents, birth_rate = rep(1, nrow(residents))) |>
    community_demography()
}

test_that("every shipped harness honours the derivative contract", {
  for (nm in names(shipped)) {
    res <- harness_check_derivatives(community_for(shipped[[nm]]), n_points = 6)
    expect_true(all(res$pass), info = nm)
    expect_setequal(unique(res$quantity), c("fitness_gradient", "fitness_hessian"))
    # the pass column is the contract; beyond it, a first difference of a
    # smooth fitness is good to ~1e-9, so the model gradient must match that
    # closely wherever the gradient is not itself near zero
    g <- res[res$quantity == "fitness_gradient" & abs(res$fd) > 1e-3, ]
    expect_lt(max(g$rel_err), 1e-6)
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
  fd_only <- harness_explicit(
    fitness = dd99_fitness,
    equilibrium = function(x_res, pars) dd99_equilibrium(x_res, pars),
    pars = list(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4),
    trait_names = "x", label = "dd99_fd")
  comm <- community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                          harness = fd_only)
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
  comm <- community_for(list(h, bounds(x = c(-2, 2)), "linear", 0.4))
  err <- tryCatch(harness_check_derivatives(comm), error = identity)
  expect_s3_class(err, "derivative_check_error")
  expect_match(conditionMessage(err), "disagree with finite differences")
  expect_match(conditionMessage(err), "fitness_gradient\\[x\\]")
  expect_true(is.data.frame(err$results))
  expect_false(all(err$results$pass))
  expect_equal(unique(err$results$quantity), "fitness_gradient")
})

test_that("the check refuses a community whose closures do not match its harness", {
  comm <- community_for(shipped$dd99_wide)
  comm$fitness_derivatives$source[] <- "finite difference"
  expect_error(harness_check_derivatives(comm), "no model-supplied closure")
})
