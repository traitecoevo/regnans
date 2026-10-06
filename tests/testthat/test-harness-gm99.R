# Tests for the Geritz, van der Meijden & Metz 1999 seed-size safe-site model
# (issue #33; the model in Daniel's MATLAB). There is no closed-form singular
# strategy, so the strong analytic oracle is the demographic-equilibrium
# invariant W_m(m) = 1 (log fitness 0); singular-strategy locations are
# regression values cross-checked against the paper's qualitative Fig. 5
# behaviour (alpha*R = 4.5 -> CSS, alpha*R = 7 -> branching at beta*R = 15).

curvature <- function(comm, at, d = 1e-5) {
  f <- comm$fitness_function
  (f(at + d) - 2 * f(at) + f(at - d)) / (d * d)
}

test_that("gm99 resident at equilibrium has log fitness zero", {
  pars <- list(R = 1, alpha = 7, beta = 15)
  for (m in c(0.1, 0.3, 0.645, 0.9)) {
    N <- gm99_equilibrium(m, pars)
    expect_gt(N, 0)
    expect_equal(gm99_fitness(m, m, N, pars), 0, tolerance = 1e-10)
    # equilibrium invariant: (R/m) s(m) (1 - e^-N)/N == 1
    s <- max(0, 1 - 2 * exp(-pars$beta * m))
    expect_equal((pars$R / m) * s * (1 - exp(-N)) / N, 1, tolerance = 1e-8)
  }
})

test_that("gm99 self-competition matches the closed form (1-e^-N)/N", {
  pars <- list(R = 1, alpha = 4.5, beta = 15)
  # W(m,m) = (R/m) s(m) g_res(N); so g_res(N) = W / ((R/m) s(m)).
  m <- 0.3; N <- 2.5
  s <- 1 - 2 * exp(-pars$beta * m)
  logW <- gm99_fitness(m, m, N, pars)
  g_res <- exp(logW) / ((pars$R / m) * s)
  expect_equal(g_res, (1 - exp(-N)) / N, tolerance = 1e-8)
})

test_that("gm99 non-viable seed sizes have zero equilibrium density", {
  pars <- list(R = 1, alpha = 6, beta = 15)
  # just above the minimum viable size m_min = ln2/beta, (R/m)s(m) < 1
  m_min <- log(2) / pars$beta
  expect_equal(gm99_equilibrium(m_min * 1.005, pars), 0, tolerance = 1e-12)
  # log fitness is -Inf below m_min (s = 0)
  expect_equal(gm99_fitness(m_min * 0.5, 0.3,
                                    gm99_equilibrium(0.3, pars), pars),
               -Inf)
})

test_that("gm99 is a CSS at low asymmetry (alpha*R = 4.5, beta*R = 15)", {
  h <- harness_gm99(alpha = 4.5, beta = 15)
  out <- community_start(bounds(x = c(0.06, 0.9)), harness = h) |>
    community_solve_singularity()
  mstar <- as.numeric(out$traits)
  expect_equal(mstar, 0.183, tolerance = 5e-3)   # numerical singular strategy
  expect_lt(curvature(out, mstar), 0)            # fitness maximum -> ESS/CSS
})

test_that("gm99 branches at high asymmetry (alpha*R = 7, beta*R = 15)", {
  h <- harness_gm99(alpha = 7, beta = 15)
  out <- community_start(bounds(x = c(0.1, 0.95)), harness = h) |>
    community_solve_singularity()
  mstar <- as.numeric(out$traits)
  expect_equal(mstar, 0.645, tolerance = 5e-3)   # numerical singular strategy
  expect_gt(curvature(out, mstar), 0)            # fitness minimum -> branching
})

test_that("gm99 depends only on the products alpha*R and beta*R", {
  # same alpha*R, beta*R via different R should give the same singular strategy
  h1 <- harness_gm99(alpha = 7, beta = 15, R = 1)
  h2 <- harness_gm99(alpha = 3.5, beta = 7.5, R = 2)
  s1 <- community_start(bounds(x = c(0.1, 0.95)), harness = h1) |>
    community_solve_singularity()
  # second uses seed sizes in units up to R = 2, so x* scales by R: compare x*/R
  s2 <- community_start(bounds(x = c(0.2, 1.9)), harness = h2) |>
    community_solve_singularity()
  expect_equal(as.numeric(s1$traits), as.numeric(s2$traits) / 2,
               tolerance = 5e-3)
})
