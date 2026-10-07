# The parallel map and the model evaluations routed through it.

# DD99 solved by the Newton equilibrium solver rather than the model's own
# equilibrium, so that where each solve starts from matters to its answer.
dd99_newton <- function(r_jacobian = 1L) {
  community_start(bounds(x = c(-2, 2)), trait_scale = "linear",
                  demography_control = demographic_step_control(list(equilibrium_solver_name = "equilibrium_solve_newton")),
                  derivative_control = derivative_control(list(r_jacobian = r_jacobian)),
                  harness = harness_dd99(sigma_C = 0.4))
}
dd99_newton_off <- function(r_jacobian = 1L) {
  dd99_newton(r_jacobian) |>
    community_add(trait_matrix(c(-0.3, 0.6), "x"), birth_rate = c(300, 300)) |>
    community_demography()
}

# ---- the sequential path is pinned -------------------------------------------
# With no plan, every routed path must be the arithmetic it was before routing,
# the warm-start chain from one stencil point to the next included. Seeding
# each point from the centre instead moves these values far beyond 1e-12 but
# within the oracle tolerances of the other tests, so they are pinned here.

test_that("the resident Jacobian, singularity solve and classification are pinned sequentially", {
  J <- community_selection_gradient_jacobian(dd99_newton_off(r_jacobian = 2L))
  expect_equal(as.numeric(J), c(-2.3204965543239702, 2.1524899760059748, 1.67464366695251, -2.6885701560685891),
               tolerance = 1e-12)
  expect_identical(attr(J, "evaluations"), 8L)

  pair <- community_solve_singularity(dd99_newton(), x0 = trait_matrix(c(-0.3, 0.6), "x"),
                                      birth_rate = c(300, 300))
  expect_true(attr(pair, "converged"))
  expect_equal(as.numeric(pair$traits), c(-0.4420268727422344, 0.44202688729495382), tolerance = 1e-12)
  expect_identical(attr(pair, "evaluations"), 10L)

  cl <- community_classify_singularity(pair)
  expect_equal(as.numeric(cl$jacobian),
               c(-2.5329623352992989, 1.9609547023982508, 1.9609546877876345, -2.5329623443894227),
               tolerance = 1e-12)
  expect_equal(as.numeric(cl$invasion_fitness), c(0.9130434786348407, 0.91304347751611148), tolerance = 1e-12)
  expect_identical(cl$evaluations, 6L)
  expect_equal(cl$classification, "branching point")
})

test_that("the parameter Jacobian and the canonical branching rate are pinned sequentially", {
  off <- dd99_newton_off()
  G <- community_selection_gradient_parameter_jacobian(off, community_parameter_map(off, c("x0", "sigma_K")))
  expect_equal(as.numeric(G), c(0.64585173892039993, 0.53607837136718239, -0.70624586226180475, 1.0608255501970212),
               tolerance = 1e-12)
  expect_identical(attr(G, "evaluations"), 4L)

  one <- dd99_newton() |>
    community_add(trait_matrix(0, "x"), birth_rate = 300) |>
    community_demography()
  br <- canonical_branch_rate(community_clear_residents(one), one, community_trait_transform(one),
                              canonical_control(), matrix(0, 1, 1), 1L, 1, 0.05, 1)
  expect_equal(br$rate, 0.01299719785847951, tolerance = 1e-12)
  expect_identical(br$evaluations, 12L)
  expect_identical(br$failed, 0L)
})
