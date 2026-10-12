# plant's own mutant gradient through the harness contract (plant#650 E2).
#
# Slow, plant-SCM tier (see AGENTS.md): one 30-year resident solved to
# equilibrium, then the reverse sweep of the invasion replay against the
# community's finite-difference reference. The plant control runs on the
# birth-date coordinate, which the sweep needs, and at a tight seed-height root
# tolerance so the difference differentiates the root the sweep does rather
# than the bisection iterate (plant's test-ff16-stand-gradient.R).

skip_if_not_installed("plant")

ad_model_support <- function(max_patch_lifetime = 30) {
  list(p = plant_default_assembly_pars(max_patch_lifetime = max_patch_lifetime),
       plant_control = plant_default_assembly_control(
         node_density_in_birth_date = TRUE,
         offspring_production_tol = 1e-12))
}

test_that("the plant FF16 harness advertises the mutant gradient and it checks out", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = ad_model_support()) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
    community_demography()
  expect_true(attr(comm, "converged"))
  expect_true(harness_provides(comm$harness, "fitness_gradient"))
  expect_identical(unname(comm$fitness_derivatives$source["fitness_gradient"]),
                   "model")

  # The selection gradient from the sweep against a three-level Richardson
  # difference of the same fitness function (replays of the same resident).
  g <- community_fitness_gradient(comm)
  f <- comm$fitness_function
  y <- derivative_points(comm$traits, "lma")
  reference <- modifyList(derivative_control(), list(r_gradient = 3L))
  g_fd <- fd_fitness_gradient(f, y, reference)
  expect_equal(unname(g[1, 1]), unname(g_fd[1, 1]), tolerance = 1e-5)
  expect_equal(attr(g, "value"), attr(g_fd, "value"), tolerance = 1e-10)

  # The contract across the bounds. 1e-3: near the lower bound the SCM's fitness
  # has kinks of the cohort schedule (plant#653) that a 1e-4 difference reads;
  # at a 1e-5 step the difference meets the sweep to 1e-7 there.
  res <- harness_check_derivatives(comm, n_points = 5, tol_rel = 1e-3,
                                   tol_abs = 1e-8)
  expect_true(all(res$pass))
  expect_equal(nrow(res), 5L)
})

test_that("a plant harness without derivatives falls back to differences", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = ad_model_support(),
                          harness = harness_plant("FF16", derivatives = FALSE)) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
    community_demography()
  expect_false(harness_provides(comm$harness, "fitness_gradient"))
  expect_identical(unname(comm$fitness_derivatives$source["fitness_gradient"]),
                   "finite difference")
  expect_identical(nrow(harness_check_derivatives(comm)), 0L)
})
