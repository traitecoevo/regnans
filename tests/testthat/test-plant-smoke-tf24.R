# TF24 through the plant harness. The resident (lma 0.0825, hmat 6, patch
# lifetime 6, plant's default TF24 environment) is viable and short enough to
# solve in seconds: its log R0 is about +3 at a birth rate of 10.

tf24_support <- function() {
  list(p = plant_default_assembly_pars(hmat = 6, max_patch_lifetime = 6, model = "TF24"),
       plant_control = plant_default_assembly_control())
}

tf24_resident <- function(fitness_gradient = NULL, birth_rate = 54) {
  community_start(bounds(lma = c(0.02, 0.6)), model_support = tf24_support(),
                  harness = harness_plant("TF24", fitness_gradient = fitness_gradient)) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = birth_rate) |>
    community_demography()
}

test_that("a TF24 resident's equilibrium is found by Newton's method (SCM)", {
  # Newton is TF24's default because fixed-point iteration from a birth rate
  # of 10 was measured to settle on a two-cycle (about 1.45 and 527: the sparse
  # stand floods the patch with seed and the dense one crashes). That is not
  # re-run here; the solver's name and the equilibrium are asserted.
  comm <- tf24_resident(birth_rate = 10)
  expect_equal(comm$demography_control$equilibrium_solver_name, "equilibrium_solve_newton")
  expect_true(attr(comm, "converged"))
  expect_equal(as.numeric(comm$birth_rate), 53.972, tolerance = 1e-3)
  expect_equal(as.numeric(comm$resident_fitness), 0, tolerance = 1e-4)
})

test_that("plant's TF24 mutant gradient is the default and agrees with the finite difference (SCM)", {
  # The sweep is TF24's default because it was measured to cost what the
  # difference does at one trait and less from two (not timed here). It agreed
  # with the finite difference of the same replay to 6e-5 on macOS; it is held
  # to 1e-3, the difference's own truncation error on other platforms.
  model <- tf24_resident()
  expect_identical(model$harness$fitness_gradient, "model")
  expect_true(harness_provides(model$harness, "fitness_gradient"))
  y <- c(0.06, 0.0825, 0.11)
  g <- community_fitness_gradient(model, y)
  fd <- community_fitness_gradient(tf24_resident("fd", birth_rate = model$birth_rate), y)
  expect_equal(as.numeric(g), as.numeric(fd), tolerance = 1e-3)
  expect_equal(attr(g, "value"), as.numeric(model$fitness_function(y)), tolerance = 1e-10)
  # the sweep's replay shares the fitness function's SCM and leaves it as found
  expect_false(environment(model$fitness_function)$scm$keep_states)
  checked <- harness_check_derivatives(model, n_points = 3, tol_rel = 1e-2)
  expect_true(all(checked$pass))
})

test_that("TF24's implicit resident Jacobian has the stencil's sign and size (SCM)", {
  # With each fixed-environment run refining its own schedule, the implicit
  # Jacobian here read +733 against the stencil's -486 at a step of 1e-3: the
  # wrong sign, so a convergence-stable point would be classified unstable.
  # Replaying the centre's discretisation it was measured at -489; held to 2e-2.
  comm <- tf24_resident()
  jacobian <- function(route) {
    comm$derivative_control <- derivative_control(list(resident_jacobian = route, d_second = 1e-3, eps_second = 1e-3))
    community_selection_gradient_jacobian(comm)[1, 1]
  }
  stencil <- jacobian("stencil")
  implicit <- jacobian("implicit")
  expect_lt(stencil, 0)
  expect_equal(unname(implicit), unname(stencil), tolerance = 2e-2)
})
