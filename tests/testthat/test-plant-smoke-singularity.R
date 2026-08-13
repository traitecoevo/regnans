# Plant SCM anchors for the machinery added/repaired in this branch.
#
# Two-tier testing (see AGENTS.md): every algorithm here is developed and
# validated against the toy harnesses, where the answers are analytic and a
# solve takes milliseconds -- test-demography-solvers.R for the equilibrium
# solvers, test-singularity.R for the N-D solver and the classifier. This file
# is the second tier: a minimal set of genuine plant-SCM runs confirming the
# real physiological model still drives them.
#
# Deliberately *reference-free*: these tests assert internal consistency
# (solvers agreeing with each other, a gradient vanishing where the solver says
# it does) rather than pinned trait values, so they stay meaningful when the
# plant parameterisation changes. The pinned plant reference values live in
# test-plant-smoke.R.
#
# These are slow (each residual evaluation re-solves a community to demographic
# equilibrium). Keep the count small.

test_that("the alternative equilibrium solvers agree with iteration (SCM)", {
  # equilibrium_solve_nleqslv / equilibrium_solve_dfsane / equilibrium_hybrid
  # had never been run against the real SCM -- only equilibrium_iteration had.
  solve_with <- function(solver) {
    community_start(bounds(lma = c(0.01, 2)),
                    model_support = assembly_model_support(),
                    demography_control = demographic_step_control(
                      list(equilibrium_solver_name = solver))) |>
      community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
      community_demography()
  }

  ref <- solve_with("equilibrium_iteration")
  expect_true(attr(ref, "converged"))

  for (solver in c("equilibrium_solve_nleqslv", "equilibrium_solve_dfsane",
                   "equilibrium_hybrid")) {
    got <- solve_with(solver)
    expect_true(attr(got, "converged"), info = solver)
    expect_equal(as.numeric(got$birth_rate), as.numeric(ref$birth_rate),
                 tolerance = 1e-3, info = solver)
    # at equilibrium the resident just replaces itself
    expect_equal(as.numeric(got$resident_fitness), 0, tolerance = 1e-3,
                 info = solver)
  }
})

test_that("community_solve_singularity and the classifier run on the SCM", {
  # A bracket wide enough to hold the lma attractor under any current FF16
  # parameterisation; the assertions below are on the solution's own
  # properties, not on where it lands.
  out <- community_start(bounds(lma = c(0.02, 0.6)),
                         model_support = assembly_model_support()) |>
    community_solve_singularity(x0 = 0.08, tol = 1e-3)

  expect_true(attr(out, "converged"))
  root <- as.numeric(attr(out, "singularity"))
  expect_gt(root, 0.02)                 # strictly inside the bracket, so this
  expect_lt(root, 0.6)                  # is a real root, not a clamped edge
  expect_equal(as.numeric(out$traits), root, tolerance = 1e-8)
  # the defining property of a singular strategy
  expect_lt(abs(as.numeric(out$selection_gradient)), 0.5)

  cl <- community_classify_singularity(out, d = 1e-2, eps = 1e-2)
  expect_s3_class(cl, "singularity_classification")
  expect_equal(dim(cl$hessian), c(1L, 1L))
  expect_equal(dim(cl$jacobian), c(1L, 1L))
  expect_true(is.finite(cl$hessian))
  expect_true(is.finite(cl$jacobian))
  expect_equal(names(cl$traits), "lma")

  # The plant lma singularity is an evolutionary *attractor*: the selection
  # gradient is positive below it and negative above (test-plant-smoke.R makes
  # that check directly), so the Jacobian must be negative.
  expect_lt(as.numeric(cl$jacobian), 0)
  expect_true(cl$convergence_stable)
  expect_true(cl$classification %in%
                c("CSS", "branching point", "repeller", "Garden of Eden"))
  # and the verdict must agree with the matrices it was derived from
  expect_equal(cl$evolutionarily_stable, as.numeric(cl$hessian) < 0)
  expect_equal(cl$classification,
               if (cl$evolutionarily_stable) "CSS" else "branching point")
})
