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

  for (solver in c("equilibrium_solve_newton", "equilibrium_solve_nleqslv",
                   "equilibrium_solve_dfsane", "equilibrium_hybrid")) {
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
  # Coarser second-derivative steps than the default: the SCM's fitness is only
  # as smooth as its cohort schedule, so a 1e-3 relative step sits in the noise.
  # The selection gradient itself jumps: on the log scale it falls at about
  # -740 per unit log(lma) near the root, but steps up by ~0.23 every ~2e-5 in
  # relative lma (reproduced from cold solves, so the model's, not the warm
  # starts'; plant#653). A residual tolerance below the jumps can leave the search on one
  # of them, which #74 rightly reports as not converged; 0.3 sits above them
  # and still places the root to ~4e-4 of log(lma).
  out <- community_start(bounds(lma = c(0.02, 0.6)),
                         model_support = assembly_model_support(),
                         derivative_control = list(d_second = 1e-2,
                                                   eps_second = 1e-2)) |>
    community_solve_singularity(x0 = 0.08, tol = 0.3)

  expect_true(attr(out, "converged"))
  root <- as.numeric(attr(out, "singularity"))
  expect_gt(root, 0.02)                 # strictly inside the bracket, so this
  expect_lt(root, 0.6)                  # is a real root, not a clamped edge
  expect_equal(as.numeric(out$traits), root, tolerance = 1e-8)

  cl <- community_classify_singularity(out)
  expect_s3_class(cl, "singularity_classification")
  expect_equal(dim(cl$hessian), c(1L, 1L))
  expect_equal(dim(cl$jacobian), c(1L, 1L))
  expect_true(is.finite(cl$hessian))
  expect_true(is.finite(cl$jacobian))
  expect_equal(names(cl$traits), "lma")

  # The defining property of a singular strategy, stated scale-free: the Newton
  # step the remaining gradient implies, |g / J|, is a negligible fraction of
  # the trait value. (An absolute threshold on the gradient would not survive a
  # reparameterisation -- the curvature here spans several orders of magnitude
  # between FF16 parameterisations.)
  newton_step <- abs(as.numeric(cl$selection_gradient) / as.numeric(cl$jacobian))
  expect_lt(newton_step, 0.01 * root)

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

  # Under a plan the Jacobian's stencil points are solved by workers, each
  # from the root's densities as they are without one, so the answer is the
  # same. Only the community without its solve goes to a worker; the solved
  # one holds the SCM, which no other process can use, so this runs under a
  # plan whose workers are separate processes too.
  expect_gt(regnans_pointers(out), 0L)
  expect_identical(regnans_pointers(community_clear_residents(out)), 0L)
  for (type in c("multicore", "multisession")) {
    parallel <- with_two_workers(type, community_classify_singularity(out))
    expect_identical(parallel$jacobian, cl$jacobian)
    expect_identical(parallel$evaluations, cl$evaluations)
    expect_identical(parallel$classification, cl$classification)
  }
})

test_that("a singularity is followed along a plant parameter, its sensitivity matching the roots (SCM)", {
  # hmat, the height at maturation, through a caller-written p -> community.
  # Reference-free: the sensitivity at the middle point must agree with the
  # central difference of the roots the corrector found either side. The
  # parameter step is coarse for the reason the second-derivative steps are
  # (plant#653); the sensitivity moves by 0.3% between 1e-2 and 5e-2, and the
  # roots' difference agrees with it to ~1%.
  by_hmat <- function(hmat) {
    community_start(bounds(lma = c(0.02, 0.6)),
                    model_support = list(
                      p = plant_default_assembly_pars(hmat = hmat, max_patch_lifetime = 30),
                      plant_control = plant_default_assembly_control()),
                    derivative_control = list(d_second = 1e-2, eps_second = 1e-2,
                                              d_parameter = 1e-2, eps_parameter = 1e-2))
  }
  start <- community_solve_singularity(by_hmat(9), x0 = 0.08, tol = 0.3)
  path <- community_continue_singularity(start, by_hmat, c(9, 10, 11), tol = 0.3)

  expect_null(path$stopped)
  x <- as.numeric(path$traits)
  expect_true(all(x > 0.02 & x < 0.6))
  expect_equal(unname(path$sensitivity[2, 1]), (x[3] - x[1]) / 2, tolerance = 0.05)
  # the predictor is first order: its miss is small beside the step it took
  miss <- abs(path$predicted[-1, 1] - x[-1])
  expect_true(all(miss < 0.25 * abs(diff(x))))
  expect_true(all(path$classification %in% c("CSS", "branching point")))
})

test_that("the resident Jacobian's stencil is solved on its centre's discretisation (SCM)", {
  # Each stencil point would otherwise refine its own node schedule and ODE
  # steps, and the change of steps over the step swamps the derivative: at
  # lma 0.15 the free difference of the selection gradient read -1619, 193
  # and 107 at steps 1e-2, 1e-3 and 1e-4. On the centre's discretisation it is
  # one number to the difference's truncation error.
  comm <- community_start(bounds(lma = c(0.02, 0.6)),
                          model_support = list(p = plant_default_assembly_pars(hmat = 10, max_patch_lifetime = 30),
                                               plant_control = plant_default_assembly_control())) |>
    community_add(trait_matrix(0.15, "lma"), birth_rate = 10) |>
    community_demography()
  disc <- community_discretisation(comm)
  expect_named(disc, c("node_schedule_times", "ode_times", "ode_step_sizes"))
  expect_gt(length(disc$ode_times), 0L)
  frozen <- stencil_community(comm, 1e-3, disc)
  expect_identical(frozen$model_support$frozen, disc)
  expect_equal(frozen$model_support$plant_control$offspring_production_tol, 1e-9)
  # a frozen discretisation survives a change of residents, so a stencil point
  # built by community_add() replays it
  moved <- community_add(community_clear_residents(frozen), trait_matrix(0.151, "lma"), birth_rate = 80)
  expect_identical(plant_community_parameters(moved)$ode_times, disc$ode_times)

  jacobian <- function(d) {
    comm$derivative_control <- derivative_control(list(d_second = d, eps_second = d))
    community_selection_gradient_jacobian(comm)[1, 1]
  }
  coarse <- jacobian(1e-2)
  fine <- jacobian(1e-3)
  expect_lt(coarse, 0)
  expect_lt(fine, 0)
  expect_equal(coarse, fine, tolerance = 0.25)
})

test_that("plant's own mutant gradient agrees with the finite difference (SCM)", {
  # plant#650 E2: the mutants' replay in the recorded canopy, differentiated
  # by plant's reverse sweep through the hyperparameter function. The finite
  # difference of the same replay agrees to its truncation error: 7e-5 and
  # 9e-4 away from the resident, 5e-3 at it, where the landscape is steepest
  # (a selection gradient of -1e4); the Richardson reference of the contract
  # check agrees to 1e-2.
  build <- function(source) {
    community_start(bounds(lma = c(0.02, 0.6)),
                    model_support = list(p = plant_default_assembly_pars(hmat = 10, max_patch_lifetime = 30),
                                         plant_control = plant_default_assembly_control()),
                    harness = harness_plant(fitness_gradient = source)) |>
      community_add(trait_matrix(0.15, "lma"), birth_rate = 10) |>
      community_demography()
  }
  model <- build("model")
  expect_true(harness_provides(model$harness, "fitness_gradient"))
  expect_equal(unname(model$fitness_derivatives$source[["fitness_gradient"]]), "model")
  y <- c(0.1, 0.15, 0.3)
  g <- community_fitness_gradient(model, y)
  fd <- community_fitness_gradient(build("fd"), y)
  expect_equal(as.numeric(g), as.numeric(fd), tolerance = 1e-2)
  expect_equal(attr(g, "value"), as.numeric(model$fitness_function(y)), tolerance = 1e-10)
  checked <- harness_check_derivatives(model, n_points = 3, tol_rel = 1e-2)
  expect_true(all(checked$pass))
  expect_error(harness_plant("TF24", fitness_gradient = "model"), "FF16's fitness only")
})
