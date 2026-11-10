# Plant SCM anchors for the coexistence boundaries (R/boundary.R): the label
# and the scouting map, which carries only numbers across to its workers, and
# the square Newton solve for a boundary point from a bracketed seed.
# Reference-free, as test-plant-smoke-singularity.R: the assertions are on the
# end state's own properties, not pinned values. Slow (an assembly is a
# singularity solve, a classification and the landscape's mutant sweeps, each
# solve a few seconds), so the designs are two points and the cap one
# resident where the count is not the point.

# At file level, not inside the test: a map defined beside a solved plant
# community would carry it, and its SCM pointer, to the workers in its frame.
# The mutant gradient is the finite difference unless asked for plant's
# sweep, which on FF16 costs about three times as much (plant#729). The
# parameter step is coarse, as the singularity's continuation takes it
# (test-plant-smoke-singularity.R); the label and the scout never use it.
by_hmat <- function(hmat, fitness_gradient = "fd") {
  community_start(bounds(lma = c(0.02, 0.6)),
                  model_support = list(
                    p = plant_default_assembly_pars(hmat = hmat, max_patch_lifetime = 30),
                    plant_control = plant_default_assembly_control()),
                  harness = harness_plant(fitness_gradient = fitness_gradient),
                  derivative_control = list(d_second = 1e-2, eps_second = 1e-2,
                                            d_parameter = 1e-2, eps_parameter = 1e-2))
}
# residual tolerances on plant sit above the steps in its gradient
# (plant#653), as in test-plant-smoke-singularity.R; fitness_tol above the
# mutant-fitness noise
smoke_control <- boundary_control(list(max_residents = 1, evolve_tol = 0.3, fitness_tol = 1e-3))

test_that("the label assembles the SCM's lma community (SCM)", {
  one <- community_label(by_hmat(10), smoke_control)
  l <- attr(one, "label")
  expect_equal(l$n, 1L)
  expect_true(l$status %in% c("stable", "capped"))
  expect_gt(l$evaluations, 2L)
  # the resident is at a singular strategy, strictly inside the bounds
  x <- as.numeric(one$traits)
  expect_gt(x, 0.02)
  expect_lt(x, 0.6)
  cl <- community_classify_singularity(one)
  expect_lt(abs(as.numeric(cl$selection_gradient) / as.numeric(cl$jacobian)), 0.01 * x)
  expect_true(cl$classification %in% c("CSS", "branching point"))
  # the status is the landscape's: capped means a second species could
  # establish somewhere in the box (an ESS resident can still be invaded
  # away from itself), stable that nothing can
  peak <- attr(max_fitness(one), "fitness")
  expect_equal(l$status == "capped", peak > smoke_control$fitness_tol)
})

test_that("scouting hmat labels each point on a separate process, as sequentially (SCM)", {
  # only the community without its solve crosses; the solved one holds the
  # SCM, which no other process can use
  design <- matrix(c(9, 11), ncol = 1, dimnames = list(NULL, "hmat"))
  sequential <- boundary_scout(by_hmat, design, smoke_control)
  expect_equal(sequential$n, c(1L, 1L))
  expect_true(all(sequential$status %in% c("stable", "capped")))
  parallel <- with_two_workers("multisession", boundary_scout(by_hmat, design, smoke_control))
  expect_identical(as.data.frame(parallel), as.data.frame(sequential))
})

test_that("Newton on the augmented system lands on the lma boundary along hmat from a bracketed seed (SCM)", {
  # Along hmat the lma community goes from two species to one. Bisection
  # (bisect_tol 0.05, 98 equilibrium solves) bracketed the change between
  # 15.5 and 15.875, the midpoint's label unresolved: at the boundary the
  # newcomer establishes at fitness_tol and the rounds cycle. The bisection
  # is not re-run; its ends seed the solve, the one-species side's resident
  # and the newcomer from the two-species side. The answer is held inside
  # that bracket, 2.4% of hmat wide, and to the labels either side of it: a
  # parameterisation that moves the boundary out of the bracket fails here,
  # as it should. Through plant's sweep, so every row of the residual (the
  # resident's selection gradient, the mutant's gradient and its fitness) and
  # the gradient row of y*'s Jacobian column come from it.
  bracket <- c(15.5, 15.875)
  seed <- by_hmat(mean(bracket), fitness_gradient = "model") |>
    community_add(trait_matrix(0.0775, "lma"), birth_rate = 1.14)
  by_model <- function(hmat) by_hmat(hmat, fitness_gradient = "model")
  at <- community_solve_boundary(seed, by_model, y = 0.23, p = c(hmat = mean(bracket)), tol = 1e-3)
  b <- attr(at, "boundary")
  expect_true(b$converged)
  expect_equal(b$type, "invasion")
  expect_equal(b$free, "hmat")
  expect_gt(b$p[["hmat"]], bracket[1])
  expect_lt(b$p[["hmat"]], bracket[2])
  # every residual to its tolerance; plant's noise stops Newton near 1e-4
  expect_lt(max(abs(b$residual)), 1e-3)
  expect_equal(names(b$residual), c("g:lma", "ds/dlma(y)", "s(y)"))
  # the resident and the marginal mutant inside the bounds, the mutant a
  # strict peak away from the resident, which is itself an ESS: the
  # two-species side arises by invasion, not by branching
  x <- as.numeric(at$traits)
  expect_true(x > 0.02 && x < 0.6 && b$y > 0.02 && b$y < 0.6)
  expect_gt(b$y / x, 2)
  expect_lt(community_fitness_hessian(at, matrix(b$y, nrow = 1L))[1, 1], 0)
  expect_lt(community_fitness_hessian(at, at$traits)[1, 1], 0)
  # a few Newton steps: one Jacobian (2mk + 2 = 4 solves) at the start,
  # Broyden-updated after, and a residual per step
  expect_gte(b$iterations, 1L)
  expect_lte(b$iterations, 6L)
  expect_lte(b$evaluations, 30L)
  # 2% either side the labels are two species and one, both the model's;
  # the labels take the finite-difference gradient, where it is cheaper
  label <- function(hmat) attr(community_label(by_hmat(hmat), boundary_control(
    modifyList(smoke_control, list(max_residents = 2)))), "label")[c("n", "status")]
  expect_equal(label(b$p[["hmat"]] * 0.98), list(n = 2L, status = "stable"))
  expect_equal(label(b$p[["hmat"]] * 1.02), list(n = 1L, status = "stable"))
})
