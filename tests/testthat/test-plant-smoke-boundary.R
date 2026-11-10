# Plant SCM anchor for the coexistence label (R/boundary.R): the real model
# drives the assembly, and the scouting map carries only numbers across to
# its workers. Reference-free, as test-plant-smoke-singularity.R: the
# assertions are on the end state's own properties, not pinned values. Slow
# (an assembly is a singularity solve, a classification and the landscape's
# mutant sweeps, each solve a few seconds), so the design is two points and
# the cap one resident.

# At file level, not inside the test: a map defined beside a solved plant
# community would carry it, and its SCM pointer, to the workers in its frame.
by_hmat <- function(hmat) {
  community_start(bounds(lma = c(0.02, 0.6)),
                  model_support = list(
                    p = plant_default_assembly_pars(hmat = hmat, max_patch_lifetime = 30),
                    plant_control = plant_default_assembly_control()),
                  derivative_control = list(d_second = 1e-2, eps_second = 1e-2))
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

# The same community through plant's own mutant gradient, so the residual's
# ds/dy(y) row comes from the sweep, and a coarser parameter step, as the
# singularity's continuation takes (test-plant-smoke-singularity.R).
by_hmat_model <- function(hmat) {
  community_start(bounds(lma = c(0.02, 0.6)),
                  model_support = list(
                    p = plant_default_assembly_pars(hmat = hmat, max_patch_lifetime = 30),
                    plant_control = plant_default_assembly_control()),
                  harness = harness_plant(fitness_gradient = "model"),
                  derivative_control = list(d_second = 1e-2, eps_second = 1e-2,
                                            d_parameter = 1e-2, eps_parameter = 1e-2))
}

test_that("Newton on the augmented system lands on the lma boundary along hmat from a bracketed seed (SCM)", {
  # Along hmat the lma community goes from two species to one: bisection
  # (bisect_tol 0.05, 98 equilibrium solves) bracketed the change between
  # 15.5 and 15.875 on 2026-11-10, with the midpoint's label unresolved (the
  # rounds cycled: at the boundary the newcomer's establishment is at
  # fitness_tol). The bisection is not re-run; its ends seed the solve, the
  # one-species side's resident and the newcomer from the two-species side.
  # Reference-free: the answer is held inside the bracket, with the labels
  # either side of it, not to a value.
  bracket <- c(15.5, 15.875)
  seed <- by_hmat_model(mean(bracket)) |>
    community_add(trait_matrix(0.0775, "lma"), birth_rate = 1.14)
  at <- community_solve_boundary(seed, by_hmat_model, y = 0.23, p = c(hmat = mean(bracket)), tol = 1e-3)
  b <- attr(at, "boundary")
  expect_true(b$converged)
  expect_equal(b$type, "invasion")
  expect_gt(b$p[["hmat"]], bracket[1])
  expect_lt(b$p[["hmat"]], bracket[2])
  # every residual to its tolerance; plant's noise stops Newton near 1e-4
  expect_lt(max(abs(b$residual)), 1e-3)
  expect_equal(names(b$residual), c("g:lma", "ds/dlma(y)", "s(y)"))
  # the marginal mutant is a strict peak away from the resident
  expect_lt(community_fitness_hessian(at, matrix(b$y, nrow = 1L))[1, 1], 0)
  expect_gt(b$y / as.numeric(at$traits), 2)
  # a few Newton steps: each costs one solve and a Jacobian 2mk + 2 = 4
  expect_lte(b$iterations, 6L)
  expect_lte(b$evaluations, 30L)
  # 2% either side the labels are two species and one, both the model's
  label <- function(hmat) attr(community_label(by_hmat_model(hmat), boundary_control(list(
    max_residents = 2, evolve_tol = 0.3, fitness_tol = 1e-3))), "label")[c("n", "status")]
  expect_equal(label(b$p[["hmat"]] * 0.98), list(n = 2L, status = "stable"))
  expect_equal(label(b$p[["hmat"]] * 1.02), list(n = 1L, status = "stable"))
})
