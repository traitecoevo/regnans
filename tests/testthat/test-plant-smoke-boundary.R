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
