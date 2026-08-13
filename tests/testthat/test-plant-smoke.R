# Plant SCM smoke tests (issue #27).
#
# The model-agnostic pipeline (demography, selection gradient, fitness
# landscape, attractor solving, assembly, the alternative equilibrium solvers,
# util_nlsolve, inviable-strategy checks) is exercised fast on the DD99 toy
# harness elsewhere. This file keeps a *minimal* set of genuine plant-SCM runs
# as the integration anchor: it is the only place that confirms the real
# physiological model still drives the regnans pipeline, and it pins the
# reference values for the installed plant (2.0.0.9001).
#
# The reference values below were re-blessed when #45 was fixed. Until then
# plant_default_assembly_pars() set a_l1/a_l2/hmat as flat strategy fields, which
# since plant #410 are silently ignored, so these tests had been pinning runs
# with plant's default allometry and hmat rather than the intended assembly
# parameterisation. The values moved a long way (the 1D attractor from lma ~0.14
# to ~0.045); that is the bug being corrected, not drift in plant.
#
# These are the slow tests in the suite (a few seconds per SCM solve,
# max_patch_lifetime = 30 via helper-assembly.R). Keep the count small, and
# prefer cheap gradient checks over full iterative solves where they give the
# same coverage (see the attractor test below).

test_that("empty community: fundamental niche and fitness peak (SCM)", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support())

  # fundamental niche: the positive-fitness interval, strictly inside the bounds
  vb <- community_viable_fitness_1D(comm)
  expect_equal(rownames(vb), "lma")
  expect_equal(colnames(vb), c("lower", "upper"))
  expect_lt(vb[, "lower"], vb[, "upper"])
  # reference values for the installed plant (2.0.0.9001)
  expect_equal(unname(vb[, "lower"]), 0.013091, tolerance = 1e-3)
  expect_equal(unname(vb[, "upper"]), 1.644165, tolerance = 1e-3)

  ff <- plant_community_update_fitness_function(comm)$fitness_function
  expect_gt(ff(exp(mean(log(vb)))), 1)          # clearly positive inside
  expect_lt(abs(ff(vb[, "lower"])), 0.01)       # ~0 at the edges
  expect_lt(abs(ff(vb[, "upper"])), 0.01)

  # max_fitness / max_growth_rate against the same fitness function
  mx <- max_fitness(comm)
  expect_equal(names(mx), "lma")
  expect_gt(as.numeric(mx), 0.01)
  expect_lt(as.numeric(mx), 2)
  g <- max_growth_rate(comm, c(0.05, 0.0825, 0.2))
  expect_equal(g, ff(c(0.05, 0.0825, 0.2)), tolerance = 1e-8)
  grid <- max_growth_rate(comm, seq_log_range(c(0.01, 2), 7))
  expect_gte(attr(mx, "fitness"), max(grid) - 1e-3)
})

test_that("single resident solves to equilibrium, gradient and inviable check (SCM)", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support()) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
    community_demography()

  expect_true(attr(comm, "converged"))
  # equilibrium birth rate, reference for plant 2.0.0.9001 (the genuinely
  # iterative equilibrium_iteration fixed point, robust to the starting rate)
  expect_equal(as.numeric(comm$birth_rate), 112.657365, tolerance = 1e-3)
  expect_equal(comm$resident_fitness, 0, tolerance = 1e-4)

  # selection gradient is finite at the (non-singular) resident
  comm <- community_selection_gradient(comm)
  expect_length(comm$selection_gradient, 1L)
  expect_true(is.finite(comm$selection_gradient))

  # plant_community_check_for_inviable_strategies: the viable resident is kept
  op <- community_check_for_inviable_strategies(comm)
  expect_false(any(attr(op, "drop")))
})

test_that("selection gradient brackets the 1D attractor (SCM)", {
  # The model-agnostic singularity root-finder (community_solve_singularity_1D)
  # is exercised fast on the DD99 toy harness. The full iterative solve against
  # the real SCM cost ~45s (it re-solves a resident to equilibrium on every
  # gradient evaluation). Here we only confirm the genuine plant model produces
  # a *convergent* evolutionary attractor in [0.04, 0.05] (~lma 0.045 for plant
  # 2.0.0.9001): the selection gradient changes sign across the singularity,
  # which brackets it with two equilibrium solves instead of a dozen.
  grad_at <- function(x) {
    comm <- community_start(bounds(lma = c(0.01, 2)),
                            model_support = assembly_model_support()) |>
      community_add(trait_matrix(x, "lma"), birth_rate = 200) |>
      community_demography() |>
      community_selection_gradient()
    as.numeric(comm$selection_gradient)
  }

  g_lo <- grad_at(0.04)
  g_hi <- grad_at(0.05)
  expect_true(is.finite(g_lo))
  expect_true(is.finite(g_hi))
  expect_gt(g_lo, 0)   # below the singularity selection pushes lma up,
  expect_lt(g_hi, 0)   # above it selection pushes lma down -> attractor between
})
