# Tests for the assembly loop (assembler_control/start/run), output tidying
# (tidy_assembly) and the mutation helper (mutational_vcv_proportion).

# ---- assembler_control (fast, no SCM) --------------------------------------

test_that("assembler_control has sane defaults", {
  ctrl <- assembler_control()
  expect_equal(ctrl$run_type, "to_equilibrium")
  expect_equal(ctrl$birth_type, "maximum")
  expect_true(ctrl$compute_viable_fitness)
})

test_that("assembler_control switches run_type for stochastic births", {
  ctrl <- assembler_control(list(birth_type = "stochastic", vcv = diag(1)))
  expect_equal(ctrl$run_type, "single_step")
})

test_that("assembler_control validates its arguments", {
  expect_error(assembler_control(list(not_a_param = 1)),
               "Unknown control parameters")
  expect_error(assembler_control(list(birth_type = "stochastic")),
               "vcv must be provided")
  expect_error(
    assembler_control(list(birth_type = "maximum", run_type = "single_step")),
    "Must use 'to_equilibrium' run type")
})

# ---- mutational_vcv_proportion (fast, no SCM) ------------------------------

test_that("mutational_vcv_proportion builds a diagonal vcv on log scale", {
  b <- bounds(lma = c(0.01, 2))
  v <- mutational_vcv_proportion(b, p = 0.01)
  expect_true(is.matrix(v))
  expect_equal(dim(v), c(1L, 1L))
  # diagonal entry = p * (log(upper) - log(lower))
  expect_equal(v[1, 1], 0.01 * (log(2) - log(0.01)))

  # off-diagonals are zero for multiple independent traits
  v2 <- mutational_vcv_proportion(bounds(lma = c(0.01, 2), rho = c(1, 1000)),
                                  p = 0.02)
  expect_equal(dim(v2), c(2L, 2L))
  expect_equal(v2[1, 2], 0)
  expect_equal(v2[2, 1], 0)
})

test_that("mutational_vcv_proportion accepts a community and rejects infinite bounds", {
  comm <- community_start(bounds(lma = c(0.01, 2)))
  expect_equal(mutational_vcv_proportion(comm, p = 0.01),
               mutational_vcv_proportion(comm$bounds, p = 0.01))
  expect_error(mutational_vcv_proportion(bounds_infinite("lma")),
               "must be finite")
})

test_that("mutational_vcv_proportion uses the community's trait scale", {
  # A bare bounds matrix has no scale attached -> log range (historical default).
  expect_equal(mutational_vcv_proportion(bounds(lma = c(0.01, 2)), p = 0.01)[1, 1],
               0.01 * (log(2) - log(0.01)))
  # A log-scale community matches the bare-bounds default.
  log_comm <- community_start(bounds(lma = c(0.01, 2)), trait_scale = "log")
  expect_equal(mutational_vcv_proportion(log_comm, p = 0.01)[1, 1],
               0.01 * (log(2) - log(0.01)))
  # A linear-scale community uses the linear range instead (issue: stochastic
  # births should not be hardcoded to log).
  lin_comm <- community_start(bounds(x = c(-3, 3)), trait_scale = "linear")
  expect_equal(mutational_vcv_proportion(lin_comm, p = 0.01)[1, 1], 0.01 * 6)
})

# ---- merging residents (fast, no SCM) ---------------------------------------

test_that("two residents within eps_too_close of each other are one, either side of zero", {
  # GK98's backward boundary sweep at d = 1 ended on this pair, collapsed onto
  # the one singular strategy at 0 (test-boundary.R). Two finite-difference
  # steps (relative, so 4.6e-8 here) never reach across zero; the distance at
  # which a birth beside a resident is that resident does.
  harness <- harness_gk98(sigma = 1, mu = c(-1, 0, 1), K = c(1, 1, 1))
  pair <- function(x, b = bounds(x = c(-4, 4))) {
    community_start(b, trait_scale = "linear", harness = harness) |>
      community_add(trait_matrix(x, "x"), birth_rate = c(2, 1))
  }
  collapsed <- pair(c(-2.32e-5, 2.32e-5))
  expect_equal(evolve_merged(collapsed), integer(0))
  expect_equal(evolve_merged(collapsed, eps_too_close = 1e-3), 2L)
  # 0.2 apart is 2.5% of the width
  expect_equal(evolve_merged(pair(c(-0.1, 0.1)), eps_too_close = 1e-3), integer(0))
  # with no finite width only the steps decide
  expect_equal(evolve_merged(pair(c(-2.32e-5, 2.32e-5), bounds_infinite("x")), eps_too_close = 1e-3),
               integer(0))
})

test_that("only a pair the canonical route ends on is merged by eps_too_close", {
  # A root Newton solved keeps a close pair (a pitchfork's is genuinely that
  # close); a canonical run's end may still be closing. Each route is stood in
  # for by one that ends on the collapsed pair.
  collapsed <- community_start(bounds(x = c(-4, 4)), trait_scale = "linear",
                               harness = harness_gk98(sigma = 1, mu = c(-1, 0, 1), K = c(1, 1, 1))) |>
    community_add(trait_matrix(c(-2.32e-5, 2.32e-5), "x"), birth_rate = c(2, 1))
  ends_on <- function(route) {
    function(community, ctrl) {
      list(community = community, converged = TRUE, classification = NA_character_, reason = NULL,
           evaluations = 0L, route = route)
    }
  }
  obj <- list(control = assembler_control(list(evolve = "newton")), community = collapsed, evaluations = 0L)
  local_mocked_bindings(evolve_newton = ends_on("newton"))
  expect_equal(nrow(assembler_evolve(obj)$community$traits), 2L)
  local_mocked_bindings(evolve_newton = ends_on("canonical"))
  kept <- assembler_evolve(obj)$community
  # the commoner resident stays, its density solved again alone
  expect_equal(as.numeric(kept$traits), -2.32e-5)
})

# ---- assembler loop (integration) ------------------------------------------
# Driven on the fast DD99 harness (issue #27) rather than the SCM, on a linear
# trait scale (DD99's optimum sits at x0 = 0). The stochastic-births path now
# mutates on whichever trait scale the community uses, so it runs on these
# negative-spanning bounds too (previously it was hardcoded to log).

dd99_assembler_community <- function() {
  community_start(
    bounds(x = c(-3, 3)),
    harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4),
    trait_scale = "linear",
    fitness_control = list(method = "grid", n_evals = 40))
}

test_that("assembler_run drives maximum-fitness assembly", {
  control <- assembler_control(list(run_type = "to_equilibrium",
                                    birth_type = "maximum",
                                    birth_move_tol = 1,
                                    compute_viable_fitness = FALSE))

  obj <- assembler_start(dd99_assembler_community(), control = control)
  expect_s3_class(obj, "assembler")
  expect_length(obj, 1L)                # history holds the initial community

  obj <- assembler_run(obj, nsteps = 3)
  expect_length(obj, 4L)                # initial + 3 steps
  expect_gte(length(obj$community), 1L)

  # all residents lie within the trait bounds
  x <- obj$community$traits[, "x"]
  expect_true(all(x >= -3 & x <= 3))
  expect_true(all(obj$community$birth_rate > 0))
})

test_that("an assembler saved before the evolution and refusal controls restores and runs", {
  control <- assembler_control(list(birth_move_tol = 1, compute_viable_fitness = FALSE))
  old <- assembler_start(dd99_assembler_community(), control = control)
  old$control[c("evolve", "evolve_tol", "branch_offset", "max_residents", "edge_tol")] <- NULL
  old$evaluations <- NULL
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(old, path)
  restored <- assembler_start(dd99_assembler_community(), control = control, filename = path)
  expect_identical(restored$control$max_residents, Inf)
  expect_identical(restored$control$evolve, "none")
  expect_identical(restored$evaluations, 0L)
  out <- assembler_run(restored, nsteps = 1)
  expect_length(out, 2L)
  expect_gt(out$evaluations, 0L)
})

test_that("assembler_run drives stochastic assembly on a linear trait scale", {
  set.seed(42)
  control <- assembler_control(list(run_type = "single_step",
                                    birth_type = "stochastic",
                                    vcv = diag(1) * 0.1))

  obj <- assembler_start(dd99_assembler_community(), control = control)
  obj <- assembler_run(obj, nsteps = 5)
  expect_s3_class(obj, "assembler")
  expect_length(obj, 6L)                # initial + 5 steps

  if (length(obj$community) > 0L) {
    x <- obj$community$traits[, "x"]
    expect_true(all(x >= -3 & x <= 3))
  }
})

# ---- tidy_assembly ---------------------------------------------------------

test_that("tidy_assembly tidies an assembly run into one row per resident-step", {
  control <- assembler_control(list(run_type = "to_equilibrium",
                                    birth_type = "maximum",
                                    birth_move_tol = 1,
                                    compute_viable_fitness = FALSE))
  obj <- assembler_start(dd99_assembler_community(), control = control) |>
    assembler_run(nsteps = 2)

  ta <- tidy_assembly(obj)
  expect_s3_class(ta, "tbl_df")
  expect_true(all(c("step", "strategy_id", "resident", "births", "traits",
                    "fitness_landscape") %in% names(ta)))
  expect_type(ta$strategy_id, "character")
  expect_true(all(ta$resident))
  expect_true(all(ta$births > 0))
  # traits and fitness_landscape are list columns of tibbles
  expect_true(is.list(ta$traits))
  expect_s3_class(ta$traits[[1]], "tbl_df")
  expect_s3_class(ta$fitness_landscape[[1]], "tbl_df")
})

test_that("maximum-fitness births refuse a community without bounds instead of stopping in the debugger", {
  control <- assembler_control(list(birth_type = "maximum"))
  expect_error(regnans:::community_new_types_maximum_fitness(list(bounds = NULL, trait_names = "x"), control),
               "need trait bounds")
})
