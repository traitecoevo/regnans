# plant_default_assembly_pars() sets biological parameters through
# `strategy_default$pars` (plant #410). Assigning the flat name -- e.g.
# `p$strategy_default$a_l1 <- 2.17` -- does not error; it attaches an inert list
# element to the R-side object, so the value never reaches the strategy and runs
# silently use plant's defaults (#45). These assertions are cheap and catch that
# regression, which is otherwise invisible.

test_that("plant_default_assembly_pars() actually sets the strategy parameters", {
  defaults <- plant::base_parameters("FF16")$strategy_default$pars

  p <- plant_default_assembly_pars()

  # The values the helper exists to set, read back through $pars.
  expect_equal(p$strategy_default$pars$a_l1, 2.17)
  expect_equal(p$strategy_default$pars$a_l2, 0.5)
  expect_equal(p$strategy_default$pars$hmat, 10)

  # ...and they must differ from plant's defaults, or the test proves nothing.
  expect_false(isTRUE(all.equal(defaults$a_l1, 2.17)))
  expect_false(isTRUE(all.equal(defaults$a_l2, 0.5)))
  expect_false(isTRUE(all.equal(defaults$hmat, 10)))

  # hmat is an argument, so it must track it.
  expect_equal(plant_default_assembly_pars(hmat = 15)$strategy_default$pars$hmat, 15)

  # Parameters-level fields (not strategy ones) were always fine.
  expect_equal(p$max_patch_lifetime, 60)
  expect_equal(plant_default_assembly_pars(max_patch_lifetime = 30)$max_patch_lifetime, 30)
})

test_that("plant_default_assembly_pars(model = \"TF24\") builds TF24 and keeps its allometry", {
  defaults <- plant::base_parameters("TF24")$strategy_default$pars
  p <- plant_default_assembly_pars(hmat = 6, max_patch_lifetime = 6, model = "TF24")
  expect_s3_class(p, "Parameters<TF24,TF24_Env>")
  expect_equal(p$strategy_default$pars$hmat, 6)
  expect_equal(p$max_patch_lifetime, 6)
  expect_equal(max(p$node_schedule_times_default), max(plant::node_schedule_times_default(6)))
  # FF16's assembly allometry is not TF24's
  expect_equal(p$strategy_default$pars$a_l1, defaults$a_l1)
  expect_equal(p$strategy_default$pars$a_l2, defaults$a_l2)
})

test_that("a community without a harness takes the plant model its parameters were built for", {
  # The harness's defaults differ by model, so a TF24 community given FF16's
  # would iterate an equilibrium that cycles.
  tf24 <- list(p = plant_default_assembly_pars(hmat = 6, max_patch_lifetime = 6, model = "TF24"),
               plant_control = plant_default_assembly_control())
  comm <- community_start(bounds(lma = c(0.02, 0.6)), model_support = tf24)
  expect_identical(comm$harness$model, "TF24")
  expect_identical(comm$harness$fitness_gradient, "model")
  expect_identical(comm$demography_control$equilibrium_solver_name, "equilibrium_solve_newton")
  ff16 <- community_start(bounds(lma = c(0.02, 0.6)), model_support = assembly_model_support())
  expect_identical(ff16$harness$model, "FF16")
  expect_identical(ff16$demography_control$equilibrium_solver_name, "equilibrium_iteration")
  expect_error(community_start(bounds(lma = c(0.02, 0.6)), model_support = tf24, harness = harness_plant("FF16")),
               "was given TF24 parameters")
  expect_null(plant_parameters_model(NULL))
  expect_null(plant_parameters_model(list(a = 1)))
})

test_that("plant_default_assembly_pars(fixed_RA = TRUE) sets the allocation pars", {
  p <- plant_default_assembly_pars(fixed_RA = TRUE)
  expect_equal(p$strategy_default$pars$a_f1, 0.5)
  expect_equal(p$strategy_default$pars$a_f2, 0)

  # Off by default.
  q <- plant_default_assembly_pars()
  expect_equal(q$strategy_default$pars$a_f1,
               plant::base_parameters("FF16")$strategy_default$pars$a_f1)
})

test_that("the parameters reach a generated strategy, not just the defaults", {
  # The end-to-end check: what the SCM would actually run with.
  p <- plant_default_assembly_pars()
  s <- plant::generate_strategy(p, trait_matrix(0.0825, "lma"), birth_rate = 1)[[1]]
  expect_equal(s$pars$a_l1, 2.17)
  expect_equal(s$pars$a_l2, 0.5)
  expect_equal(s$pars$hmat, 10)
})
