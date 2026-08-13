# plant_default_assembly_pars() sets biological parameters through
# `strategy_default$pars` (plant #410). Assigning the flat name -- e.g.
# `p$strategy_default$a_l1 <- 2.17` -- does not error; it attaches an inert list
# element to the R-side object, so the value never reaches the strategy and runs
# silently use plant's defaults (#45). These assertions are cheap and catch that
# regression, which is otherwise invisible.

test_that("plant_default_assembly_pars() actually sets the strategy parameters", {
  defaults <- plant::scm_base_parameters("FF16")$strategy_default$pars

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

test_that("plant_default_assembly_pars(fixed_RA = TRUE) sets the allocation pars", {
  p <- plant_default_assembly_pars(fixed_RA = TRUE)
  expect_equal(p$strategy_default$pars$a_f1, 0.5)
  expect_equal(p$strategy_default$pars$a_f2, 0)

  # Off by default.
  q <- plant_default_assembly_pars()
  expect_equal(q$strategy_default$pars$a_f1,
               plant::scm_base_parameters("FF16")$strategy_default$pars$a_f1)
})

test_that("the parameters reach a generated strategy, not just the defaults", {
  # The end-to-end check: what the SCM would actually run with.
  p <- plant_default_assembly_pars()
  s <- plant::generate_strategy(p, trait_matrix(0.0825, "lma"), birth_rate = 1)[[1]]
  expect_equal(s$pars$a_l1, 2.17)
  expect_equal(s$pars$a_l2, 0.5)
  expect_equal(s$pars$hmat, 10)
})
