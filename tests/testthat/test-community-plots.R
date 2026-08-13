# community_plot_fitness_landscape.
#
# The function was previously unusable: it demanded a
# `fitness_surrogate_function` (which only the bayesopt landscape method
# creates), and read a column literally named "x" from fitness_points, whose
# first column is actually named after the trait. Neither held for a grid
# landscape. These tests run on the fast DD99 harness and force the plot to be
# built, since a ggplot object only evaluates its aesthetics at build time --
# constructing one proves nothing.

dd99_comm <- function(trait_name = "x", scale = "linear",
                      trait_bounds = c(-2, 2), x0 = 0) {
  b <- matrix(trait_bounds, nrow = 1,
              dimnames = list(trait_name, c("lower", "upper")))
  community_start(b, trait_scale = scale,
                  harness = harness_dd99(x0 = x0, sigma_K = 1, sigma_C = 0.4,
                                         trait_name = trait_name),
                  fitness_control = list(n_evals = 20))
}

built <- function(p) {
  expect_s3_class(p, "ggplot")
  ggplot2::ggplot_build(p)
}

test_that("community_plot_fitness_landscape plots a grid landscape", {
  comm <- dd99_comm() |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography() |>
    community_fitness_landscape()

  b <- built(community_plot_fitness_landscape(comm))
  # zero line, landscape line, sampled points, residents
  expect_equal(length(b$data), 4L)
  # the resident layer carries the single resident, at fitness ~0
  residents <- b$data[[4]]
  expect_equal(nrow(residents), 1L)
  expect_equal(residents$x, 0.5, tolerance = 1e-8)
  expect_lt(abs(residents$y), 1e-6)
})

test_that("community_plot_fitness_landscape computes the landscape if needed", {
  comm <- dd99_comm() |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100)
  expect_null(comm$fitness_points)
  b <- built(community_plot_fitness_landscape(comm))
  expect_gt(nrow(b$data[[2]]), 1L)
})

test_that("community_plot_fitness_landscape handles an empty community", {
  b <- built(community_plot_fitness_landscape(dd99_comm()))
  expect_equal(length(b$data), 3L)   # no resident layer
})

test_that("community_plot_fitness_landscape works for a trait not called x", {
  # the old version indexed fitness_points[["x"]], which does not exist here
  comm <- dd99_comm(trait_name = "lma", scale = "log",
                    trait_bounds = c(0.05, 2), x0 = 0.5) |>
    community_add(trait_matrix(0.4, "lma"), birth_rate = 100) |>
    community_demography()
  p <- community_plot_fitness_landscape(comm)
  expect_equal(p$labels$x, "lma")
  b <- built(p)
  expect_true(all(is.finite(b$data[[2]]$y)))
})

test_that("community_plot_fitness_landscape takes a label and limits", {
  comm <- dd99_comm() |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography()
  b <- built(community_plot_fitness_landscape(comm, label = "step 3",
                                              xlim = c(-1, 1),
                                              ylim = c(-0.5, 0.5)))
  expect_equal(length(b$data), 5L)   # + the annotation layer
  # coord_cartesian zooms rather than filtering, so no point is dropped
  expect_gte(nrow(b$data[[2]]), 20L)
})

test_that("community_plot_fitness_landscape rejects a multi-trait community", {
  comm <- community_start(bounds(x1 = c(-2, 2), x2 = c(-2, 2)),
                          trait_scale = "linear", harness = harness_dd99_nd())
  expect_error(community_plot_fitness_landscape(comm), "single-trait")
})
