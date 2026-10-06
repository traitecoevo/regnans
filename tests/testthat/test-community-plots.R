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

# ---- plot_community --------------------------------------------------------
#
# plot_community reads the tidy_assembly schema (traits in a list column) and
# plots one step: birth rate against the trait for one trait, the trait plane
# for two. Each plotted point must be a resident of that step.

assembled <- function(comm, nsteps = 4) {
  control <- assembler_control(list(birth_type = "maximum",
                                    compute_viable_fitness = FALSE))
  assembler_start(comm, control) |> assembler_run(nsteps) |> tidy_assembly()
}

residents_at <- function(tidy, s) {
  tidyr::unnest(tidy[tidy$step == s, c("births", "traits")], "traits")
}

test_that("plot_community plots birth rate against one trait", {
  ta <- assembled(community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                                  harness = harness_dd99(x0 = 0, sigma_K = 1,
                                                         sigma_C = 0.4),
                                  fitness_control = list(n_evals = 40)))
  expect_equal(attr(ta, "trait_scale"), "linear")
  last <- max(as.integer(ta$step))
  expected <- residents_at(ta, last)
  expect_gt(nrow(expected), 1L)

  p <- plot_community(ta)
  expect_equal(p$labels$x, "x")
  b <- built(p)
  pts <- b$data[[1]]
  expect_equal(nrow(pts), nrow(expected))
  # linear x axis (negative traits survive); births on a log10 y axis
  expect_equal(sort(pts$x), sort(expected$x))
  expect_equal(sort(pts$y), sort(log10(expected$births)))
  expect_match(b$data[[2]]$label, paste0("Step = ", last))
})

test_that("plot_community picks a step, sharing axes across steps", {
  ta <- assembled(community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                                  harness = harness_dd99(x0 = 0, sigma_K = 1,
                                                         sigma_C = 0.4),
                                  fitness_control = list(n_evals = 40)))
  steps <- sort(unique(as.integer(ta$step)))
  first <- built(plot_community(ta, step = steps[1]))
  last <- built(plot_community(ta, step = max(steps)))
  expect_equal(nrow(first$data[[1]]), nrow(residents_at(ta, steps[1])))
  expect_equal(first$layout$panel_params[[1]]$x.range,
               last$layout$panel_params[[1]]$x.range)
  expect_error(plot_community(ta, step = 99), "not in tidy")
})

test_that("plot_community follows a log trait scale", {
  ta <- assembled(community_start(bounds(x = c(0.1, 5)), trait_scale = "log",
                                  harness = harness_dd99(x0 = 1, sigma_K = 1,
                                                         sigma_C = 0.4),
                                  fitness_control = list(n_evals = 40)))
  expect_equal(attr(ta, "trait_scale"), "log")
  last <- residents_at(ta, max(as.integer(ta$step)))
  pts <- built(plot_community(ta))$data[[1]]
  expect_equal(sort(pts$x), sort(log10(last$x)))
})

test_that("plot_community places two traits in the trait plane", {
  h <- harness_dd99_nd(x0 = c(0, 0), sigma_K = c(1, 1), sigma_C = c(0.6, 0.6))
  ta <- assembled(community_start(bounds(x1 = c(-3, 3), x2 = c(-3, 3)),
                                  harness = h, trait_scale = "linear",
                                  fitness_control = list(n_evals = 40)))
  expected <- residents_at(ta, max(as.integer(ta$step)))
  p <- plot_community(ta)
  expect_equal(c(p$labels$x, p$labels$y), c("x1", "x2"))
  pts <- built(p)$data[[1]]
  expect_equal(nrow(pts), nrow(expected))
  ord <- order(pts$x, pts$y)
  expect_equal(cbind(pts$x, pts$y)[ord, ],
               cbind(expected$x1, expected$x2)[order(expected$x1,
                                                     expected$x2), ])
  # birth rate is mapped to colour, so distinct rates get distinct colours
  expect_equal(length(unique(pts$colour)), length(unique(expected$births)))
})

test_that("plot_community refuses what it cannot plot", {
  ta <- assembled(community_start(bounds(x = c(-3, 3)), trait_scale = "linear",
                                  harness = harness_dd99(x0 = 0, sigma_K = 1,
                                                         sigma_C = 0.4),
                                  fitness_control = list(n_evals = 40)),
                  nsteps = 2)
  stripped <- ta
  attr(stripped, "trait_scale") <- NULL
  expect_error(plot_community(stripped), "trait_scale")
  expect_s3_class(plot_community(stripped, trait_scale = "linear"), "ggplot")
  expect_error(plot_community(ta[0, ]), "no residents")

  h3 <- harness_dd99_nd(x0 = c(0, 0, 0), sigma_K = c(1, 1, 1),
                        sigma_C = c(0.6, 0.6, 0.6),
                        trait_names = c("x1", "x2", "x3"))
  ta3 <- assembled(community_start(bounds(x1 = c(-3, 3), x2 = c(-3, 3),
                                          x3 = c(-3, 3)),
                                   harness = h3, trait_scale = "linear",
                                   fitness_control = list(n_evals = 20)),
                   nsteps = 2)
  expect_error(plot_community(ta3), "one or two traits")
})
