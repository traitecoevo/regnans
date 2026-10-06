#' Plot an invasion-fitness landscape
#'
#' Plots the sampled fitness landscape held on a community
#' (\code{community$fitness_points}, from
#' \code{\link{community_fitness_landscape}}), with the residents marked and the
#' zero-fitness line that separates trait values that can invade from those that
#' cannot. If the landscape has not been computed yet it is computed here.
#'
#' Where the community carries a Gaussian-process surrogate (the
#' \code{"bayesopt"} landscape method) the surrogate's smooth prediction is
#' drawn over the sampled points; otherwise the sampled points are joined
#' directly. The x axis follows the community's trait scale --- log for the
#' strictly positive plant traits, linear for traits spanning zero.
#'
#' @param community A \code{community} object.
#' @param label Optional label drawn in the top-left corner (e.g. an assembly
#' step number).
#' @param xlim Length-2 trait range for the x axis. Defaults to the community's
#' bounds.
#' @param ylim Optional length-2 fitness range. Applied with
#' \code{coord_cartesian()} so points outside it are clipped rather than
#' dropped; \code{NULL} (the default) lets the data set the range.
#'
#' @return A \code{ggplot} object.
#' @export
community_plot_fitness_landscape <- function(community, label = NA,
                                             xlim = NULL, ylim = NULL) {

  trait_names <- community$trait_names
  if (length(trait_names) != 1L) {
    stop("community_plot_fitness_landscape needs a single-trait community; ",
         "this one has ", length(trait_names), " traits")
  }

  ## Previously this required community$fitness_surrogate_function, which only
  ## exists after the bayesopt landscape method, and read a column literally
  ## named "x" -- so it failed for every grid landscape and for any trait not
  ## called "x". Compute the landscape if needed and plot whatever is there.
  if (is.null(community$fitness_points)) {
    community <- community_fitness_landscape(community)
  }

  ## Build the plotting frame under fixed column names: fitness_points names its
  ## first column after the trait, which is what broke the old hard-coded "x".
  landscape <- tibble::tibble(
    x = as.numeric(community$fitness_points[[trait_names]]),
    fitness = as.numeric(community$fitness_points[["fitness"]])
  )

  if (is.null(xlim)) {
    xlim <- as.numeric(community$bounds[1, ])
  }
  xlim <- range(as.numeric(xlim))

  tf <- community_trait_transform(community)
  log_scale <- identical(tf$scale, "log")

  p <- ggplot2::ggplot(landscape,
                       ggplot2::aes(x = .data[["x"]],
                                    y = .data[["fitness"]])) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed")

  if (is.null(community$fitness_surrogate_function)) {
    ## Grid landscape: join the sampled points.
    p <- p + ggplot2::geom_line(colour = "grey60")
  } else {
    ## Surrogate available: draw its smooth prediction over the samples.
    xs <- tf$inv(seq(tf$fwd(xlim[1]), tf$fwd(xlim[2]), length.out = 500))
    surrogate <- tibble::tibble(
      x = xs,
      fitness = as.numeric(community$fitness_surrogate_function(xs))
    )
    p <- p + ggplot2::geom_line(data = surrogate, colour = "blue")
  }

  p <- p + ggplot2::geom_point(size = 1)

  if (nrow(community$traits) > 0L) {
    residents <- tibble::tibble(
      x = as.numeric(community$traits[, trait_names]),
      fitness = as.numeric(community$fitness_function(
        community$traits[, trait_names]))
    )
    p <- p + ggplot2::geom_point(data = residents, colour = "red", size = 3)
  }

  p <- p +
    ggplot2::xlab(trait_names) +
    ggplot2::ylab("Invasion fitness") +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 16),
                   legend.position = "none")

  p <- p + if (log_scale) {
    ggplot2::scale_x_log10()
  } else {
    ggplot2::scale_x_continuous()
  }

  ## coord_cartesian zooms rather than filters, so the line is not broken by
  ## points falling outside the requested window.
  p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  if (!is.na(label)) {
    p <- p + ggplot2::annotate("text", x = xlim[1], y = Inf, label = label,
                               vjust = "inward", hjust = "inward", size = 5)
  }

  p
}

#' Plot the residents of an assembled community
#'
#' Plots the resident strategies at one step of an assembly, from the output of
#' \code{\link{tidy_assembly}}. With one trait, each resident's birth rate is
#' plotted against its trait value; with two, the residents are placed in the
#' trait plane and coloured by birth rate.
#'
#' Axis limits are taken over every step in \code{tidy}, so plots of successive
#' steps share their axes and can be set side by side or animated:
#' \code{lapply(unique(tidy$step), function(s) plot_community(tidy, step = s))}.
#'
#' @title Plot the residents of an assembled community
#' @param tidy Output of \code{\link{tidy_assembly}}.
#' @param step The assembly step to plot. Defaults to the last.
#' @param trait_scale \code{"log"} or \code{"linear"} trait axes. Defaults to
#'   the scale of the assembled community, which \code{tidy_assembly} records.
#' @param xlim,ylim Optional length-2 axis ranges, replacing those taken from
#'   the whole assembly. Applied with \code{coord_cartesian()}, so points
#'   outside them are clipped rather than dropped.
#'
#' @return A \code{ggplot} object.
#' @export
plot_community <- function(tidy, step = NULL,
                           trait_scale = attr(tidy, "trait_scale"),
                           xlim = NULL, ylim = NULL) {

  if (nrow(tidy) == 0L) {
    stop("tidy has no residents to plot")
  }
  if (is.null(trait_scale)) {
    stop("tidy does not record its trait scale; pass ",
         "trait_scale = \"log\" or \"linear\"")
  }
  trait_scale <- match.arg(trait_scale, c("log", "linear"))

  residents <- tidyr::unnest(tidy[c("step", "births", "traits")], "traits")
  residents$step <- as.integer(residents$step)
  trait_names <- setdiff(names(residents), c("step", "births"))
  if (length(trait_names) > 2L) {
    stop("plot_community plots one or two traits; this assembly has ",
         length(trait_names), ": ", paste(trait_names, collapse = ", "))
  }

  if (is.null(step)) {
    step <- max(residents$step)
  }
  if (!(step %in% residents$step)) {
    stop("step ", step, " is not in tidy, whose steps run ",
         min(residents$step), " to ", max(residents$step))
  }

  ## The y axis is birth rate for one trait and the second trait for two.
  x_name <- trait_names[1]
  y_name <- if (length(trait_names) == 1L) "births" else trait_names[2]
  if (is.null(xlim)) xlim <- range(residents[[x_name]])
  if (is.null(ylim)) ylim <- range(residents[[y_name]])

  log_traits <- trait_scale == "log"

  p <- ggplot2::ggplot(residents[residents$step == step, ],
                       ggplot2::aes(x = .data[[x_name]], y = .data[[y_name]])) +
    (if (log_traits) ggplot2::scale_x_log10() else ggplot2::scale_x_continuous())

  if (length(trait_names) == 1L) {
    p <- p +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_y_log10() +
      ggplot2::ylab("Birth rate")
  } else {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = .data[["births"]]), size = 2) +
      (if (log_traits) ggplot2::scale_y_log10() else ggplot2::scale_y_continuous()) +
      ggplot2::ylab(y_name) +
      ggplot2::scale_colour_viridis_c(transform = "log10", name = "Birth rate")
  }

  p +
    ggplot2::xlab(x_name) +
    ggplot2::theme_classic() +
    ggplot2::theme(text = ggplot2::element_text(size = 16)) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::annotate("text", x = xlim[1], y = Inf,
                      label = paste0("Step = ", step),
                      vjust = "inward", hjust = "inward", size = 5)
}
