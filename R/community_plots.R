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



#' Returns birth rates for residents in community
#'
#' @param tidy_community Community or history object from tidied assembly
#' @param ... additional arguments passed to \code{plot_community_1d}.
#'
#' @return Returns one plot if community of, if history, plots length of history timeseries
#' @export
#'
plot_community <- function(tidy_community, ...){
  if(is_tibble(tidy_community)){plot_community_1d(tidy_community, ...) -> p}
  if(!is_tibble(tidy_community)){purrr::imap(tidy_community, ~plot_community_1d(tidy_community = .x, step = .y), ...)-> p}
  
  invisible(p)
  
}

plot_community_1d <- function(tidy_community, step = NA, xlim = c(0.01, 1), ylim = c(1e-4, 5)){
  tidy_community %>%
    select(-births, -invader) %>%
    names() -> traits
  tidy_community %>%
    ggplot(aes_string(x = traits[1], y = "births")) + 
    geom_point(aes(colour = invader), size = 2) +
    xlab(traits[1]) +
    ylab("Birth rate") +
    theme_classic() + 
    theme(text = element_text(size = 16),
          legend.position = "none") +
    scale_x_log10(limits = xlim) +
    scale_y_log10(limits = ylim) -> p 
  
  if(!is.na(step)){
    p +
      geom_text(aes(x = xlim[1], y = ylim[2], label = paste0("Step = ",step), vjust = "inward", hjust = "inward"), size = 5) -> p
  }
  
  return(p)
}

#' Plot pair-wise trait combinations of resident community 
#'
#' @param tidy_community Community or history object from tidied assembly with two traits
#' @param ... additional arguments passed to \code{plot_community_2d_internal}.
#'
#' @return Returns one plot if community of, if history, plots length of history timeseries
#' @export
#'
plot_community_2d <- function(tidy_community, ...){
  if(is_tibble(tidy_community)){
    
    ylim_max <- max(tidy_community$hmat)
    ylim_min <- min(tidy_community$hmat)
    ylim <- c(ylim_min, ylim_max)
    
    xlim_max <- max(tidy_community$lma)
    xlim_min <- min(tidy_community$lma)
    xlim <- c(xlim_min, xlim_max)
    
    plot_community_2d_internal(tidy_community, xlim = xlim, ylim = ylim, ...) -> p
  }
  if(!is_tibble(tidy_community)){
    ylim_max <- max(map_dbl(tidy_community, ~max(.x$hmat)))
    ylim_min <- min(map_dbl(tidy_community, ~min(.x$hmat)))
    ylim <- c(ylim_min, ylim_max)
    
    xlim_max <- max(map_dbl(tidy_community, ~max(.x$lma)))
    xlim_min <- min(map_dbl(tidy_community, ~min(.x$lma)))
    xlim <- c(xlim_min, xlim_max)
    
    purrr::imap(tidy_community, ~plot_community_2d_internal(tidy_community = .x, step = .y, xlim = xlim, ylim = ylim), ...)-> p
  }
  
  invisible(p)
  
}

plot_community_2d_internal <- function(tidy_community, step = NA, xlim = c(0.01, 1), ylim = c(1e-4, 5)){
  tidy_community %>%
    select(-births, -invader) %>%
    names() -> traits
  tidy_community %>%
    ggplot(aes_string(x = traits[1], y = traits[2])) + 
    geom_point(aes(colour = births), size = 2) +
    theme_classic() + 
    theme(text = element_text(size = 16),
          legend.position = "none") +
    scale_x_log10(limits = xlim) +
    scale_y_log10(limits = ylim) -> p 
  
  if(!is.na(step)){
    p +
      geom_text(aes(x = xlim[1], y = ylim[2], label = paste0("Step = ",step), vjust = "inward", hjust = "inward"), size = 5) -> p
  }
  
  return(p)
}
