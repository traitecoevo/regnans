## Functions *related* to assembly, but not actually doing it.

#' Find evolutionary attractor for single species and trait.
#'
#' Find evolutionary attractor for single species and trait.
#' This is point at which selection gradient equals zero.
#' Currently solved using \code{uniroot}
#' @param community A \code{community} object to search within.
#' @param bounds a vector containing the end-points of the
#' interval to be searched for the root.
#' @param ... set verbpse=TRUE for verbose output, birth_rate=value
#' gives starting values when solving for demographic equilibrium
#' @param tol the desired accuracy (convergence tolerance).
#' @param edge_ok Is it (not) an error if we end up on the edge of the
#' viable bounds?
#' @author Daniel Falster
#' @export
#' @return A species object, with trait, seed rain and cohort schedule
#' information.
community_solve_singularity_1D <- function(community, bounds = NULL, tol = 1e-04, ...,
                                edge_ok = TRUE) {

  plant_log_assembler(
    sprintf("Solving 1D attractor for %s", community$trait_names))
                                  
  f <- function(x) {
    out <- 
      community %>%
      community_add(trait_matrix(x, community$trait_names)) %>%
      community_demography() %>%
      community_selection_gradient()
    
    ret <- out$selection_gradient
    
    # Add extra details so we can access these later
    attr(ret, "community") <- out

    ret
  }
    
  if(is.null(bounds)) 
    bounds <- community$bounds

  lower <- bounds[[1]]
  upper <- bounds[[2]]

  f_lower <- f(lower)
  f_upper <- f(upper)

  ## This is the exact condition used by uniroot:
  failed <- !isTRUE(as.vector(sign(f_lower) * sign(f_upper) <= 0))

  if (failed) {
    msg <- paste("Bounds do not include attractor: taking",
                 if (f_lower < 0) "lower" else "upper")
    if (edge_ok) {
      warning(msg)
    } else {
      stop(msg)
    }

    if (f_lower < 0) {
      res <- list(root = lower, f.root = f_lower)
    } else {
      res <- list(root = upper, f.root = f_upper)
    }
  } else {
    res <- uniroot(f,
      lower = lower, upper = upper,
      f.lower = f_lower, f.upper = f_upper, tol = tol
      )
  }

  # We're actualy going to 
  community_out <- attr(res$f.root, "community")
  if(community_out$traits != res$root) {
    stop("community_solve_singularity_1D: Hmm, these values should be equla. Better quit now.")
  }

  plant_log_assembler(
    sprintf("Solved! 1D attractor for %s is %s", community$trait_names, as.numeric(res$root))
  )

  community_out
}

#' Selection gradient of the residents
#'
#' The derivative of invasion fitness with respect to the mutant trait,
#' evaluated at each resident's own traits: the direction in which selection
#' pushes that resident. Solve the community to demographic equilibrium first
#' (\code{\link{community_demography}}). Derivative settings come from
#' \code{\link{derivative_control}}.
#' @param community A \code{community} solved to demographic equilibrium.
#' @author Daniel Falster
#' @export
#' @return The community with \code{selection_gradient} (a vector for a single
#' resident, otherwise a residents-by-traits matrix) and \code{resident_fitness}
#' set.
community_selection_gradient <- function(community) {

  msg <- sprintf("Calculating selection gradient for [%s] = [%s]",
    paste(community$trait_names, collapse = ", "), 
    paste(community$traits, collapse = ", ")
  )
  plant_log_assembler(msg)

  g <- community_fitness_gradient(community)
  ret <- if (nrow(g) == 1L) as.vector(g[1, ]) else g

  msg <- sprintf("Solved! Selection gradient for [%s] = [%s] is [%s]", 
    paste(community$trait_names, collapse = ", "), 
    paste(community$traits, collapse = ", "),
    paste(ret, collapse = ", ")
  )
  plant_log_assembler(msg)

  community[["resident_fitness"]] <- attr(g, "value")
  community[["selection_gradient"]] <- ret
  community
}

## NOTE: max_fitness() and max_growth_rate() now live in
## R/community_fitness_solve_max.R, reimplemented on the community machinery.
## The previous plant-style versions here called the removed plant
## fitness_landscape_empty()/fundamental_fitness() and have been deleted.
