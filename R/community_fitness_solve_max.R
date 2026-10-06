##' Fitness of trait value(s) evaluated against a community.
##'
##' Returns the invasion fitness (log net reproduction ratio) of one or more
##' mutant trait values in the community's environment. For an *empty* community
##' this is the fundamental fitness (growth with no competitors).
##'
##' Reimplemented on the community machinery: it evaluates
##' \code{community$fitness_function} (built by
##' \code{plant_community_update_fitness_function}), replacing plant's removed
##' \code{fundamental_fitness()} / \code{fitness_landscape()}.
##'
##' @title Fitness of trait values against a community
##' @param community A \code{community} object.
##' @param values A vector (1D) or matrix (multi-trait) of trait values, as for
##'   \code{\link{trait_matrix}}.
##' @return A numeric vector of fitnesses, one per trait value.
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_growth_rate <- function(community, values) {
  if (is.null(community$fitness_function)) {
    community <- community_update_fitness_function(community)
  }
  community$fitness_function(values)
}

##' Find the trait value of maximum fitness within some bounds.
##'
##' Searches \code{bounds} for the trait value(s) that maximise invasion fitness
##' into the community's environment (the fundamental niche peak for an empty
##' community), on the community's trait scale (see
##' \code{community_start(trait_scale = )}). One trait is searched by golden
##' section over the whole interval; more traits by L-BFGS-B from the centre of
##' the bounds, which finds the maximum nearest it.
##'
##' @title Find point of maximum fitness within some range
##' @param community A \code{community} object.
##' @param bounds Bounds matrix (\code{lower}/\code{upper} per trait). Defaults
##'   to \code{community$bounds}.
##' @param tol Tolerance on the trait scale for the one-trait search; with more
##'   traits L-BFGS-B uses its own convergence test.
##' @return The maximising trait value(s), named by trait, with the achieved
##'   fitness in attribute \code{"fitness"}.
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_fitness <- function(community, bounds = NULL, tol = 1e-3) {
  if (is.null(bounds)) {
    bounds <- community$bounds
  }
  bounds <- check_bounds(bounds)

  if (is.null(community$fitness_function)) {
    community <- community_update_fitness_function(community)
  }

  tf <- community_trait_transform(community)
  lb <- tf$fwd(bounds)
  if (!all(is.finite(lb))) {
    stop("max_fitness needs finite bounds on the trait scale")
  }
  centre <- tf$inv(rowMeans(lb))
  fit <- maximize_scaled(community$fitness_function, centre, bounds, tf,
                         tol = tol)
  structure(as.numeric(fit$par), names = rownames(bounds), fitness = fit$value)
}
