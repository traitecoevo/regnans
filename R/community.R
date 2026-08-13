##' Initialise a community object
##'
##' Used to store full description of community. Is a list
##' with elements parameters, bounds, birth_rate_initial,
##' trait_names, traits, birth_rate, fitness_control.
##'
##' @title Initialise a community object
##' @param bounds A set of bounds, as specified in \code{plant}.
##' @param birth_rate_initial A vector of birth rates.
##' @param demography_control List controlling the demography/equilibrium
##' solver. See \code{\link{demographic_step_control}}.
##' @param fitness_control List of parameters controlling
##' how invasion-fitness landscapes are generated. See
##' \code{\link{fitness_landscape_control}}. Values given here override the
##' defaults; \code{NULL} uses them all.
##' @param model_support Model support data passed to the harness (e.g. the
##' \code{plant} parameters and control objects).
##' @param harness A model harness wiring the community to a demographic
##' backend (see \code{\link{harness_plant}}); defaults to the \code{plant} SCM.
##' @param trait_scale How trait space is spaced/searched during assembly:
##' \code{"log"} for strictly positive traits or \code{"linear"} for traits
##' spanning zero.
##' @return An \code{community} object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
## TODO: Put birth_rate_initial into parameters and set up
## appropriately?  We use 1 in a couple of places, no?
community_start <- function(bounds,
                            birth_rate_initial = 1e-3,
                            demography_control = demographic_step_control(),
                            fitness_control = NULL,
                            model_support = NULL,
                            harness = NULL,
                            trait_scale = c("log", "linear")
                            ) {

  if (is.character(bounds)) {
    bounds <-  bounds_infinite(bounds)
  }

  ## A harness wires the community to a model backend (see R/harness.R). Default
  ## to the plant SCM, so existing callers that pass only `model_support` keep
  ## working unchanged.
  if (is.null(harness)) {
    harness <- harness_plant()
  }

  ## How trait space is spaced/searched during assembly. "log" suits strictly
  ## positive biological traits (plant: lma, height, seed mass); "linear" suits
  ## traits that span zero or are naturally additive (e.g. DD99/GK98, whose
  ## optima sit at 0). See community_trait_transform().
  trait_scale <- match.arg(trait_scale)

  ## Fitness-landscape options. Previously left NULL, which made
  ## community_fitness_landscape() fail on `fitness_control$method` unless the
  ## caller knew to supply one; now always a complete control list, with any
  ## user-supplied values overriding the defaults (method = "grid").
  fitness_control <- fitness_landscape_control(fitness_control)

  ret <- list(
    bounds = check_bounds(bounds),
    birth_rate_initial = birth_rate_initial,
    demography_control = demography_control,
    fitness_control = fitness_control,
    model_support = model_support,
    harness = harness,
    trait_scale = trait_scale
  )
  ret$trait_names = rownames(ret$bounds)
  ret$traits <- trait_matrix(numeric(0), ret$trait_names)
  ret$birth_rate <- numeric(0)
  
  class(ret) <- "community"
  ret
}

##' Add one or more strategies (residents) to a \code{community}.
##'
##' @title Add strategies to a community
##' @param obj A \code{community} object (see \code{\link{community_start}}).
##' @param traits A trait matrix (see \code{\link{trait_matrix}}) giving the
##'   traits of the strategies to add.
##' @param birth_rate Optional birth rate(s) for the added strategies; if
##'   \code{NULL} a default is used.
##' @return The updated \code{community} object.
##' @export
community_add <- function(obj, traits, birth_rate=NULL) {

  if (is.null(birth_rate)) {
    birth_rate <- obj$birth_rate_initial
  }
  if (!is.matrix(traits)) {
    stop("traits must be a matrix")
  }
  if (length(birth_rate) == 1) {
    birth_rate <- rep_len(birth_rate, nrow(traits))
  }
  if (length(birth_rate) != nrow(traits)) {
    stop("Incompatible length birth rate")
  }
  if (ncol(traits) != ncol(obj$traits)) {
    stop("Incorrect size trait matrix")
  }
  if (!is.null(colnames(traits)) &&
      !identical(colnames(traits), obj$trait_names)) {
    stop("Incorrect traits")
  }
  if (nrow(traits) > 0L) {
    obj$traits <- rbind(obj$traits, traits)
    obj$birth_rate <- c(obj$birth_rate, birth_rate)

    obj <- community_reset(obj)
  }
  obj
}

community_drop <- function(community, which) {
  n_spp <- length(community)
  if (is.logical(which)) {
    if (length(which) != n_spp) {
      stop(sprintf("Invalid length: expected %d, recieved %d",
                   n_spp, length(which)))
    }
    keep <- !which
  } else if (is.numeric(which) || is.integer(which)) {
    which <- as.integer(which)
    if (any(which < 1 | which > n_spp)) {
      stop("Invalid indicies")
    }
    keep <- !(seq_len(n_spp) %in% which)
  } else {
    stop("Invalid index")
  }
  if (!all(keep)) {
    community$traits <- community$traits[keep,,drop=FALSE]
    community$birth_rate <- community$birth_rate[keep]
    community <- community_reset(community)
  }

  if (is.null(community$fitness_function)) {
    community <- community %>% community_demography()
  }

  community
}

community_reset <- function(community) {
  community$resident_fitness <- NULL
  community$fitness_points <- NULL
  community$fitness_function <- NULL
  
  community$model_support$node_schedule_times <- NULL

  community
}


##' Helper function to create trait matrices suitable for
##' \code{\link[plant]{generate_strategy}}.
##'
##' @title Create trait matrix
##' @param x Values
##' @param trait_name Name of a single trait
##' @export
##' @author Rich FitzJohn
trait_matrix <- function(x, trait_name) {
  m <- matrix(x, ncol = length(trait_name))
  colnames(m) <- trait_name
  m
}

##' Returns number of types in community
##'
##' @title Returns number of types in \code{community}
##' @param x A \code{community} object.
##' @return number of rows in trait matrix.
##' @author Rich FitzJohn, Daniel Falster
##' @export
length.community <- function(x) {
  nrow(x$traits)
}
