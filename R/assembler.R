##' Initialise a community assembly object
##'
##' Control options affect how the assembly proceeds.
##'
##' @title start community assembly process.
##' @param community A \code{community} object.
##' @param control A control object affecting how the assembly proceeds,
##' see \code{assembler_control} for details,
##' @param filename Location to save results.
##' @param prev History from previously run assembly, used to restart.
##' @return An \code{assembler} object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
assembler_start <- function(community, control=assembler_control(control), filename=NULL, prev=NULL) {
  
  ret <- list(history=list(),
              done=FALSE,
              control=control,
              filename=filename)

  if (!is.null(filename) && file.exists(filename)) {
    if (is.null(prev)) {
      prev <- readRDS(filename)
    } else if (!isTRUE(all.equal(readRDS(filename), prev))) {
      stop("filename present and prev given, but don't agree!")
    }
  }

  if (is.null(prev)) {
    ret <- assembler_initialise(ret, community)
  } else {
    ret <- assembler_restore(ret, community, prev)
  }

  class(ret) <- "assembler"
  ret
}

##' Controls how community assembly works.
##'
##' Returns a list. Passing in a list of value via \code{
##' control} will override the defaults. Options include
##' run_type determines whether population is stepped to
##' demographic equilibrium ("to_equilibrium") or not ("single_step").
##' "birth_type" determines sampling of new types -- "stochastic" or
##' "maximum" (on fitness peak). With "stochastic" births,
##' "n_mutants" and "n_immigrants" determine the frequency of
##' resident mutations and immigrations from global pool.
##' "vcv" is variance-covariance matrix for mutations.
##' If "birth_move_tol" is trait distance
##' within which we attempt to move an existing resident rather introduce
##' a new type (this helps reduce the number of types).
##' "compute_viable_fitness" asks whether to check bounds of viable
##' trait space. "check_positive" determines whether the fitness of an
##' invader is checked before it is introduced. If
##' "check_inviable" causes dead residents to be removed when birth rate
##' drops below "dead_birth_rate".
##' "eps_too_close" is tolerance in trait values when searching for maxima.
##'
##' @title Options controllings community assembly process.
##' @param control A list of values to modify from defaults.
##' @return A list with elements run_type, birth_type, birth_move_tol,
##' compute_viable_fitness, n_mutants, n_immigrants, check_positive,
##' vcv, check_inviable, dead_birth_rate, eps_too_close
##' @author Rich FitzJohn, Daniel Falster
##' @export
assembler_control <- function(control=NULL) {
  defaults <- list(run_type="to_equilibrium",
                   birth_type="maximum",
                   birth_move_tol=0,
                   compute_viable_fitness=TRUE,
                   ## births_stochastic_naive:
                   n_mutants=1L,
                   n_immigrants=1L,
                   check_positive=TRUE,
                   vcv=NULL, # will be needed...
                   ## community_deaths:
                   check_inviable=TRUE,
                   dead_birth_rate=1e-3,
                   ## births_maximum_fitness:
                   eps_too_close=1e-3,
                   min_birth_rate_initial = 1e-3,
                   max_birth_rate_initial = 500,
                  ## This magic number is not great, but needed to 
                  ## prevent suggesting adding 1e8 as the birth 
                  ## rate (can happen!).  This should be "quite 
                  ## large", but obviously that number depends on
                  ## the situation.
                  eps_fitness_invasion = 0.01
                  )
  if (identical(control[["birth_type"]], "stochastic")) {
    defaults$run_type <- "single_step"
  }

  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse=", "))
  }
  ret <- modifyList(defaults, control)
  if (ret$birth_type == "maximum" && ret$run_type != "to_equilibrium") {
    stop("Must use 'to_equilibrium' run type with maximum fitness births")
  }
  if (ret$birth_type == "stochastic" && is.null(ret$vcv)) {
    stop("vcv must be provided")
  }

  ret
}

assembler_initialise <- function(obj, community) {
  if (!inherits(community, "community")) {
    stop("Expecting a community object")
  }
  plant_log_assembler("Starting empty assembler")
  obj$community <- community %>% community_demography()
  
  if (isTRUE(obj$control$compute_viable_fitness)) {
    plant_log_assembler("Computing viable bounds")
    obj$community <- community_viable_bounds(obj$community)
  }
  obj$done <- is.null(obj$community$bounds)
  obj <- assembler_append_history(obj)
  obj
}


assembler_restore <- function(obj, community, prev) {
  plant_log_assembler("Restoring previous community and history")

  if (length(community) > 0L) {
    stop("community must be an empty placeholder")
  }
  ## Check that this looks legit:
  x <- last(prev$history)
  if (!isTRUE(all.equal(community$parameters, x$parameters))) {
    stop("Restoring community: parameters don't match")
  }
  if (!all(x$bounds[, "lower"] >= community$bounds[, "lower"] &
           x$bounds[, "upper"] <= community$bounds[, "upper"])) {
    stop("Restoring community: incompatible bounds")
  }
  plant_log_assembler(sprintf("...restored %d steps of history",
                              length(prev$history)))
  prev
}

assembler_append_history <- function(obj) {

  ## The fitness-landscape grid is 1D (it feeds find_max_fitness_1D and 1D
  ## plotting). For multi-trait communities the nD max search
  ## (find_max_fitness_2d) does not need it, so skip it.
  if (obj$control$birth_type == "maximum" &&
      length(obj$community$trait_names) == 1L) {
    obj$community <- community_fitness_landscape(obj$community)
  }

  obj$history <- c(obj$history, list(obj$community))
  if (isTRUE(obj$done)) {
    attr(obj$history, "done") <- TRUE
  }
  filename <- obj$filename
  if (!is.null(filename)) {
    ok <- try(saveRDS(obj, filename))
    if (inherits(ok, "try-error")) {
      warning("History saving has failed",
              immediate.=TRUE, call.=FALSE)
    }
  }
  obj
}

assembler_step <- function(obj) {
  plant_log_assembler(sprintf("step %d, (%d strategies), %s",
                  length(obj$history), length(obj$community),
                  obj$control$run_type))
  plant_log_assembler_state(obj$community)
  if (!obj$done) obj <- assembler_births(obj)
  if (!obj$done) obj$community <- community_demography(obj$community)
  if (!obj$done) obj <- assembler_deaths(obj)
  assembler_append_history(obj)
}

##' Set the initial traits (and optional birth rates) on an empty assembler.
##'
##' @title Set traits on an assembler
##' @param obj An \code{assembler} object with an empty community.
##' @param traits A trait matrix (see \code{\link{trait_matrix}}).
##' @param birth_rate Optional birth rate(s); if \code{NULL} a default is used.
##' @return The updated \code{assembler} object.
##' @export
assembler_set_traits <- function(obj, traits, birth_rate=NULL) {
  if (length(obj$community) > 0L) {
    stop("This is for an empty community only")
  }
  plant_log_assembler("Setting traits")
  obj$community <- community_add(obj$community, traits, birth_rate)
  plant_log_assembler_state(obj$community)
  obj$community <- community_demography(obj$community)
  obj <- assembler_deaths(obj)
  obj <- assembler_append_history(obj)
  obj
}

##' Take specified number of steps with assembler.
##'
##'
##' @title Take specified number of steps with assembler.
##' @param obj An  \code{assembler} object.
##' @param nsteps Number of steps
##' @return An \code{assembler} object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
assembler_run <- function(obj, nsteps) {
  for (i in seq_len(nsteps)) {
    if (obj$done) {
      plant_log_assembler(sprintf("Assembly completed after %d steps",
                                  length(obj$history)))
      break
    }
    obj <- assembler_step(obj)
  }
  obj
}

##' Returns number taken by an assembler object
##'
##'
##' @title Returns number taken by an \code{assembler} object
##' @param x An  \code{assembler} object.
##' @return Length of history in assembler object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
length.assembler <- function(x) {
  length(x$history)
}

should_move <- function(obj, to_add) {
  ret <- integer(0)
  if (nrow(to_add) == 1L &&
      length(obj$community) > 0L &&
      obj$control$run_type == "to_equilibrium") {
    ## NOTE: not rescaled by obj$bounds?
    tf <- community_trait_transform(obj$community)
    i <- closest(tf$fwd(to_add), tf$fwd(obj$community$traits))
    if (attr(i, "distance") < obj$control$birth_move_tol) {
      ret <- i
    }
  }
  ret
}

## Helper function for printing community state.
plant_log_assembler_state <- function(community) {
  str <- format_community_state(community$traits,
                                community$birth_rate,
                                sprintf("\t%s\t", c("res", seq_along(community$birth_rate))))
  msg <- paste(c("Traits:", str), collapse="\n")
  plant_log_state(msg,
                      traits=community$traits,
                      birth_rate=community$birth_rate)
}

format_community_state <- function(traits, birth_rate, prefix=NULL) {
  m <- cbind(as.data.frame(traits, stringsAsFactors=FALSE),
             birth_rate=I(prettyNum(birth_rate)))
  if (nrow(m) == 0) {
    m[1,] <- rep("<empty>", ncol(m))
    rownames(m) <- ""
  }
  str <- capture.output(print(m))
  if (!is.null(prefix)) {
    str <- paste0(prefix, str)
  }
  str
}
