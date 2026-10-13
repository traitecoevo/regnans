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
##' With \code{evolve} other than \code{"none"} the residents evolve between
##' births: after each step's deaths (and once at the start) every resident is
##' moved to the singular coalition they form together, so that each
##' maximum-fitness birth is an immigrant into an evolutionarily stationary
##' community and the assembly ends at an evolutionarily stable community
##' (\code{\link{community_label}}). \code{"newton"} solves for the coalition
##' with \code{\link{community_solve_singularity}} (to \code{evolve_tol}) and
##' classifies it with \code{\link{community_classify_singularity}}; where
##' that solve fails, or finds a root the dynamics does not reach (not
##' convergence stable, or without protected coexistence), the dynamics are
##' integrated instead from where the residents were. \code{"canonical"}
##' always integrates: \code{\link{community_canonical_equation}} to its
##' first stationary coalition (\code{branch = "none"}, polished to
##' \code{evolve_tol}), which only reaches attractors and drops residents
##' that go extinct or merge, at the cost of an ODE integration. Either
##' way an evolution that fails ends the assembly, reported in
##' \code{$evolved}. A maximum of fitness within
##' \code{eps_too_close} of a resident is then that resident at a branching
##' point, and the birth is placed \code{branch_offset} (a fraction of the
##' bounds' width, on the trait scale) away from it towards the maximum rather
##' than ending the assembly. A birth is refused, ending the assembly with
##' \code{$refused} saying why, when it would take the community beyond
##' \code{max_residents}, or (with \code{edge_tol > 0} and a resident
##' present) when it lies within \code{edge_tol} of the bounds' edge, that
##' fraction of their width on the trait scale: the box, not the model, is then
##' setting the community. A trait whose bounds are not finite on the trait
##' scale has no edge.
##'
##' @title Options controllings community assembly process.
##' @param control A list of values to modify from defaults.
##' @return A list with elements run_type, birth_type, birth_move_tol,
##' compute_viable_fitness, n_mutants, n_immigrants, check_positive,
##' vcv, check_inviable, dead_birth_rate, eps_too_close, evolve, evolve_tol,
##' branch_offset, max_residents, edge_tol
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
                  eps_fitness_invasion = 0.01,
                  ## evolution of the residents between births
                  evolve = "none",
                  evolve_tol = 1e-6,
                  branch_offset = 1e-2,
                  ## births refused: beyond this many residents, or on the
                  ## edge of the bounds (0 turns the edge rule off)
                  max_residents = Inf,
                  edge_tol = 0
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
  ret$evolve <- match.arg(ret$evolve, c("none", "newton", "canonical"))
  if (ret$evolve != "none" && ret$birth_type != "maximum") {
    stop("evolve needs birth_type = \"maximum\": the residents evolve between ",
         "maximum-fitness births")
  }
  for (nm in c("evolve_tol", "branch_offset")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a single positive number")
    }
  }
  if (!is.numeric(ret$max_residents) || length(ret$max_residents) != 1L ||
      is.na(ret$max_residents) || ret$max_residents < 1) {
    stop("max_residents must be a number of residents of at least 1")
  }
  if (!is.numeric(ret$edge_tol) || length(ret$edge_tol) != 1L || !is.finite(ret$edge_tol) || ret$edge_tol < 0) {
    stop("edge_tol must be a single non-negative number")
  }

  ret
}

assembler_initialise <- function(obj, community) {
  if (!inherits(community, "community")) {
    stop("Expecting a community object")
  }
  plant_log_assembler("Starting empty assembler")
  obj$community <- community %>% community_demography()
  obj$evaluations <- 1L

  if (isTRUE(obj$control$compute_viable_fitness)) {
    plant_log_assembler("Computing viable bounds")
    obj$community <- community_viable_bounds(obj$community)
  }
  obj$done <- is.null(obj$community$bounds)
  if (!obj$done) obj <- assembler_evolve(obj)
  obj <- assembler_append_history(obj)
  obj
}

## Move the residents to the singular coalition they form together, for an
## assembly with evolve other than "none". Two residents the evolution brings
## within two finite-difference steps of each other have merged, as at a
## pitchfork from the side where the pair no longer exists: the rarer is
## dropped and the survivors evolve again, as the canonical equation's merge
## does. A coalition the evolution cannot reach ends the assembly: $evolved
## says why.
assembler_evolve <- function(obj) {
  ctrl <- obj$control
  if (ctrl$evolve == "none" || length(obj$community) == 0L) {
    return(obj)
  }
  route <- switch(ctrl$evolve, newton = evolve_newton, canonical = evolve_canonical)
  community <- obj$community
  for (attempt in seq_len(length(community))) {
    plant_log_assembler(sprintf("Evolving %d residents to their singular coalition (%s)",
                                length(community), ctrl$evolve))
    res <- route(community, ctrl)
    obj$evaluations <- obj$evaluations + res$evaluations
    merged <- if (res$converged) evolve_merged(res$community) else integer(0)
    if (length(merged) == 0L) {
      break
    }
    plant_log_assembler(sprintf("Residents %s merged; dropping the rarer", paste(merged, collapse = ", ")))
    community <- community_drop(res$community, merged)
    obj$evaluations <- obj$evaluations + 1L
  }
  obj$community <- res$community
  obj$evolved <- res[c("converged", "classification", "reason")]
  if (!res$converged) {
    obj$done <- TRUE
  }
  obj
}

## The rarer resident of each pair within two finite-difference steps of
## each other, the converse of a branching.
evolve_merged <- function(community) {
  m <- nrow(community$traits)
  pairs <- residents_merged(as.numeric(community$traits), m, length(community$trait_names),
                            community_derivative_control(community), reach = 2)
  n <- as.numeric(community$birth_rate)
  sort(unique(apply(pairs, 1, function(p) p[which.min(n[p])])))
}

## Solve for the coalition and classify it. Newton from the residents is
## fast but local: from an immigrant's starting point it can leave the bounds,
## collapse two residents onto a single resident's root, or land on a root the
## dynamics flees (a repeller, a Garden of Eden, an unprotected coalition).
## In all of these the dynamics themselves are integrated instead, from where
## the residents were: the assembly's next birth is only meaningful into a
## community the dynamics reaches.
evolve_newton <- function(community, ctrl) {
  quietly <- function(expr) suppressWarnings(expr)
  sol <- tryCatch(quietly(community_solve_singularity(community, tol = ctrl$evolve_tol)),
                  error = function(e) e)
  evaluations <- if (inherits(sol, "error")) {
    if (is.null(sol$evaluations)) 0L else sol$evaluations
  } else {
    attr(sol, "evaluations")
  }
  why <- if (inherits(sol, "error")) {
    conditionMessage(sol)
  } else if (!isTRUE(attr(sol, "converged"))) {
    "the coalition solve did not converge"
  } else if (length(evolve_merged(sol)) > 0L) {
    ## a pair collapsed onto one point passes as a root and cannot be
    ## classified; the caller merges it
    return(list(community = sol, converged = TRUE, classification = NA_character_,
                reason = NULL, evaluations = evaluations))
  } else {
    cl <- quietly(community_classify_singularity(sol))
    evaluations <- evaluations + cl$evaluations
    if (cl$classification %in% c("CSS", "branching point")) {
      return(list(community = sol, converged = TRUE, classification = cl$classification,
                  reason = NULL, evaluations = evaluations))
    }
    sprintf("the singular coalition found is %s", cl$classification)
  }
  plant_log_assembler(paste0(why, "; integrating the canonical equation"))
  out <- evolve_canonical(community, ctrl)
  out$evaluations <- out$evaluations + evaluations
  out
}

## Integrate the canonical equation to its first stationary coalition: the
## route that only reaches attractors.
evolve_canonical <- function(community, ctrl) {
  ## the polish solves the gradients to a tenth of gradient_tol
  ce <- tryCatch(community_canonical_equation(community, x0 = community$traits, control = canonical_control(list(
    branch = "none", polish = TRUE, gradient_tol = 10 * ctrl$evolve_tol,
    max_residents = length(community)))), error = function(e) e)
  if (inherits(ce, "error")) {
    return(list(community = community, converged = FALSE, classification = NA_character_,
                reason = paste("the canonical equation failed:", conditionMessage(ce)),
                evaluations = 0L))
  }
  stable <- identical(ce$outcome, "stable")
  ## the community at the end of the trajectory, solved again at its
  ## polished traits (the run's last right-hand side predates the polish)
  final <- canonical_community(ce, max(ce$trajectory$time))
  list(community = final, converged = stable, classification = NA_character_,
       reason = if (stable) NULL else sprintf("the canonical equation ended with outcome %s", ce$outcome),
       evaluations = as.integer(ce$evaluations) + 1L)
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
  ## an assembler saved before a control field or the count of solves existed
  ## takes their defaults
  prev$control <- assembler_control(prev$control)
  if (is.null(prev$evaluations)) prev$evaluations <- 0L
  prev
}

assembler_append_history <- function(obj) {

  ## The fitness-landscape grid is 1D: it gives find_max_fitness() its starting
  ## bracket in one trait, and the 1D plots. With more traits the maximum is
  ## found by multistart from the residents instead, so skip it.
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
  if (!obj$done) {
    obj$community <- community_demography(obj$community)
    obj$evaluations <- obj$evaluations + 1L
  }
  if (!obj$done) obj <- assembler_deaths(obj)
  if (!obj$done) obj <- assembler_evolve(obj)
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
  obj$evaluations <- obj$evaluations + 1L
  obj <- assembler_deaths(obj)
  obj <- assembler_evolve(obj)
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
