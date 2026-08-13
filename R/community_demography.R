plant_log_eq <- function(...) {
  plant_log_info(..., routine = "equilibrium")
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

demographic_step_control <- function(control=NULL) {

  ## Demography / equilibrium parameters. These live on the community in
  ## community$demography_control. (The plant SCM control() lives separately in
  ## community$model_support$plant_control and is passed to run_scm().)
  ##
  ## Note: plant renamed "seed_rain" -> "birth_rate"/"offspring"; the names here
  ## follow the current plant terminology.
  defaults <- list(
    # which solver community_demography() dispatches to:
    #   "single_step", "equilibrium_iteration",
    #   "equilibrium_solve_nleqslv", "equilibrium_solve_dfsane",
    #   "equilibrium_hybrid"
    equilibrium_solver_name = "equilibrium_iteration",
    equilibrium_eps = 1e-5,
    # iteration solver
    equilibrium_nsteps = 100,
    # runner: reset the integration schedule when birth rates jump more than this
    equilibrium_large_birth_rate_change = 10,
    # birth rate below which a species is treated as extinct
    equilibrium_extinct_birth_rate = 1e-3,
    # root-finding solvers (nleqslv / dfsane)
    equilibrium_min_offspring_arriving = 1e-10,
    equilibrium_solver_logN = TRUE,
    equilibrium_solver_try_keep = TRUE,
    # hybrid solver
    equilibrium_nattempts = 5
  )

  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse=", "))
  }
  ret <- modifyList(defaults, control)

  ret
}

##' Update demography of community according to specified rules
##'
##' @title Update demography of community
##' @param community A \code{community} object.
##' @return  A \code{community} object.
##' @export
##' @author Rich FitzJohn, Daniel Falster
community_demography <- function(community){

  solver <- community$demography_control$equilibrium_solver_name

  plant_log_assembler(sprintf("Updating demography using %s", solver))
  plant_log_assembler_state(community)
  
  if(nrow(community$traits) > 0L) {
    community <- 
      switch(solver,
         single_step = demography_single_step(community),
         equilibrium_iteration = demography_solve_equilibrium_iteration(community),
         equilibrium_hybrid = demography_solve_equilibrium_hybrid(community),
         equilibrium_solve_nleqslv = demography_solve_equilibrium_solve(community, solver = "nleqslv"),
         equilibrium_solve_dfsane = demography_solve_equilibrium_solve(community, solver = "dfsane"),
         stop("Unknown solver ", solver))
    ## community is now the updated community returned by the solver via
    ## community_demography_runner_cleanup() (birth_rate, schedule times and
    ## fitness_points are all set there).
  }

  community_update_fitness_function(community)
}

## This is the simplest update: it simply takes a single step
demography_single_step <- function(community) {

  runner <- community_make_demography_runner(community)
  ## A single step always "converges" (no iteration to assess).
  runner(community$birth_rate)
  community_demography_runner_cleanup(community, runner, converged = TRUE)
}

  
## This is the simplest equilbrium solver: it simply iterates the outgoing offspring
## produced as incoming offspring arrival.  No attempt at projection is made.
demography_solve_equilibrium_iteration <- function(community) {
  
  check <- function(x_in, x_out, eps) {
    achange <- x_out - x_in
    rchange <- 1 - x_out / x_in
    ## eps > 0 && # <- this was a precondition - seems odd.
    converged <- all(abs(achange) < eps | abs(rchange) < eps)
    if (converged) {
      achange <- max(abs(achange))
      rchange <- max(abs(rchange))
      fmt <- "Reached target accuracy (delta %2.5e, %2.5e < %2.5e eps)"
      plant_log_eq(sprintf(fmt, achange, rchange, eps),
                   stage="converged", achange=achange, rchange=rchange)
    }
    converged
  }

  ctrl <- community$demography_control
  birth_rates <- community$birth_rate

  runner <- community_make_demography_runner(community)
  
  for (i in seq_len(ctrl$equilibrium_nsteps)) {
    offspring_production <- runner(birth_rates)
    converged <- check(birth_rates, offspring_production, ctrl$equilibrium_eps)
    birth_rates <- offspring_production
    if (converged) {
      break
    }
  }

  community_demography_runner_cleanup(community, runner, converged)
}

demography_solve_equilibrium_solve <- function(community,
                                         solver="nleqslv") {
  ctrl <- community$demography_control
  
  try_keep <- ctrl$equilibrium_solver_try_keep
  logN <- ctrl$equilibrium_solver_logN
  min_offspring_arriving <- ctrl$equilibrium_min_offspring_arriving

  plant_log_eq(paste("Solving offspring arrival using", solver),
               stage="start", solver=solver)

  birth_rates <- community$birth_rate
  runner <- community_make_demography_runner(community)

  ## First, we exclude species that have offspring arrivals below some minimum
  ## level.
  to_drop <- birth_rates < min_offspring_arriving
  if (any(to_drop)) {
    i_keep <- which(!to_drop)
    msg <- sprintf("Species %s extinct: excluding from search",
                   paste(which(to_drop), collapse=" & "))
    plant_log_eq(msg, stage="drop species", drop=which(to_drop))
    offspring_arriving_full <- birth_rates
    birth_rates <- offspring_arriving_full[i_keep]

    runner_full <- runner
    runner <- function(x) {
      x_full <- rep(0, length(offspring_arriving_full))
      x_full[i_keep] <- x
      runner_full(x_full)[i_keep]
    }
  }

  ## Then see if any species should be retained:
  if (try_keep) {
    ans <- runner(birth_rates)
    keep <- unname(ans >= birth_rates)

    msg <- sprintf("Keeping species %s",
                   paste(which(!to_drop)[keep], collapse=", "))
    plant_log_eq(msg, stage="keep species", keep=which(!to_drop)[keep])
  } else {
    keep <- rep(FALSE, length(birth_rates))
  }

  ## TODO: This is annoying, but turns out to be a problem for getting
  ## the solution working nicely.
  max_offspring_arriving <- pmax(birth_rates * 100, 10000)
  target <- demography_solve_equilibrium_solve_target(runner, keep, logN,
                                               min_offspring_arriving, max_offspring_arriving)
  x0 <- if (logN) log(birth_rates) else birth_rates

  tol <- ctrl$equilibrium_eps
  ## NOTE: Hard coded minimum of 100 steps here.
  maxit <- max(100, ctrl$equilibrium_nsteps)
  sol <- util_nlsolve(x0, target, tol = tol, maxit = maxit, solver = solver)
  
  res <- community_demography_runner_cleanup(community, runner, attr(sol, "converged"))
  
  attr(res, "sol") <- sol
  res
}

## The idea is to use rounds of iteration to try and push the
## system into the basin of attraction of the stable equilibrium.  The
## final approach is slow so use a root-finding approach there.
## However, if we are *not* in the basin of attraction the root finder
## will happily select zero offspring arrivals for species that are not
## zeros.  So after running a round with the solver, check any species
## that were zeroed to make sure they're really dead.
demography_solve_equilibrium_hybrid <- function(community) {

  ctrl <- community$demography_control

  attempts <- ctrl$equilibrium_nattempts

  ## Then expand this out so that we can try alternating solvers
  solver <- rep(c("nleqslv", "dfsane"), length.out=attempts)

  current <- community

  for (i in seq_len(attempts)) {
    eq_solution_iteration <- demography_solve_equilibrium_iteration(current)

    converged_it <- isTRUE(attr(eq_solution_iteration, "converged"))
    msg <- sprintf("Iteration %d %s",
                   i, if (converged_it) "converged" else "did not converge")
    plant_log_eq(msg, step="equilibrium_iteration", converged = converged_it, iteration=i)

    ## Carry the iteration's progress forward so the next round starts from
    ## where this one left off rather than from the original birth rates.
    current <- eq_solution_iteration

    eq_solution <- try(
      demography_solve_equilibrium_solve(
        eq_solution_iteration,
        solver = solver[[i]]
        ),
      silent = TRUE
      )

    converged_sol <- !failed(eq_solution) && isTRUE(attr(eq_solution, "converged"))

    msg <- sprintf("Solve %d %s",
                    i, if (converged_sol) "converged" else "did not converge")
    plant_log_eq(msg, step="solve", converged=converged_sol, iteration=i)

    if (converged_sol) {

      ## Check that species the solver sent to zero really are unviable.
      ##
      ## NOTE: `eq_solution` is a `community`, not a plant `Parameters` object:
      ## the equilibrium state lives in `community$birth_rate`, and the model is
      ## reached through the harness connectors, never through run_scm()
      ## directly. (Previously this block read `eq_solution$strategies` and
      ## called run_scm(), which made the whole solver plant-only and broken
      ## even there.)
      ## NOTE: equilibrium_extinct_birth_rate is an *absolute* birth rate, so
      ## what counts as extinct depends on the model's units -- tune it in
      ## demographic_step_control() when moving between models. It is used the
      ## same way here as in community_check_for_inviable_strategies(): both the
      ## threshold below which a species is treated as gone, and the level at
      ## which it is re-introduced to test that.
      eps <- ctrl$equilibrium_extinct_birth_rate
      extinct <- eq_solution$birth_rate < eps

      if (any(extinct)) {
        plant_log_eq(sprintf("Checking species %s driven to extinction",
                             paste(which(extinct), collapse = ", ")),
                     step = "check extinct", species = which(extinct))

        ## Re-introduce the extinct species at a very low birth rate and take a
        ## single demographic step: if any of them more than replaces itself,
        ## the solver has driven a viable species extinct.
        check <- eq_solution
        check$birth_rate[extinct] <- eps
        runner <- community_make_demography_runner(check)
        offspring_production <- runner(check$birth_rate)

        # `next` breaks the loop iterating over solutions and does not return `eq_solution`
        if (any(offspring_production[extinct] > eps)) {
          plant_log_eq("Solver drove viable species extinct: rejecting")
          next
        }
      }
      plant_log_eq("Accepting solution via solver")
      return(eq_solution)
    }
  }

  ## This one should be a warning?
  plant_log_eq("Repeated rounds failed to find optimum; returning solution from demography_solve_equilibrium_iteration")
  return(eq_solution_iteration)
}

## Another layer of runner for the solver code:
demography_solve_equilibrium_solve_target <- function(runner, keep, logN,
                                               min_offspring_arriving, max_offspring_arriving
                                               ) {
  force(runner)
  force(keep)
  force(logN)
  force(min_offspring_arriving)
  force(max_offspring_arriving)
  function(x, ...) {
    if (logN) {
      x <- exp(x)
    }
    ## TODO: most of the plant_log_eq things here should be DEBUG not INFO?
    ## Avoid negative offspring arrivals:
    x[x < min_offspring_arriving & keep] <- min_offspring_arriving
    x[x < min_offspring_arriving & !keep] <- 0.0
    if (!any(x > 0)) {
      plant_log_eq("All species extinct?")
    }
    too_high <- x > max_offspring_arriving
    if (any(too_high)) {
      plant_log_eq(sprintf("Truncating offspring arrival of species %s",
                           paste(which(too_high), collapse=", ")))
      if (length(max_offspring_arriving) == 1L) {
        max_offspring_arriving <- rep(max_offspring_arriving, length.out=length(x))
      }
      x[too_high] <- max_offspring_arriving[too_high]
    }
    xout <- unname(runner(x))

    xout[keep] <- xout[keep] / x[keep] - 1.0
    i <- !keep & x > 0
    xout[i] <- xout[i] - x[i]

    ## !keep & x <= 0 will now  be zero
    if (any(xout[!keep & x <= 0] != 0.0)) {
      warning("This is not expected", immediate.=TRUE)
    }
    xout
  }
}
