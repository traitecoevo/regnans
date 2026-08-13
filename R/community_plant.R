##' Default \code{plant} FF16 parameters tuned for community-assembly runs.
##'
##' @title Default plant parameters for assembly
##' @param hmat Height at maturation.
##' @param max_patch_lifetime Maximum patch lifetime.
##' @param fixed_RA If \code{TRUE}, use a fixed reproductive allocation.
##' @return A \code{plant} \code{Parameters} object.
##' @export
plant_default_assembly_pars <- function(hmat = 10, max_patch_lifetime = 60, fixed_RA = FALSE) {
  p <- scm_base_parameters("FF16")
  ## Biological parameters live under `$pars` (plant #410). Assigning the flat
  ## name instead attaches an inert list element to the R-side object without
  ## erroring, so the value never reaches the strategy -- see #45.
  p$strategy_default$pars$a_l1 <- 2.17
  p$strategy_default$pars$a_l2 <- 0.5
  p$strategy_default$pars$hmat <- hmat
  p$max_patch_lifetime <- max_patch_lifetime
  ## Regenerate the default node introduction times to match the (possibly
  ## lowered) patch lifetime. scm_base_parameters() builds these for its own
  ## default max_patch_lifetime (~105), so without this the default schedule
  ## spans times beyond max_patch_lifetime and the SCM errors with
  ## "time_max must be greater than (or equal to) current time" whenever the
  ## equilibrium runner resets the schedule to its defaults.
  p$node_schedule_times_default <-
    plant::node_schedule_times_default(max_patch_lifetime)
  if(fixed_RA) {
    # Fixed allocation to reproduction
    p$strategy_default$pars$a_f1 <- 0.5
    p$strategy_default$pars$a_f2 <- 0
  }
  p
}

##' Default plant SCM control for community assembly.
##'
##' Thin wrapper around \code{plant::control()}. Kept as a named helper so the
##' plant control used for assembly lives in one place and can be tuned in
##' future (e.g. schedule tolerances) without touching every call site.
##'
##' @title Default plant control for community assembly
##' @param ... Passed to \code{plant::control()} to override defaults.
##' @return A plant \code{Control} object.
##' @export
##' @author Daniel Falster
plant_default_assembly_control <- function(...) {
  plant::control(...)
}

plant_community_hyperpar <- function(community) {
 if (is.null(community$model_support$hyperpar)) {
   hyperpar <- plant::param_hyperpar(community$model_support$p)
 } else {
   hyperpar <- community$model_support$hyperpar
 }
 hyperpar
}


plant_community_parameters <- function(community) {
  
  p <- community$model_support$p

  hyperpar <- plant_community_hyperpar(community)

  if(nrow(community$traits) > 0)
    p$strategies <- plant::generate_strategy(p, community$traits, hyperpar = hyperpar, birth_rate = community$birth_rate)

  if (!is.null(community$model_support$node_schedule_times)) {
    p$node_schedule_times <- community$model_support$node_schedule_times
  }
  if (!is.null(community$model_support$node_schedule_ode_times)) {
    p$ode_times <- community$model_support$node_schedule_ode_times
  }

  p
}


## Support code:
plant_community_make_demography_runner <- function(community) {
  pretty_num_collapse <- function(x, collapse = ", ") {
    paste0("{", paste(prettyNum(x), collapse = collapse), "}")
  }

  p <- plant_community_parameters(community)
  ctrl <- community$model_support$plant_control  # plant SCM control() for run_scm()

  # default is about 10 ind.m-2
  large_offspring_arriving_change <-
    community$demography_control$equilibrium_large_birth_rate_change

  # traverse list of strategies and pull birth_rates
  last_arrival_rates <- purrr::map_dbl(
    p$strategies,
    ~ purrr::pluck(., "birth_rate_y")
  )

  default_schedule_times <- rep(
    list(p$node_schedule_times_default),
    length(last_arrival_rates)
  )

  last_schedule_times <- p$node_schedule_times
  last_offspring_production <- NULL
  history <- NULL
  counter <- 1L

  function(birth_rates) {
    # if a runner has diverged significantly in the last iteration then
    # reset the schedule to it's defaults and rebuild
    if (any(abs(birth_rates - last_arrival_rates) > large_offspring_arriving_change)) {
      p$node_schedule_times <- default_schedule_times
    }

    # set birth rates
    for (i in seq_along(p$strategies)) {
      p$strategies[[i]]$birth_rate_y <- birth_rates[i]
    }

    # update schedule - starts from prev. schedule so should be fast for fine scale resolution
    # plant migration: build_schedule() removed; run_scm(refine_schedule = TRUE)
    # returns the SCM object carrying both the refined parameters and
    # offspring_production (previously a build_schedule() attribute).
    scm <- run_scm(p, ctrl = ctrl, refine_schedule = TRUE)
    p_new <- scm$parameters
    offspring_production <- scm$offspring_production


    # (ANDREW) Double arrow modifies counter in parent environment..
    # I don't love this approach to counting, as it introduces a side effect
    # to an otherwise functional programming design

    ## These all write up to the containing environment:
    p <<- p_new
    last_schedule_times <<- p_new$node_schedule_times
    last_offspring_production <<- offspring_production
    history <<- c(history, list(c("in" = birth_rates, "out" = offspring_production)))

    msg <- sprintf(
      "eq> %d: %s -> %s (delta = %s)", counter,
      pretty_num_collapse(birth_rates),
      pretty_num_collapse(offspring_production),
      pretty_num_collapse(offspring_production - birth_rates)
    )
    plant_log_eq(msg,
      stage = "runner",
      iteration = counter,
      birth_rate = birth_rates,
      offspring_production = offspring_production
    )

    counter <<- counter + 1L

    # TODO: check why returning node schedule as an attribute
    # attr(offspring_production, "schedule_times") <- last_schedule_times

    offspring_production
  }
}

# TODO: revisit 'gross' behaviour in cleanup utility
plant_community_demography_runner_cleanup <- function(community, runner, converged = TRUE) {
  # this pulls the history of the runner from the parent environment
  e <- environment(runner)
  if (is.function(e$runner_full)) {
    runner <- e$runner_full
    e <- environment(runner)
  }

  # the final solution has a recently built integration schedule
  p <- e$p
  p$node_schedule_times <- e$last_schedule_times

  ## The equilibrium birth rate is the offspring production from the final
  ## runner evaluation. Previously this travelled as attr(p, "offspring_production")
  ## (a build_schedule() side-channel that plant has since removed); the runner
  ## now stashes it explicitly in last_offspring_production.
  community$birth_rate <- e$last_offspring_production
  community$model_support$node_schedule_times <- p$node_schedule_times
  community$model_support$node_schedule_ode_times <- p$ode_times
  community$fitness_points <- NULL

  attr(community, "progress") <- util_rbind_list(e$history)
  attr(community, "converged") <- converged
  community
}


plant_community_viable_bounds <- function(community) {
  if (length(community) > 0) {
    stop("You don't want to run this on an existing community")
  }
  community$bounds <- community_viable_fitness_1D(community)
  community
}

plant_community_update_fitness_function <- function(community) {
  
  p <- plant_community_parameters(community)
  hyperpar <- plant_community_hyperpar(community)

  ctrl <- community$model_support$plant_control
  ctrl$save_RK45_cache <- T
  
  if (length(p$strategies) > 0L) {
    # if there's a resident, use the saved environment to calculate mutant fitness
    scm <- run_scm(p, ctrl = ctrl,
      use_ode_times = length(p$ode_times) > 0)
    community$resident_fitness <- log(scm$net_reproduction_ratios)
    } else {
    community$resident_fitness <- numeric()
  }

  community$fitness_function <- 
    function(x) {
      traits <- trait_matrix(x, community$trait_names)

      p_mutants <- plant::add_mutant(p, traits, hyperpar,
        birth_rate = rep(0, nrow(traits)))
      
      if (length(p$strategies) > 0L) {
        scm$run_mutant(p_mutants)
        ret <- log(scm$net_reproduction_ratios)
      } else {
        # otherwise just run mutants with zero birth rate
        scm <- run_scm(p_mutants, ctrl = ctrl, use_ode_times = 0)
        ret <- log(scm$net_reproduction_ratios)
      }
      
    ret
  }

  community
}

##' Check low-abundance strategies for viability.
##'
##' @title Check low-abundance strategies for viability
##' @param community A community object
##' @export
plant_community_check_for_inviable_strategies <- function(community) {
  # TODO: Remove direct depoency o n plant, rewrite to use existing community fitness functions

  p <- community_parameters(community)

  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a offspring arrival of less than
  ## `eps_test * max(p$birth_rate)` will be tested.  By default
  ##  this is 1 100th of the maximum offspring arrival.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- community$demography_control$equilibrium_extinct_birth_rate
  ## TODO: This was ctrl$equilibrium_inviable_test, but I think
  ## that birth offspring arrival actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  birth_rate <- sapply(p$strategies, function(s) s$birth_rate_y, simplify = TRUE)
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing offspring
  ## rain.
  runner <- community_make_demography_runner(community)
  offspring_production <- runner(birth_rate)

  test <- which(offspring_production < birth_rate &
    birth_rate < max(offspring_production) * eps_test)
  test <- test[order(offspring_production[test])]

  drop <- logical(length(offspring_production))

  for (i in test) {
    plant_log_inviable(paste("Testing species", i),
      stage = "testing", species = i
    )
    x <- offspring_production
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      plant_log_inviable(paste("Removing species", i),
        stage = "removing", species = i
      )
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      offspring_production <- res
    }
  }

  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | offspring_production < eps

  attr(offspring_production, "drop") <- drop
  offspring_production
}

# The model-agnostic connectors (community_parameters, community_*_runner*,
# community_viable_bounds, community_check_for_inviable_strategies,
# community_update_fitness_function) live in R/harness.R and forward to the
# harness on the community. harness_plant() (also in R/harness.R) points them
# back at the plant_community_* implementations above.
