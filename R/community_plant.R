##' Default \code{plant} parameters tuned for community-assembly runs.
##'
##' For FF16 the leaf-area allometry is set to the assembly values
##' (\code{a_l1 = 2.17}, \code{a_l2 = 0.5}); TF24 keeps its own. Pass the same
##' \code{model} to \code{\link{harness_plant}}.
##'
##' @title Default plant parameters for assembly
##' @param hmat Height at maturation.
##' @param max_patch_lifetime Maximum patch lifetime.
##' @param fixed_RA If \code{TRUE}, use a fixed reproductive allocation.
##' @param model plant physiological model: \code{"FF16"} or \code{"TF24"}.
##' @return A \code{plant} \code{Parameters} object.
##' @export
plant_default_assembly_pars <- function(hmat = 10, max_patch_lifetime = 60, fixed_RA = FALSE,
                                        model = c("FF16", "TF24")) {
  model <- match.arg(model)
  ## plant's testing set: today's defaults with the frozen trade-off functions
  ## attached as p$hyperpar (plant#731). The bare base_parameters() carries
  ## none, so a trait written through it would reach pars without its trade-offs
  ## (lma without k_l, say) and the assembly would run a different model.
  p <- default_parameters(model, quiet = TRUE)
  ## Biological parameters live under `$pars` (plant #410). Assigning the flat
  ## name instead attaches an inert list element to the R-side object without
  ## erroring, so the value never reaches the strategy -- see #45.
  if (model == "FF16") {
    p$strategy_default$pars$a_l1 <- 2.17
    p$strategy_default$pars$a_l2 <- 0.5
  }
  p$strategy_default$pars$hmat <- hmat
  p$max_patch_lifetime <- max_patch_lifetime
  ## Regenerate the default node introduction times to match the (possibly
  ## lowered) patch lifetime. base_parameters() builds these for its own
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
##' The size density is carried on the birth-date coordinate
##' (\code{node_density_coordinate = "birth_date"}) unless overridden: plant's
##' gradients (\code{plant::stand_gradient()}) are defined on it alone.
##'
##' @title Default plant control for community assembly
##' @param ... Passed to \code{plant::control()} to override defaults.
##' @return A plant \code{Control} object.
##' @export
##' @author Daniel Falster
plant_default_assembly_control <- function(...) {
  args <- list(...)
  if (is.null(args$node_density_coordinate)) {
    args$node_density_coordinate <- "birth_date"
  }
  do.call(plant::control, args)
}

## The plant model a Parameters object was built for ("FF16", "TF24", ...),
## read off its class, or NULL where p is not plant Parameters.
plant_parameters_model <- function(p) {
  pattern <- "^Parameters<([^,>]+),.*$"
  if (is.null(p) || !grepl(pattern, class(p)[[1]])) {
    return(NULL)
  }
  sub(pattern, "\\1", class(p)[[1]])
}

## The plant Environment a community's runs take, built afresh from
## `model_support$environment`, or NULL (plant's default for the model) where
## there is none. The specification is plain data, so it crosses to a worker
## where an Environment, an external pointer, cannot; and each run gets its own,
## so none starts from another's soil state. `soil` is passed to
## plant::set_tf24_soil() (TF24 only); each entry of `drivers` is a constant,
## or a list(x, y) of control points for a variable driver.
plant_community_environment <- function(community) {
  spec <- community$model_support$environment
  if (is.null(spec)) {
    return(NULL)
  }
  unknown <- setdiff(names(spec), c("soil", "drivers"))
  if (length(unknown) > 0L) {
    stop("model_support$environment takes soil and drivers; unknown: ", paste(unknown, collapse = ", "))
  }
  model <- plant_parameters_model(community$model_support$p)
  if (is.null(model)) {
    stop("model_support$environment needs model_support$p, plant Parameters, to say which model it is for")
  }
  env <- plant::Environment(model)
  if (!is.null(spec$soil)) {
    if (!startsWith(model, "TF24")) {
      stop("model_support$environment$soil is a TF24 soil profile; ", model, " has none")
    }
    do.call(plant::set_tf24_soil, c(list(env), spec$soil))
  }
  drivers <- as.list(spec$drivers)
  if (length(drivers) > 0L) {
    if (is.null(names(drivers)) || any(names(drivers) == "")) {
      stop("model_support$environment$drivers must be named, one entry per driver")
    }
    if (!is.function(env$extrinsic_drivers_get_names)) {
      stop("model_support$environment$drivers: ", model, " has no extrinsic drivers")
    }
    ## plant accepts any name, making a driver nothing reads, so a misspelt
    ## one would run the default silently
    known <- env$extrinsic_drivers_get_names()
    unknown <- setdiff(names(drivers), known)
    if (length(unknown) > 0L) {
      stop("Unknown driver for ", model, ": ", paste(unknown, collapse = ", "),
           "; its drivers are ", paste(known, collapse = ", "))
    }
  }
  for (name in names(drivers)) {
    driver <- drivers[[name]]
    if (is.list(driver)) {
      if (!is.numeric(driver$x) || !is.numeric(driver$y) || length(driver$x) != length(driver$y) ||
          length(driver$x) < 2L) {
        stop("Driver ", name, ": a variable driver is list(x, y), numeric and of one length, at least 2")
      }
      env$extrinsic_drivers_set_variable(name, driver$x, driver$y)
    } else {
      if (!is.numeric(driver) || length(driver) != 1L || !is.finite(driver)) {
        stop("Driver ", name, ": a constant driver is one finite number, a variable one list(x, y)")
      }
      env$extrinsic_drivers_set_constant(name, driver)
    }
  }
  env
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
  # times and step sizes together replay the equilibrium run exactly; times
  # alone would let the controller choose its own sub-steps. The soil's
  # substep records go with the sizes: plant refuses sizes without them.
  ms <- community$model_support
  if (!is.null(ms$node_schedule_ode_times)) {
    p$ode_times <- ms$node_schedule_ode_times
    if (!is.null(ms$node_schedule_ode_step_sizes)) {
      p$ode_step_sizes <- ms$node_schedule_ode_step_sizes
      p$ode_soil_input_slopes <- ms$node_schedule_ode_soil_input_slopes
      p$ode_soil_substep_ends <- ms$node_schedule_ode_soil_substep_ends
    }
  }
  # a frozen discretisation (a stencil point's) overrides both, and survives
  # community_add(); it describes as many residents as the centre it was read
  # from had
  frozen <- community$model_support$frozen
  if (!is.null(frozen) && length(p$strategies) > 0L) {
    if (length(frozen$node_schedule_times) != length(p$strategies)) {
      stop("A frozen discretisation for ", length(frozen$node_schedule_times),
           " residents cannot be replayed with ", length(p$strategies))
    }
    p$node_schedule_times <- frozen$node_schedule_times
    p$ode_times <- frozen$ode_times
    p$ode_step_sizes <- frozen$ode_step_sizes
    p$ode_soil_input_slopes <- frozen$ode_soil_input_slopes
    p$ode_soil_substep_ends <- frozen$ode_soil_substep_ends
  }

  p
}

## Drop the ODE schedule p carries so a run integrates freely. The four fields
## go together: plant refuses step sizes without the soil's substep records,
## and a refining run records all four into its parameters.
plant_clear_ode_schedule <- function(p) {
  p$ode_times <- numeric(0)
  p$ode_step_sizes <- numeric(0)
  p$ode_soil_input_slopes <- list()
  p$ode_soil_substep_ends <- list()
  p
}

## The discretisation a solved plant community was integrated on -- its node
## schedule and the ODE steps of its equilibrium run, and the patch lifetime
## they span -- as plain numbers, or NULL where it records none.
plant_community_discretisation <- function(community) {
  ms <- community$model_support
  if (length(ms$node_schedule_ode_times) == 0L || length(ms$node_schedule_ode_step_sizes) == 0L ||
      is.null(ms$node_schedule_times)) {
    return(NULL)
  }
  list(node_schedule_times = ms$node_schedule_times, ode_times = ms$node_schedule_ode_times,
       ode_step_sizes = ms$node_schedule_ode_step_sizes,
       ode_soil_input_slopes = ms$node_schedule_ode_soil_input_slopes,
       ode_soil_substep_ends = ms$node_schedule_ode_soil_substep_ends,
       max_patch_lifetime = ms$p$max_patch_lifetime)
}

## Solve `community` on a recorded discretisation: every run replays its node
## schedule and ODE steps, with no refinement, so a finite difference between
## two such solves differs only in the model. plant's offspring_production_tol
## (the tolerance of its root-finds for a newborn's height and the light
## compensation point) is held to the cube of the step, as the stencil's
## equilibrium is (stencil_community()); at plant's default (1e-8) the
## replayed differences still scatter at the 1e-3 level.
##
## A discretisation spans one patch lifetime. A community with another (a
## parameter Jacobian's column in max_patch_lifetime) cannot replay it and
## refines its own, as both points of that column then do. Without a step (an
## equilibrium solve's own rounds) the tolerance is left as it is; without
## `replay_steps` only the node schedule is frozen, and each run integrates
## freely on it.
plant_community_freeze_discretisation <- function(community, discretisation, step = NULL,
                                                  replay_steps = TRUE) {
  if (!identical(discretisation$max_patch_lifetime, community$model_support$p$max_patch_lifetime)) {
    discretisation <- NULL
  }
  if (!is.null(discretisation) && !replay_steps) {
    discretisation[c("ode_times", "ode_step_sizes")] <- list(numeric(0), numeric(0))
    discretisation[c("ode_soil_input_slopes", "ode_soil_substep_ends")] <- list(list(), list())
  }
  community$model_support$frozen <- discretisation
  if (!is.null(step)) {
    ctrl <- community$model_support$plant_control
    ctrl$offspring_production_tol <- min(ctrl$offspring_production_tol, step^3)
    community$model_support$plant_control <- ctrl
  }
  community
}


## Support code:
plant_community_make_demography_runner <- function(community) {
  pretty_num_collapse <- function(x, collapse = ", ") {
    paste0("{", paste(prettyNum(x), collapse = collapse), "}")
  }

  p <- plant_community_parameters(community)
  ctrl <- community$model_support$plant_control  # plant SCM control() for run()

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
  # on a frozen discretisation every run replays it, never refining
  frozen <- !is.null(community$model_support$frozen)

  function(birth_rates) {
    # if a runner has diverged significantly in the last iteration then
    # reset the schedule to it's defaults and rebuild
    if (!frozen && any(abs(birth_rates - last_arrival_rates) > large_offspring_arriving_change)) {
      p$node_schedule_times <- default_schedule_times
    }

    # set birth rates
    for (i in seq_along(p$strategies)) {
      p$strategies[[i]]$birth_rate_y <- birth_rates[i]
    }

    # update schedule - starts from prev. schedule so should be fast for fine scale resolution
    # run() replays any ODE schedule p carries, and the one recorded at
    # the previous birth rates is stale, so integrate freely.
    if (frozen) {
      # a node schedule frozen without its steps: each run integrates freely,
      # never on the steps the previous run recorded
      if (length(community$model_support$frozen$ode_times) == 0L) {
        p <- plant_clear_ode_schedule(p)
      }
      res <- run(p, env = plant_community_environment(community), control = ctrl)
    } else {
      p <- plant_clear_ode_schedule(p)
      res <- run(p, env = plant_community_environment(community), control = ctrl, refine_schedule = TRUE)
    }
    p_new <- res$scenario$p
    offspring_production <- res$offspring_production


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
  community$model_support$node_schedule_ode_step_sizes <- p$ode_step_sizes
  community$model_support$node_schedule_ode_soil_input_slopes <- p$ode_soil_input_slopes
  community$model_support$node_schedule_ode_soil_substep_ends <- p$ode_soil_substep_ends
  community$fitness_points <- NULL

  attr(community, "progress") <- util_rbind_list(e$history)
  attr(community, "converged") <- converged
  community
}


plant_community_viable_bounds <- function(community) {
  if (length(community) > 0) {
    stop("You don't want to run this on an existing community")
  }
  community$bounds <- community_viable_fitness(community)
  community
}

plant_community_update_fitness_function <- function(community) {
  
  p <- plant_community_parameters(community)
  hyperpar <- plant_community_hyperpar(community)

  ctrl <- community$model_support$plant_control

  if (length(p$strategies) > 0L) {
    # if there's a resident, use the saved environment to calculate mutant
    # fitness; run() replays the resident's schedule when p carries one
    scm <- run(p, env = plant_community_environment(community), control = ctrl)$solver
    community$resident_fitness <- log(scm$patch$net_reproduction_ratios)
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
        ret <- log(scm$patch$net_reproduction_ratios)
      } else {
        # otherwise just run mutants with zero birth rate, integrating freely
        # (run() would otherwise take any schedule p carries)
        p_mutants <- plant_clear_ode_schedule(p_mutants)
        res <- run(p_mutants, env = plant_community_environment(community), control = ctrl)
        ret <- log(res$net_reproduction_ratios)
      }
      
    ret
  }

  ## plant's rare-mutant derivative (plant#650 E2), where the harness asks for
  ## it and there is a resident to hold the canopy at; harness_fd() fills the
  ## finite difference otherwise
  community$fitness_derivatives <- list()
  if (identical(community$harness$fitness_gradient, "model") && length(p$strategies) > 0L) {
    f <- community$fitness_function
    trait_names <- community$trait_names
    community$fitness_derivatives$gradient <- function(y, control) {
      plant_mutant_gradient(scm, p, hyperpar, y, trait_names, f, control)
    }
  }

  community
}

## The gradient of log invasion fitness at the mutants y (one per row) in the
## resident's recorded canopy, from plant's reverse sweep: the mutants share
## one replay, each the invader numbered by its row, and each costs one sweep.
## Taken through the hyperparameter function, so a column is the total
## derivative by the trait as the fitness function sets it (without it, lma's
## column is a partial, 5x off on FF16). A row plant refuses (a mutant with
## fewer than two cohorts, say) or answers non-finitely is a finite
## difference of the fitness function instead. plant refuses the height
## coordinate itself, by name.
plant_mutant_gradient <- function(scm, p, hyperpar, y, trait_names, f, control) {
  n <- nrow(y)
  k <- length(trait_names)
  p_mutants <- plant::add_mutant(p, trait_matrix(y, trait_names), hyperpar, birth_rate = rep(0, n))
  # scm is the fitness function's too: a replay left recording would keep
  # every later mutant's states
  recording <- scm$keep_states
  on.exit(scm$keep_states <- recording)
  scm$keep_states <- TRUE
  scm$run_mutant(p_mutants)
  metrics <- paste0(seq_len(n), ".net_reproduction_ratio")
  columns <- paste0(rep(seq_len(n), each = k), ".", rep(trait_names, n))
  sg <- plant::stand_gradient(scm, metrics = metrics, traits = columns, hyperpar = hyperpar)
  r0 <- as.numeric(sg$value[metrics])
  out <- matrix(NA_real_, n, k)
  for (i in seq_len(n)) {
    out[i, ] <- sg$gradient[metrics[i], paste0(i, ".", trait_names)] / r0[i]
  }
  value <- log(r0)
  bad <- plant::stand_gradient_refused(sg)[metrics] | !is.finite(rowSums(out)) | !is.finite(value)
  if (any(bad)) {
    fd <- fd_fitness_gradient(f, y[bad, , drop = FALSE], control)
    out[bad, ] <- fd
    value[bad] <- as.numeric(f(y[bad, , drop = FALSE]))
  }
  attr(out, "value") <- value
  out
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
