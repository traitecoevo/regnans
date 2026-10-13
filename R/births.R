assembler_births <- function(obj) {
  to_add <- assembler_new_types(obj)
  obj$births <- to_add
  obj$refused <- NULL

  ## A birth the control refuses ends the assembly, with the reason kept:
  ## the community would outgrow max_residents, or the newcomer sits on the
  ## edge of the bounds, where the box rather than the model put it.
  if (nrow(to_add) > 0) {
    refused <- births_refused(obj$community, to_add, obj$control)
    if (!is.null(refused)) {
      plant_log_births(paste0("Refusing birth: ", refused$reason))
      obj$refused <- refused
      obj$done <- TRUE
      obj$community$done <- TRUE
      return(obj)
    }
  }

  ## Now we add things:
  if (nrow(to_add) > 0) {
    plant_log_births_new_types(to_add)
    i <- should_move(obj, to_add)
    if (length(i) == 1L) {
      obj <- assembler_births_try_move(obj, to_add, i)
    } else {
      rain <- initial_birth_rate(attr(to_add, "fitness"), obj)
      obj$community <- community_add(obj$community, to_add, rain)
    }
  }

  ## This is a bit ugly: we only really want to respond to no births
  ## with being done if we are running to equilibrium.  If we're using
  ## the simple step, there is no way of finishing here unless we look
  ## to a limited amount of change in the community since last step.
  if (has_attr(to_add, "done")) {
    obj$done <- (obj$control$run_type == "to_equilibrium" &&
                 attr(to_add, "done"))
    obj$community$done <- TRUE
  }
  obj
}

assembler_births_try_move <- function(obj, to_add, i) {
  community <- obj$community
  plant_log_births(paste0("Trying to move resident ", i))
  x_prev <- community$traits[i, , drop=FALSE]

  test <- community
  test$traits[i, ] <- to_add
  test <- community_demography(test)

  ## Now, try adding the previous case back in.  First, check
  ## the fitness:
  w_prev <- test$fitness_function(x_prev)
  if (w_prev > obj$control$eps_fitness_invasion) {
    ## Here, there is no need to run anything: the original
    ## resident has positive fitness and will increase
    plant_log_births("Adding original resident back")
    community <- community_add(test, x_prev, initial_birth_rate(w_prev, obj))
  } else {
    plant_log_births("Move was successful")
    community <- test
  }
  obj$community <- community
  obj
}

assembler_new_types <- function(obj) {
  community_new_types(obj$community, obj$control)
}

community_new_types <- function(community, control) {
  if (control$birth_type == "maximum") {
    community_new_types_maximum_fitness(community, control)
  } else if (control$birth_type == "stochastic") {
    community_new_types_stochastic(community, control)
  } else {
    stop("Unknown birth type ", control$birth_type)
  }
}

## First order birth rate approximation based on fitness.
initial_birth_rate <- function(fitness, obj) {
  if (is.null(fitness) || obj$control$run_type == "single_step") {
    birth_rate <- NULL
  } else {
    birth_rate <- pmin(exp(fitness), obj$control$max_birth_rate_initial)
    birth_rate <- pmax(birth_rate, obj$control$min_birth_rate_initial)
  }
  birth_rate
}

plant_log_births_new_types <- function(to_add) {
  str <- format_community_state(to_add, "-", sprintf("\t%s\t", c("i", seq_along(to_add))))
  plant_log_births(paste(c("Proposed new type(s):", str), collapse="\n"),
                   to_add=to_add)
}

## NULL, or a list(status, reason) saying why the control refuses these births.
births_refused <- function(community, to_add, control) {
  m <- length(community)
  if (m + nrow(to_add) > control$max_residents) {
    return(list(status = "capped", reason = sprintf(
      "a birth would take the community beyond %s residents", format(control$max_residents))))
  }
  if (control$edge_tol > 0 && m > 0L) {
    tf <- community_trait_transform(community)
    lb <- tf$fwd(community$bounds)
    z <- matrix(tf$fwd(to_add), nrow = nrow(to_add))
    width <- rep(lb[, 2] - lb[, 1], each = nrow(z))
    ## the edge is a fraction of the width, so a trait whose bounds are not
    ## finite on the trait scale has none
    on_edge <- is.finite(width) & (z - rep(lb[, 1], each = nrow(z)) < control$edge_tol * width |
      rep(lb[, 2], each = nrow(z)) - z < control$edge_tol * width)
    if (any(on_edge)) {
      return(list(status = "edge", reason = "the newcomer's maximum lies on the edge of the bounds"))
    }
  }
  NULL
}
