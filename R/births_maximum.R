community_new_types_maximum_fitness <- function(sys, control) {
  
  ## This bit of weirdness exists so that extra information from the
  ## fitness search, such as approxiate fitness points, gets copied
  ## back so we can work with it.
  empty <- function(sys, m=NULL) {
    ret <- trait_matrix(numeric(0), sys$trait_names)
    ret <- copy_attributes(m, ret, exclude=c("dim", "dimnames"))
    attr(ret, "done") <- TRUE
    ret
  }

  if (is.null(sys$bounds)) {
    stop("Maximum-fitness births need trait bounds on the community to search within")
  }

  ret <- find_max_fitness(sys, control)

  if (attr(ret, "fitness") < control$eps_fitness_invasion) {
    plant_log_max_fitness(paste0("Best point had nonpositive fitness: ",
                                 attr(ret, "fitness")))
    return(empty(sys, ret))
  }

  if (length(sys) > 0L) {
    tf <- community_trait_transform(sys)
    i <- closest(tf$fwd(drop(ret)), tf$fwd(sys$traits), tf$fwd(sys$bounds))
    if (attr(i, "distance") < control$eps_too_close) {
      plant_log_max_fitness("Best point too close to existing")
      return(empty(sys, ret))
    }
  }

  ret
}

## This is fundamentally a really hard problem because we want to
## chase a lot of data back to the main object.  In particular:
##
## * We might not have an up-to-date cohort schedule or ode times: we
##   want to set them.
## * We want to get the full approximate fitness landscape
find_max_fitness <- function(sys, control) {

  plant_log_assembler("Finding maximum in fitness landscape")

  bounds <- check_bounds(sys$bounds)
  if (nrow(bounds) == 1L) {
    find_max_fitness_1D(sys, control$eps_too_close)
  } else {
    find_max_fitness_2d(sys, control$eps_too_close)
  }
}

## Simple function - takes existing maximum and then runs optim. Risks missing maximum
find_max_fitness_1D <- function(sys, eps_too_close) {

  # option 1 - use existing points
  i <- which.max(sys$fitness_points$fitness)
  xx <- sys$fitness_points[, 1, drop = TRUE]

  # todo: option 2 - fit 1D GP and use this to find max
  # xx <- seq_log_range(sys$bounds, 500)
  # yy <- fitness_approximate(xx)

  # Polish root by optnmising with actual fitness function
  ## first find range over which to look 
  r <- xx[c(max(1, i - 1), min(i + 1, length(xx)))]
  f <- sys$fitness_function
  
  opt <- optimise(f, r, maximum = TRUE)

  ret <- trait_matrix(opt$maximum, sys$trait_names)
  attr(ret, "fitness") <- opt$objective
  ret
}


find_max_fitness_2d <- function(sys, eps_too_close, tol=1e-2) {
  tf <- community_trait_transform(sys)
  do_fit <- function(p) {
    f <- function(x) {
      w <- sys$fitness_function(x)
      if (is.finite(w)) w else -1e6   # keep the optimiser away from blow-ups
    }
    fit <- maximize_scaled(f, p, sys$bounds, tol, tf)
    ret <- trait_matrix(fit$par, sys$trait_names)
    attr(ret, "fitness") <- if (is.finite(fit$value)) fit$value else -Inf
    ret
  }
  check <- function(fit, X) {
    j <- closest(tf$fwd(drop(fit)), tf$fwd(X), tf$fwd(sys$bounds))
    w <- attr(fit, "fitness")
    d <- attr(j, "distance")
    plant_log_max_fitness(sprintf("\t...fitness: %s, distance: %s from %d",
                                  prettyNum(w), prettyNum(d), j))
    isTRUE(is.finite(w) && d > eps_too_close && w > 0.0)
  }

  centre <- tf$inv(rowMeans(tf$fwd(sys$bounds)))

  if (length(sys) == 0L) {
    ret <- do_fit(centre)
  } else {
    ## Multistart hill-climb on invasion fitness: from each resident and from
    ## the centre of trait space. Prefer the highest-fitness optimum that is
    ## viable and not too close to an existing resident; otherwise return the
    ## global best so the caller (community_new_types_maximum_fitness) can decide
    ## it is "done".
    X <- sys$traits
    starts <- rbind(X, matrix(centre, nrow = 1))
    fits <- lapply(seq_len(nrow(starts)), function(i) do_fit(starts[i, ]))
    w <- vnapply(fits, function(f) attr(f, "fitness"))
    ord <- order(w, decreasing = TRUE)

    ret <- fits[[ord[[1]]]]          # global best optimum (default)
    for (i in ord) {
      if (check(fits[[i]], X)) {      # viable and distinct from residents
        ret <- fits[[i]]
        break
      }
    }
  }
  ret
}
