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

## The highest-fitness trait combination that could found a new resident.
##
## Each candidate is a local maximisation of invasion fitness on the trait
## scale (maximize_scaled). Where the community holds a sampled fitness
## landscape, one search is confined to the samples either side of its best
## point; otherwise the searches start from each resident and from the centre
## of trait space. The highest-fitness candidate that is positive and not within
## eps_too_close of a resident is returned; failing that the global best, so
## community_new_types_maximum_fitness() can decide the assembly is done.
find_max_fitness <- function(sys, control) {

  plant_log_assembler("Finding maximum in fitness landscape")

  bounds <- check_bounds(sys$bounds)
  tf <- community_trait_transform(sys)
  search <- function(start, region) {
    fit <- maximize_scaled(sys$fitness_function, start, region, tf)
    ret <- trait_matrix(fit$par, sys$trait_names)
    attr(ret, "fitness") <- if (is.finite(fit$value)) fit$value else -Inf
    ret
  }

  points <- sys$fitness_points
  if (!is.null(points)) {
    xx <- points[, 1, drop = TRUE]
    i <- which.max(points$fitness)
    region <- bounds
    region[1, ] <- xx[c(max(1, i - 1), min(i + 1, length(xx)))]
    fits <- list(search(xx[i], region))
  } else {
    centre <- tf$inv(rowMeans(tf$fwd(bounds)))
    starts <- rbind(sys$traits, centre, deparse.level = 0)
    fits <- lapply(seq_len(nrow(starts)),
                   function(j) search(starts[j, ], bounds))
  }

  w <- vnapply(fits, function(fit) attr(fit, "fitness"))
  ord <- order(w, decreasing = TRUE)
  if (length(sys) == 0L) {
    return(fits[[ord[[1]]]])
  }
  for (j in ord) {
    k <- closest(tf$fwd(drop(fits[[j]])), tf$fwd(sys$traits), tf$fwd(bounds))
    d <- attr(k, "distance")
    plant_log_max_fitness(sprintf("\t...fitness: %s, distance: %s from %d",
                                  prettyNum(w[[j]]), prettyNum(d), k))
    if (is.finite(w[[j]]) && w[[j]] > 0 && d > control$eps_too_close) {
      return(fits[[j]])
    }
  }
  fits[[ord[[1]]]]
}
