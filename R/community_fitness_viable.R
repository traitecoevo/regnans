plant_log_viable <- function(...) {
  plant_log_info(..., routine = "viable")
}

plant_log_inviable <- function(...) {
  plant_log_info(..., routine = "inviable")
}

##' Compute the region of positive (fundamental) fitness for a community.
##'
##' Finds the trait region over which a strategy has positive invasion fitness
##' into the community's environment, as the smallest box of trait bounds that
##' contains it. This is normally run on an \emph{empty} community, in which
##' case it is the fundamental niche: where a lone strategy can persist with no
##' competitors.
##'
##' Along each trait the box runs as far as the region does: it is the interval
##' over which the \emph{profile} of fitness, its maximum over the other traits
##' with this one held fixed, stays positive. That interval is bracketed and
##' refined by root-finding outwards from the fitness maximum, on the
##' community's trait scale. With one trait the profile is the fitness itself;
##' with \code{k} traits each profile value is a \code{k - 1} trait
##' maximisation, which dominates the cost.
##'
##' @title Compute region of positive fitness for a community
##' @param community A \code{community} object (usually empty).
##' @param x Trait values inside the region to search out from. Defaults to the
##'   fitness maximum within the community's bounds (\code{\link{max_fitness}}),
##'   which is also used when fitness at \code{x} is negative.
##' @param dx Initial step outwards when bracketing each edge, on the trait
##'   scale; it doubles until fitness turns negative or a bound is reached.
##' @return A bounds matrix (\code{lower}/\code{upper} columns, one row per
##'   trait), or \code{NULL} if fitness is nowhere positive within the bounds.
##'   An edge that reaches the community's bounds is returned as that bound.
##' @export
##' @author Rich FitzJohn, Daniel Falster
community_viable_fitness <- function(community, x = NULL, dx = 1) {
  bounds <- check_bounds(community$bounds)
  traits <- rownames(bounds)
  k <- length(traits)

  if (is.null(community$fitness_function)) {
    community <- community_update_fitness_function(community)
  }
  fitness <- function(y) community$fitness_function(trait_matrix(y, traits))

  tf <- community_trait_transform(community)
  lb <- tf$fwd(bounds)
  if (!all(is.finite(lb))) {
    stop("community_viable_fitness needs finite bounds on the trait scale")
  }

  if (!is.null(x)) {
    x <- check_point(x, bounds)
  }
  if (is.null(x) || fitness(x) < 0) {
    x <- max_fitness(community, bounds = bounds)
    w <- attr(x, "fitness")
    plant_log_viable(sprintf("\t...searching out from max fitness at %s (w=%2.5f)",
                             paste(formatC(x), collapse = ", "), w))
    if (w < 0) {
      return(NULL)
    }
  }
  x <- as.numeric(x)

  ret <- t(vapply(seq_len(k), function(i) {
    h <- viable_profile(fitness, i, x, bounds, tf)
    tf$inv(positive_1d(function(z) h(tf$inv(z)), tf$fwd(x[[i]]), dx,
                       lower = lb[i, 1], upper = lb[i, 2]))
  }, numeric(2)))
  dimnames(ret) <- list(traits, c("lower", "upper"))
  ret
}

## Fitness profile along trait i: t -> the maximum of fitness over the other
## traits with trait i held at t. Each maximisation starts from the previous
## one's answer, since successive values of t along a bracket sit close
## together. With one trait there is nothing to maximise over.
viable_profile <- function(fitness, i, x, bounds, tf) {
  if (length(x) == 1L) {
    return(function(t) fitness(t))
  }
  others <- x[-i]
  function(t) {
    g <- function(y_others) {
      y <- numeric(length(x))
      y[i] <- t
      y[-i] <- y_others
      fitness(y)
    }
    fit <- maximize_scaled(g, others, bounds[-i, , drop = FALSE], tf)
    others <<- fit$par
    fit$value
  }
}

## --- pure numerical root helpers (reimplemented, no plant dependency) -------
## Find the interval around `x` (which must have f(x) >= 0) where f crosses zero,
## i.e. the bounds of the positive-fitness region. `positive_1d_bracket` brackets
## the sign change in each direction; `positive_1d` then refines with uniroot.

positive_1d <- function(f, x, dx, lower = -Inf, upper = Inf, tol = 1e-3) {
  root <- function(b, type) {
    xs <- b[[type]]$x
    fx <- b[[type]]$fx
    if (prod(fx[1:2]) < 0) {
      ## suppressWarnings: -Inf replaced by maximally negative value, which is OK.
      suppressWarnings(uniroot(f, xs,
                               f.lower = fx[[1]], f.upper = fx[[2]],
                               tol = tol)$root)
    } else {
      if (type == "lower") xs[[1]] else xs[[2]]
    }
  }
  b <- positive_1d_bracket(f, x, dx, lower, upper)
  c(root(b, "lower"), root(b, "upper"))
}

positive_1d_bracket <- function(f, x, dx, lower, upper, grow = 2) {
  fx <- f(x)
  if (fx < 0) {
    stop("Don't yet support doing this with no positive values")
  }

  bracket <- function(x, dx, bound) {
    cleanup <- function(x, x_next, fx, fx_next) {
      if (dx < 0) {
        x <- c(x_next, x)
        fx <- c(fx_next, fx)
      } else {
        x <- c(x, x_next)
        fx <- c(fx, fx_next)
      }
      list(x = x, fx = fx)
    }
    hit_bounds <- FALSE
    repeat {
      x_next <- x + dx
      if ((dx < 0 && x_next < bound) || (dx > 0 && x_next > bound)) {
        x_next <- bound
        hit_bounds <- TRUE
      }
      fx_next <- f(x_next)
      if (fx_next < 0 || hit_bounds) {
        return(cleanup(x, x_next, fx, fx_next))
      } else {
        x <- x_next
        fx <- fx_next
        dx <- dx * grow
      }
    }
  }

  list(lower = bracket(x, -dx, lower),
       upper = bracket(x, dx, upper))
}
