# Pairwise invasibility: the invasion-fitness surface s(y; x) over residents x
# and mutants y, its zero contours, and the two-resident pictures built on it.
#
# The resident axis is the expensive one (each resident is a demographic
# equilibrium solve); mutants against a solved resident are vectorised and
# cheap. So community_pip() spends its effort there. Equilibrium density varies
# smoothly along the trait axis, so a coarse set of residents is solved first,
# in parallel, and every other resident starts its solve from the density
# interpolated between them -- close to the answer, and independent of its
# neighbours, so the whole set runs in parallel too. The resident grid is then
# refined only where the zero contours bend. Along the cheap mutant axis the
# contours are located exactly, by Newton on the model's fitness gradient,
# rather than read off a sign grid.
#
# The pairwise invasibility plot is the sign of the surface; the mutual
# invasibility plot compares s(y; x) with s(x; y) and needs no new solves; the
# trait-evolution plot (community_tep) solves each mutually invasible pair as a
# dimorphic community and reports the selection gradient on both residents.

##' Control how a pairwise invasibility surface is computed.
##'
##' @title Pairwise-invasibility settings
##' @param control A list of values to modify from the defaults:
##' \describe{
##'   \item{\code{n_resident}}{initial number of residents across the bounds,
##'     spaced on the trait scale. Each is an equilibrium solve.}
##'   \item{\code{n_mutant}}{number of mutant values per resident. Cheap; this
##'     sets the resolution of the shaded surface and of contour detection.}
##'   \item{\code{seed}}{how each resident's equilibrium solve is started.
##'     \code{"interpolate"} (default): \code{n_coarse} residents are solved
##'     first and the rest start from the density interpolated between them, on
##'     the trait scale and in log density. \code{"neighbour"}: residents are
##'     solved in order, each from the previous one's density, in one contiguous
##'     chunk per worker. \code{"cold"}: every solve starts from the community's
##'     default birth rate.}
##'   \item{\code{n_coarse}}{number of residents in the first pass of
##'     \code{seed = "interpolate"}.}
##'   \item{\code{refine}}{rounds of resident refinement: a resident is inserted
##'     between neighbours where the zero contours bend --- a contour appears or
##'     disappears between them, or its second difference across three
##'     consecutive residents exceeds \code{refine_tol} of the mutant range.}
##'   \item{\code{refine_tol}}{contour curvature (second difference of a
##'     crossing across three residents, as a fraction of the mutant range on
##'     the trait scale) that triggers refinement.}
##'   \item{\code{tol}}{convergence tolerance of the contour root-finding, on
##'     the trait scale.}
##'   \item{\code{parallel}}{evaluate residents through the active
##'     \code{future::plan()} when one is set.}
##' }
##' @return A list of the settings above.
##' @author Daniel Falster
##' @export
pip_control <- function(control = NULL) {
  defaults <- list(
    n_resident = 41L,
    n_mutant   = 201L,
    seed       = "interpolate",
    n_coarse   = 9L,
    refine     = 2L,
    refine_tol = 0.02,
    tol        = 1e-8,
    parallel   = TRUE
  )
  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse = ", "))
  }
  ret <- modifyList(defaults, control)
  ret$seed <- match.arg(ret$seed, c("interpolate", "neighbour", "cold"))
  for (nm in c("n_resident", "n_mutant", "n_coarse", "refine")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || v < 0 || v != round(v)) {
      stop(nm, " must be a whole number")
    }
    ret[[nm]] <- as.integer(v)
  }
  if (ret$n_resident < 2L || ret$n_mutant < 3L || ret$n_coarse < 2L) {
    stop("n_resident and n_coarse must be at least 2 and n_mutant at least 3")
  }
  for (nm in c("refine_tol", "tol")) {
    v <- ret[[nm]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0) {
      stop(nm, " must be a single positive number")
    }
  }
  if (!is.logical(ret$parallel) || length(ret$parallel) != 1L || is.na(ret$parallel)) {
    stop("parallel must be TRUE or FALSE")
  }
  ret
}

## Even spacing on the trait scale between the bounds of a single trait.
pip_grid <- function(tf, bounds, n) {
  lo <- tf$fwd(bounds[1, 1])
  hi <- tf$fwd(bounds[1, 2])
  if (!is.finite(lo) || !is.finite(hi)) {
    stop("community_pip needs finite bounds on the trait scale")
  }
  tf$inv(seq(lo, hi, length.out = n))
}

## Newton on f with derivative df inside a bracket [lo, hi] that contains a
## sign change, from the start y0. A Newton step that leaves the bracket, or
## fails to shrink it, is replaced by a bisection step, so this converges
## whatever the start and takes the quadratic steps where it can.
newton_bracketed <- function(f, df, lo, hi, y0, tol, maxit = 50L) {
  flo <- f(lo)
  fhi <- f(hi)
  if (flo == 0) return(lo)
  if (fhi == 0) return(hi)
  if (sign(flo) == sign(fhi)) {
    stop("newton_bracketed: no sign change in the bracket")
  }
  y <- min(max(y0, lo), hi)
  for (i in seq_len(maxit)) {
    fy <- f(y)
    if (fy == 0) return(y)
    if (sign(fy) == sign(flo)) { lo <- y; flo <- fy } else { hi <- y; fhi <- fy }
    if (hi - lo <= tol) return((lo + hi) / 2)
    d <- df(y)
    step <- if (is.finite(d) && d != 0) fy / d else NA_real_
    candidate <- y - step
    if (!is.finite(candidate) || candidate <= lo || candidate >= hi) {
      candidate <- (lo + hi) / 2
    }
    y <- candidate
  }
  (lo + hi) / 2
}

## Zero crossings of s(y; x) in the mutant y for one solved resident x, on the
## trait scale. The diagonal root y = x is known; dividing it out exposes any
## second root that sits close to it (near a singular strategy the two are
## separated by about 2 g / H, which a grid cannot resolve), and that root is
## started from exactly that estimate. Remaining roots are bracketed by sign
## changes of the deflated values on the mutant grid.
pip_crossings <- function(solved, z_mutant, fitness, tf, tol) {
  x <- as.numeric(solved$traits[1, 1])
  zx <- tf$fwd(x)
  f <- function(z) as.numeric(solved$fitness_function(tf$inv(z)))
  ## d/dz of s on the trait scale: ds/dy * dy/dz
  df <- function(z) {
    y <- tf$inv(z)
    g <- as.numeric(community_fitness_gradient(solved, y))
    if (identical(tf$scale, "log")) g * y else g
  }
  deflated <- function(z) f(z) / (z - zx)
  d_deflated <- function(z) {
    s <- f(z)
    (df(z) * (z - zx) - s) / (z - zx)^2
  }

  roots <- zx
  away <- abs(z_mutant - zx) > tol
  zz <- z_mutant[away]
  vals <- fitness[away] / (zz - zx)
  ok <- is.finite(vals)
  zz <- zz[ok]
  vals <- vals[ok]
  if (length(zz) < 2L) {
    return(roots)
  }

  ## the near-diagonal root from the local quadratic s ~ g (y-x) + H (y-x)^2 / 2
  g <- df(zx)
  H <- {
    h <- as.numeric(community_fitness_hessian(solved))
    if (identical(tf$scale, "log")) h * x^2 + g else h
  }
  z_near <- if (is.finite(g) && is.finite(H) && H != 0) zx - 2 * g / H else NA_real_

  flips <- which(sign(vals[-1]) != sign(vals[-length(vals)]) & vals[-1] != 0)
  for (j in flips) {
    lo <- zz[j]
    hi <- zz[j + 1L]
    start <- if (is.finite(z_near) && z_near > lo && z_near < hi) z_near else (lo + hi) / 2
    r <- tryCatch(newton_bracketed(deflated, d_deflated, lo, hi, start, tol),
                  error = function(e) NA_real_)
    if (is.finite(r)) roots <- c(roots, r)
  }
  ## a near-diagonal root that fell inside the same grid cell as the diagonal
  ## itself has no bracket above; polish it from the quadratic estimate
  if (is.finite(z_near) && abs(z_near - zx) > tol &&
      z_near > min(z_mutant) && z_near < max(z_mutant) &&
      !any(abs(roots - z_near) <= max(tol, abs(2 * g / H)) & roots != zx)) {
    cell <- diff(range(z_mutant)) / (length(z_mutant) - 1L)
    lo <- if (z_near < zx) max(z_near - cell, min(z_mutant)) else zx + tol
    hi <- if (z_near < zx) zx - tol else min(z_near + cell, max(z_mutant))
    if (hi > lo && sign(deflated(lo)) != sign(deflated(hi))) {
      r <- tryCatch(newton_bracketed(deflated, d_deflated, lo, hi, z_near, tol),
                    error = function(e) NA_real_)
      if (is.finite(r) && !any(abs(roots - r) <= tol)) roots <- c(roots, r)
    }
  }
  ## a contour through the diagonal (the singular strategy itself) is found
  ## again by the deflated search as a root within tol of zx
  roots <- sort(roots)
  roots[c(TRUE, diff(roots) > tol)]
}

## Solve one resident and collect everything the surface needs from it. Runs
## on a worker, so it returns plain numbers only (a plant harness's fitness
## closure holds an external pointer that cannot cross to another process).
## n_evals is the number of demography-runner calls the equilibrium took: the
## cost of the solve, which the seeding strategy is there to reduce.
pip_resident <- function(base, x, z_mutant, tf, control, birth_rate = NULL) {
  solved <- base |>
    community_add(trait_matrix(x, base$trait_names), birth_rate = birth_rate) |>
    community_demography()
  fitness <- as.numeric(solved$fitness_function(tf$inv(z_mutant)))
  crossings <- pip_crossings(solved, z_mutant, fitness, tf, control$tol)
  list(resident = x, birth_rate = as.numeric(solved$birth_rate),
       n_evals = NROW(attr(solved, "progress")),
       fitness = fitness, crossings = crossings)
}

## Starting densities for new residents, interpolated in log density along the
## trait scale between the residents already solved (density is smooth in the
## trait, so this lands close to the equilibrium).
pip_interpolate_seeds <- function(solved, z_new, tf) {
  if (length(solved) == 0L) return(rep(NA_real_, length(z_new)))
  z <- tf$fwd(vapply(solved, `[[`, numeric(1), "resident"))
  n <- vapply(solved, `[[`, numeric(1), "birth_rate")
  ok <- is.finite(n) & n > 0
  if (sum(ok) == 0L) return(rep(NA_real_, length(z_new)))
  if (sum(ok) == 1L) return(rep(n[ok], length(z_new)))
  exp(stats::approx(z[ok], log(n[ok]), xout = z_new, rule = 2, ties = mean)$y)
}

## Solve a set of residents in parallel chunks. With seeds supplied (one per
## resident, NA where unknown) each solve is independent; without them, each
## chunk runs in order and warm-starts from its previous resident.
pip_solve_residents <- function(base, residents, seeds, z_mutant, tf, control) {
  chunks <- regnans_chunks(length(residents), control$parallel)
  solved <- regnans_map(chunks, function(idx) {
    out <- vector("list", length(idx))
    previous <- NULL
    for (j in seq_along(idx)) {
      i <- idx[j]
      seed <- if (is.finite(seeds[i]) && seeds[i] > 0) seeds[i] else previous
      out[[j]] <- pip_resident(base, residents[i], z_mutant, tf, control, seed)
      br <- out[[j]]$birth_rate
      previous <- if (control$seed != "cold" && length(br) == 1L && is.finite(br) && br > 0) br else NULL
    }
    out
  }, parallel = control$parallel)
  unlist(solved, recursive = FALSE)
}

## The two-pass solve of the initial resident grid for seed = "interpolate":
## a coarse, evenly spread subset first, then everything else seeded from it.
pip_solve_grid <- function(base, resident, z_mutant, tf, control) {
  n <- length(resident)
  none <- rep(NA_real_, n)
  if (control$seed != "interpolate" || control$n_coarse >= n) {
    return(pip_solve_residents(base, resident, none, z_mutant, tf, control))
  }
  coarse <- unique(round(seq(1, n, length.out = control$n_coarse)))
  solved <- pip_solve_residents(base, resident[coarse], none[coarse], z_mutant, tf, control)
  rest <- setdiff(seq_len(n), coarse)
  seeds <- pip_interpolate_seeds(solved, tf$fwd(resident[rest]), tf)
  c(solved, pip_solve_residents(base, resident[rest], seeds, z_mutant, tf, control))
}

## Which intervals between consecutive residents deserve a new resident: those
## where a contour appears or disappears (the crossing counts differ), and
## those on either side of a resident at which a contour bends -- the second
## difference of a crossing across three consecutive residents exceeds
## refine_tol of the mutant range. A straight contour is never refined, however
## steep; a count change always is, since that is where a contour enters,
## leaves, or meets the diagonal at a singular strategy.
pip_refine_intervals <- function(crossings, span, control) {
  n <- length(crossings)
  if (n < 2L) return(integer(0))
  counts <- lengths(crossings)
  mark <- counts[-1] != counts[-n]
  for (i in seq_len(n)[-c(1L, n)]) {
    if (counts[i - 1L] == counts[i] && counts[i] == counts[i + 1L]) {
      bend <- abs(crossings[[i + 1L]] - 2 * crossings[[i]] + crossings[[i - 1L]])
      if (any(bend > control$refine_tol * span)) mark[c(i - 1L, i)] <- TRUE
    }
  }
  which(mark)
}

##' Pairwise invasibility surface and its zero contours.
##'
##' For each resident trait value the community is reduced to that single
##' resident, solved to demographic equilibrium, and the invasion fitness of
##' every mutant value is evaluated against it. The sign of the resulting
##' surface \eqn{s(y; x)} is the pairwise invasibility plot: it shows where
##' selection points, where the singular strategies are, and whether each is
##' an ESS or a branching point.
##'
##' Residents are the expensive axis (one equilibrium solve each) and are
##' handled accordingly: a coarse subset is solved first, through the active
##' \code{future::plan()}, and every other resident starts its solve from the
##' equilibrium density interpolated between them, so each lands near its
##' answer and all of them run in parallel; after the initial grid the resident
##' axis is refined where neighbouring residents disagree about their zero
##' contours (see \code{\link{pip_control}}). \code{residents$n_evals}
##' records what each solve cost. Along the cheap mutant axis the contours
##' \eqn{s(y; x) = 0} are located exactly by Newton's method on the model's
##' fitness gradient, with the diagonal root \eqn{y = x} divided out so that the
##' second contour is found even where it hugs the diagonal near a singular
##' strategy.
##'
##' @title Pairwise invasibility plot
##' @param community A single-trait \code{community}; any residents it carries
##' are discarded.
##' @param resident Resident trait values for the initial grid. Defaults to
##' \code{control$n_resident} values spaced on the trait scale across
##' \code{bounds}.
##' @param mutant Mutant trait values. Defaults to \code{control$n_mutant}
##' values across \code{bounds}.
##' @param bounds Trait bounds for the default grids; defaults to the
##' community's.
##' @param control See \code{\link{pip_control}}.
##' @return A \code{pip} object: a list with \code{surface} (a tibble of
##' \code{resident}, \code{mutant}, \code{fitness}), \code{contours} (a tibble
##' of \code{resident}, \code{mutant} on the zero set, including the diagonal),
##' \code{residents} (each resident's equilibrium \code{birth_rate}, the number
##' of demography evaluations its solve took, and its number of crossings), the
##' trait name and scale, the control used, and \code{elapsed} seconds.
##' \code{as.data.frame()} returns the surface. See \code{\link{plot.pip}},
##' \code{\link{pip_mutual}}, \code{\link{community_tep}}.
##' @author Daniel Falster
##' @export
community_pip <- function(community, resident = NULL, mutant = NULL,
                          bounds = community$bounds, control = pip_control()) {
  trait_names <- community$trait_names
  if (length(trait_names) != 1L) {
    stop("community_pip needs a single-trait community; this one has ",
         length(trait_names), " traits")
  }
  control <- pip_control(control)
  tf <- community_trait_transform(community)
  if (is.null(resident)) resident <- pip_grid(tf, bounds, control$n_resident)
  if (is.null(mutant)) mutant <- pip_grid(tf, bounds, control$n_mutant)
  resident <- sort(unique(as.numeric(resident)))
  mutant <- sort(unique(as.numeric(mutant)))
  z_mutant <- tf$fwd(mutant)
  span <- diff(range(z_mutant))

  plant_log_assembler(sprintf(
    "Pairwise invasibility for %s: %d residents x %d mutants, %d refinement rounds",
    trait_names, length(resident), length(mutant), control$refine))

  started <- proc.time()[["elapsed"]]
  base <- community_clear_residents(community)
  solved <- pip_solve_grid(base, resident, z_mutant, tf, control)

  for (round in seq_len(control$refine)) {
    x <- vapply(solved, `[[`, numeric(1), "resident")
    o <- order(x)
    solved <- solved[o]
    x <- x[o]
    z <- tf$fwd(x)
    n <- length(solved)
    where <- pip_refine_intervals(lapply(solved, `[[`, "crossings"), span, control)
    if (length(where) == 0L) break
    new_z <- (z[where] + z[where + 1L]) / 2
    plant_log_assembler(sprintf("  refinement round %d: %d residents added", round, length(new_z)))
    seeds <- if (control$seed == "cold") rep(NA_real_, length(new_z)) else pip_interpolate_seeds(solved, new_z, tf)
    solved <- c(solved, pip_solve_residents(base, tf$inv(new_z), seeds, z_mutant, tf, control))
  }

  x <- vapply(solved, `[[`, numeric(1), "resident")
  solved <- solved[order(x)]
  x <- sort(x)

  surface <- tibble::tibble(
    resident = rep(x, each = length(mutant)),
    mutant = rep(mutant, times = length(x)),
    fitness = unlist(lapply(solved, `[[`, "fitness"))
  )
  contours <- tibble::tibble(
    resident = rep(x, times = vapply(solved, function(s) length(s$crossings), integer(1))),
    mutant = tf$inv(unlist(lapply(solved, `[[`, "crossings")))
  )
  residents <- tibble::tibble(
    resident = x,
    birth_rate = vapply(solved, `[[`, numeric(1), "birth_rate"),
    n_evals = vapply(solved, `[[`, integer(1), "n_evals"),
    n_crossings = vapply(solved, function(s) length(s$crossings), integer(1))
  )

  structure(list(surface = surface, contours = contours, residents = residents,
                 trait_name = trait_names, trait_scale = tf$scale,
                 control = control,
                 elapsed = proc.time()[["elapsed"]] - started),
            class = "pip")
}

##' @export
as.data.frame.pip <- function(x, ...) as.data.frame(x$surface, ...)

##' @export
print.pip <- function(x, ...) {
  cat(sprintf("<pip: %s, %d residents x %d mutants, %d contour points; %d demography evaluations in %.2fs>\n",
              x$trait_name, nrow(x$residents), length(unique(x$surface$mutant)),
              nrow(x$contours), sum(x$residents$n_evals), x$elapsed))
  invisible(x)
}

## s(y; x_i) for arbitrary mutant values y, interpolated along the mutant axis
## of resident i's column. The mutant grid is fine and fitness is smooth in the
## mutant, so this is accurate wherever it is used (resident-against-resident
## values for the mutual-invasibility view).
pip_fitness_at <- function(pip, i, y) {
  col <- pip$surface[pip$surface$resident == pip$residents$resident[i], ]
  tf <- if (identical(pip$trait_scale, "log")) log else identity
  fin <- is.finite(col$fitness)
  stats::approx(tf(col$mutant[fin]), col$fitness[fin], xout = tf(y), rule = 2)$y
}

##' Mutually invasible pairs of residents.
##'
##' Two strategies can form a protected dimorphism when each invades the
##' other's equilibrium: \eqn{s(x_2; x_1) > 0} and \eqn{s(x_1; x_2) > 0}. Both
##' are read from the surface, so no further model evaluations are needed.
##'
##' @title Mutual invasibility
##' @param pip A \code{pip} from \code{\link{community_pip}}.
##' @return A tibble with one row per ordered pair of residents: \code{x1},
##' \code{x2}, \code{s12} (fitness of \code{x2} invading \code{x1}),
##' \code{s21}, and \code{mutual}.
##' @author Daniel Falster
##' @export
pip_mutual <- function(pip) {
  x <- pip$residents$resident
  n <- length(x)
  S <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) S[i, ] <- pip_fitness_at(pip, i, x)
  pairs <- expand.grid(i = seq_len(n), j = seq_len(n))
  out <- tibble::tibble(
    x1 = x[pairs$i], x2 = x[pairs$j],
    s12 = S[cbind(pairs$i, pairs$j)],
    s21 = S[cbind(pairs$j, pairs$i)]
  )
  out$mutual <- out$x1 != out$x2 & out$s12 > 0 & out$s21 > 0
  out
}

##' Trait-evolution plot: selection on protected dimorphisms.
##'
##' For every mutually invasible pair of residents the two are introduced
##' together, the dimorphic community is solved to demographic equilibrium
##' (seeded from the pair's monomorphic densities), and the selection gradient
##' on each resident is recorded. Drawn as arrows over the coexistence region
##' this is the trait-evolution plot of Geritz et al. (1998): it shows where a
##' dimorphism formed by branching is carried next, and whether a dimorphic
##' singular coalition (both gradients zero) exists. Pairs where one strategy is
##' driven out at the dimorphic equilibrium are dropped and counted in
##' \code{attr(., "excluded")}.
##'
##' Cost: one two-resident equilibrium solve per mutually invasible pair above
##' the diagonal, through the same parallel map as \code{\link{community_pip}}.
##'
##' @title Trait-evolution plot
##' @param community The \code{community} the surface was computed from.
##' @param pip A \code{pip} from \code{\link{community_pip}}.
##' @return A \code{tep} object: a tibble with one row per coexisting pair
##' (\code{x1 < x2}), the equilibrium state \code{n1}, \code{n2} and selection
##' gradients \code{g1}, \code{g2}. See \code{\link{plot.tep}}.
##' @author Daniel Falster
##' @export
community_tep <- function(community, pip) {
  trait_names <- community$trait_names
  if (!identical(pip$trait_name, trait_names)) {
    stop("The pip was computed for trait ", pip$trait_name,
         ", not for this community's ", trait_names)
  }
  mutual <- pip_mutual(pip)
  pairs <- mutual[mutual$mutual & mutual$x1 < mutual$x2, c("x1", "x2")]
  plant_log_assembler(sprintf("Trait-evolution plot for %s: %d coexisting pairs",
                              trait_names, nrow(pairs)))

  base <- community_clear_residents(community)
  extinct <- community$demography_control$equilibrium_extinct_birth_rate
  seed_of <- stats::setNames(pip$residents$birth_rate, pip$residents$resident)
  rows <- regnans_map(seq_len(nrow(pairs)), function(i) {
    x <- c(pairs$x1[i], pairs$x2[i])
    solved <- base |>
      community_add(trait_matrix(x, trait_names),
                    birth_rate = unname(seed_of[as.character(x)])) |>
      community_demography()
    n <- as.numeric(solved$birth_rate)
    if (any(!is.finite(n)) || any(n <= extinct)) return(NULL)
    g <- community_fitness_gradient(solved)
    tibble::tibble(x1 = x[1], x2 = x[2], n1 = n[1], n2 = n[2],
                   g1 = g[1, 1], g2 = g[2, 1])
  }, parallel = pip$control$parallel)
  kept <- !vapply(rows, is.null, logical(1))
  out <- if (any(kept)) do.call(rbind, rows[kept]) else
    tibble::tibble(x1 = numeric(0), x2 = numeric(0), n1 = numeric(0),
                   n2 = numeric(0), g1 = numeric(0), g2 = numeric(0))
  attr(out, "trait_name") <- trait_names
  attr(out, "trait_scale") <- pip$trait_scale
  attr(out, "excluded") <- sum(!kept)
  class(out) <- c("tep", class(out))
  out
}

## Cell edges for drawing an irregular grid: midpoints between neighbours on
## the trait scale, mapped back, so cells meet exactly on a log or linear axis.
pip_edges <- function(values, tf) {
  z <- tf$fwd(values)
  mid <- (z[-1] + z[-length(z)]) / 2
  step <- if (length(z) > 1L) c(mid[1] - z[1], z[length(z)] - mid[length(mid)]) else c(0.5, 0.5)
  lo <- c(z[1] - step[1], mid)
  hi <- c(mid, z[length(z)] + step[2])
  list(lo = tf$inv(lo), hi = tf$inv(hi))
}

pip_transform <- function(scale) {
  if (identical(scale, "log")) list(fwd = log, inv = exp, scale = "log")
  else list(fwd = identity, inv = identity, scale = "linear")
}

pip_axes <- function(scale, trait, xlab, ylab) {
  scales <- if (identical(scale, "log")) {
    list(ggplot2::scale_x_log10(expand = c(0, 0)), ggplot2::scale_y_log10(expand = c(0, 0)))
  } else {
    list(ggplot2::scale_x_continuous(expand = c(0, 0)), ggplot2::scale_y_continuous(expand = c(0, 0)))
  }
  c(scales, list(
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed"),
    ggplot2::coord_equal(),
    ggplot2::xlab(paste(xlab, trait)), ggplot2::ylab(paste(ylab, trait)),
    ggplot2::theme_classic(),
    ggplot2::theme(text = ggplot2::element_text(size = 16))))
}

##' Plot a pairwise or mutual invasibility surface.
##'
##' \code{type = "pip"} shades the (resident, mutant) cells where the mutant
##' can invade and overlays the exactly located zero contours;
##' \code{type = "mip"} shades the resident pairs that invade each other (see
##' \code{\link{pip_mutual}}). With \code{fill = "fitness"} the PIP shows the
##' fitness value on a diverging scale centred on zero instead of its sign.
##' Cells follow the (possibly refined) resident grid, so they are not all the
##' same width. The diagonal, where the mutant is the resident, is drawn
##' through.
##'
##' @title Plot a pairwise or mutual invasibility surface
##' @param x A \code{pip} from \code{\link{community_pip}}.
##' @param type \code{"pip"} or \code{"mip"}.
##' @param fill For \code{"pip"}: \code{"sign"} (default) or \code{"fitness"}.
##' @param contours Draw the zero contours (\code{"pip"} only).
##' @param ... Ignored.
##' @return A \code{ggplot} object.
##' @author Daniel Falster
##' @export
plot.pip <- function(x, type = c("pip", "mip"), fill = c("sign", "fitness"),
                     contours = TRUE, ...) {
  type <- match.arg(type)
  fill <- match.arg(fill)
  tf <- pip_transform(x$trait_scale)
  rx <- x$residents$resident
  ex <- pip_edges(rx, tf)

  if (type == "mip") {
    m <- pip_mutual(x)
    i <- match(m$x1, rx)
    j <- match(m$x2, rx)
    data <- tibble::tibble(xmin = ex$lo[i], xmax = ex$hi[i], ymin = ex$lo[j], ymax = ex$hi[j],
                           mutual = m$mutual)
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_rect(ggplot2::aes(xmin = .data[["xmin"]], xmax = .data[["xmax"]],
                                      ymin = .data[["ymin"]], ymax = .data[["ymax"]],
                                      fill = .data[["mutual"]])) +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "steelblue", `FALSE` = "white"),
                                 name = "coexist")
    return(p + pip_axes(x$trait_scale, x$trait_name, "resident 1", "resident 2"))
  }

  my <- sort(unique(x$surface$mutant))
  ey <- pip_edges(my, tf)
  i <- match(x$surface$resident, rx)
  j <- match(x$surface$mutant, my)
  data <- tibble::tibble(xmin = ex$lo[i], xmax = ex$hi[i], ymin = ey$lo[j], ymax = ey$hi[j],
                         fitness = x$surface$fitness)
  rect <- function(aes_fill) {
    ggplot2::geom_rect(ggplot2::aes(xmin = .data[["xmin"]], xmax = .data[["xmax"]],
                                    ymin = .data[["ymin"]], ymax = .data[["ymax"]],
                                    fill = .data[[aes_fill]]))
  }
  if (fill == "sign") {
    data$invades <- data$fitness > 0
    p <- ggplot2::ggplot(data) + rect("invades") +
      ggplot2::scale_fill_manual(values = c(`TRUE` = "forestgreen", `FALSE` = "white"),
                                 na.value = "grey85", name = "mutant invades")
  } else {
    data$fitness[!is.finite(data$fitness)] <- NA_real_
    p <- ggplot2::ggplot(data) + rect("fitness") +
      ggplot2::scale_fill_gradient2(low = "firebrick", mid = "white", high = "forestgreen",
                                    midpoint = 0, na.value = "grey85", name = "fitness")
  }
  if (contours && nrow(x$contours) > 0L) {
    p <- p + ggplot2::geom_point(data = x$contours,
                                 ggplot2::aes(x = .data[["resident"]], y = .data[["mutant"]]),
                                 size = 0.6)
  }
  p + pip_axes(x$trait_scale, x$trait_name, "resident", "mutant")
}

##' Plot a trait-evolution plot.
##'
##' Each coexisting pair is a point at \eqn{(x_1, x_2)} with an arrow in the
##' direction of the selection gradient on the two residents, scaled so the
##' longest arrow spans \code{arrow} of the axis range on the plotting scale.
##' A dimorphic singular coalition sits where the arrows vanish.
##'
##' @title Plot a trait-evolution plot
##' @param x A \code{tep} from \code{\link{community_tep}}.
##' @param arrow Length of the longest arrow as a fraction of the axis range.
##' @param ... Ignored.
##' @return A \code{ggplot} object.
##' @author Daniel Falster
##' @export
plot.tep <- function(x, arrow = 0.05, ...) {
  tf <- pip_transform(attr(x, "trait_scale"))
  data <- tibble::as_tibble(unclass(x)[c("x1", "x2", "n1", "n2", "g1", "g2")])
  ## On a log axis the natural direction is the gradient with respect to
  ## log(x), which is x times the gradient.
  d1 <- if (identical(tf$scale, "log")) data$g1 * data$x1 else data$g1
  d2 <- if (identical(tf$scale, "log")) data$g2 * data$x2 else data$g2
  span <- diff(range(tf$fwd(c(data$x1, data$x2))))
  longest <- max(sqrt(d1^2 + d2^2), .Machine$double.eps)
  k <- arrow * span / longest
  data$xend <- tf$inv(tf$fwd(data$x1) + k * d1)
  data$yend <- tf$inv(tf$fwd(data$x2) + k * d2)

  ggplot2::ggplot(data, ggplot2::aes(x = .data[["x1"]], y = .data[["x2"]])) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data[["xend"]], yend = .data[["yend"]]),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.15, "cm")),
                          colour = "grey30") +
    ggplot2::geom_point(size = 0.8) +
    pip_axes(tf$scale, attr(x, "trait_name"), "resident 1", "resident 2")
}
