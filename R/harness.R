# Model harness (backend) abstraction for regnans
# ================================================
#
# A *harness* wires a `community` to a specific demographic model. The community
# pipeline (community_demography, community_selection_gradient,
# community_solve_singularity_1D, the assembler, ...) is model-agnostic: it only
# ever calls a small set of connector functions, each of which forwards to the
# harness's own implementation:
#
#   community_parameters()                    build model parameters
#   community_make_demography_runner()        closure: birth_rates -> offspring
#   community_demography_runner_cleanup()     write equilibrium state back
#   community_viable_bounds()                 viable trait region (empty community)
#   community_check_for_inviable_strategies() flag residents to drop
#   community_update_fitness_function()       build the invasion-fitness closure
#
# A harness is a list carrying a `$fns` table of these six functions plus any
# model state; its class is used for printing/dispatch. Two families ship:
#
#   harness_plant()    - the full `plant` SCM (the default), specialised by the
#                        physiological model (FF16 / TF24) and plant package
#                        version. Thin wrapper over the existing plant_community_*
#                        code.
#   harness_explicit() - fast toy models whose invasion fitness and equilibrium
#                        are supplied as explicit (C++) functions instead of
#                        emerging from the SCM. Used to develop and validate the
#                        assembly / attractor algorithms (issue #33). The shipped
#                        instances are named by author/year:
#                          harness_dd99()  Dieckmann & Doebeli 1999
#                          harness_gk98()  Geritz, Kisdi, Meszena & Metz 1998
#                          harness_gm99()  Geritz, van der Meijden & Metz 1999
#                          harness_jj12()  Johansson & Jonzen 2012 (bird arrival)
#                        ("explicit" is about the mechanism, not a guarantee that
#                        every quantity is closed-form: gm99's fitness is a
#                        Poisson series and its singular strategy is numerical.)

# ---- connectors: the only entry points the pipeline uses --------------------

community_parameters <- function(community) {
  community$harness$fns$parameters(community)
}
community_make_demography_runner <- function(community) {
  community$harness$fns$make_demography_runner(community)
}
community_demography_runner_cleanup <- function(community, runner, converged = TRUE) {
  community$harness$fns$demography_runner_cleanup(community, runner, converged)
}
community_viable_bounds <- function(community) {
  community$harness$fns$viable_bounds(community)
}
community_check_for_inviable_strategies <- function(community) {
  community$harness$fns$check_for_inviable_strategies(community)
}
community_update_fitness_function <- function(community) {
  community$harness$fns$update_fitness_function(community)
}

# ---- plant harness ----------------------------------------------------------

##' The default model harness: the full `plant` SCM.
##'
##' Forwards each connector to the existing `plant_community_*` implementation,
##' so the working plant path is unchanged. It is identified by the plant
##' physiological model (`FF16`, the current default, or `TF24`) and the plant
##' package `version`, mirroring how the toy harnesses are named by author/year.
##'
##' The actual `plant` `Parameters` are still supplied via
##' `community_start(model_support = ...)`; `model`/`version` here are recorded
##' metadata and a hook for model-specific dispatch (TF24 support depends on the
##' installed `plant`).
##'
##' @title plant model harness
##' @param model plant physiological model: "FF16" or "TF24"
##' @param version plant package version (defaults to the installed version)
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_plant <- function(model = c("FF16", "TF24"),
                          version = as.character(utils::packageVersion("plant"))) {
  model <- match.arg(model)
  h <- list(
    model   = model,
    version = version,
    fns = list(
      parameters                    = plant_community_parameters,
      make_demography_runner        = plant_community_make_demography_runner,
      demography_runner_cleanup     = plant_community_demography_runner_cleanup,
      viable_bounds                 = plant_community_viable_bounds,
      check_for_inviable_strategies = plant_community_check_for_inviable_strategies,
      update_fitness_function       = plant_community_update_fitness_function
    )
  )
  class(h) <- c(paste0("harness_plant_", tolower(model)), "harness_plant", "harness")
  h
}

# ---- generic explicit harness -----------------------------------------------

##' A harness for fast toy models with explicitly-supplied fitness/equilibrium.
##'
##' Implements all six connectors generically in terms of two model primitives:
##' a vectorised invasion-fitness function and a resident-equilibrium solver
##' (both typically backed by C++). Concrete instances --- \code{harness_dd99},
##' \code{harness_gk98}, \code{harness_gm99}, \code{harness_jj12} --- supply
##' these. "Explicit" refers to the mechanism (the fitness/equilibrium are
##' computed directly, not by running the plant SCM), NOT a claim that every
##' quantity is closed-form.
##'
##' Conventions:
##' \itemize{
##'   \item \code{community$birth_rate} stores resident abundance/density (the
##'     state the equilibrium solver returns), not a plant offspring rate.
##'   \item \code{fitness} returns invasion fitness with the resident value ~0
##'     (a log ratio for jj12/gk98/gm99; a per-capita rate for dd99).
##'   \item demography is solved by returning the model's equilibrium directly;
##'     the standard \code{equilibrium_iteration} loop then converges in a couple
##'     of steps.
##' }
##'
##' @title Explicit (toy-model) harness
##' @param fitness function(x_mut, x_res, n_res, pars) -> numeric vector of
##'   invasion fitness, one per mutant in \code{x_mut}.
##' @param equilibrium function(x_res, pars) -> numeric vector of resident
##'   equilibrium densities, one per resident.
##' @param pars named list of model parameters.
##' @param trait_names character vector naming the evolving trait(s).
##' @param label short model label (used in the harness class).
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_explicit <- function(fitness, equilibrium, pars, trait_names,
                             label = "explicit") {
  h <- list(
    pars        = pars,
    trait_names = trait_names,
    label       = label,
    # close pars into the model primitives so connectors call them with (x, ...)
    fitness     = function(x_mut, x_res, n_res) fitness(x_mut, x_res, n_res, pars),
    equilibrium = function(x_res) equilibrium(x_res, pars),
    fns = list(
      parameters                    = explicit_community_parameters,
      make_demography_runner        = explicit_community_make_demography_runner,
      demography_runner_cleanup     = explicit_community_demography_runner_cleanup,
      viable_bounds                 = explicit_community_viable_bounds,
      check_for_inviable_strategies = explicit_community_check_for_inviable_strategies,
      update_fitness_function       = explicit_community_update_fitness_function
    )
  )
  class(h) <- c(paste0("harness_", label), "harness_explicit", "harness")
  h
}

explicit_community_parameters <- function(community) {
  list(traits     = community$traits,
       birth_rate = community$birth_rate,
       pars       = community$harness$pars)
}

## Resident traits in the shape a model primitive expects: a plain vector for a
## one-trait model, the full trait matrix for a multi-trait (nD) model.
explicit_resident_traits <- function(community) {
  if (length(community$trait_names) == 1L) community$traits[, 1] else community$traits
}

explicit_community_make_demography_runner <- function(community) {
  h <- community$harness
  last_offspring_production <- NULL
  history <- list()
  function(birth_rates) {
    x_res <- explicit_resident_traits(community)
    out <- h$equilibrium(x_res)
    last_offspring_production <<- out
    history[[length(history) + 1L]] <<- list(`in` = birth_rates, out = out)
    out
  }
}

explicit_community_demography_runner_cleanup <- function(community, runner,
                                                         converged = TRUE) {
  e <- environment(runner)
  ## demography_solve_equilibrium_solve() wraps the runner when it excludes
  ## extinct species; the model state lives in the wrapped runner (mirrors
  ## plant_community_demography_runner_cleanup()).
  if (is.function(e$runner_full)) {
    runner <- e$runner_full
    e <- environment(runner)
  }
  community$birth_rate <- e$last_offspring_production
  community$fitness_points <- NULL
  attr(community, "converged") <- converged
  attr(community, "progress") <- e$history
  community
}

explicit_community_update_fitness_function <- function(community) {
  h <- community$harness
  k <- length(community$trait_names)
  x_res <- explicit_resident_traits(community)
  n_res <- community$birth_rate

  community$resident_fitness <-
    if (nrow(community$traits) > 0L) h$fitness(x_res, x_res, n_res) else numeric()

  community$fitness_function <- function(x) {
    if (k == 1L) {
      if (is.matrix(x)) x <- x[, 1]
      h$fitness(as.numeric(x), x_res, n_res)
    } else {
      # mutants as rows of a k-column matrix (a bare vector is one mutant)
      xm <- if (is.matrix(x)) x else matrix(x, nrow = 1L)
      h$fitness(xm, x_res, n_res)
    }
  }
  community
}

explicit_community_viable_bounds <- function(community) {
  # Toy models use the bounds supplied to community_start(); we do not bracket a
  # fundamental-fitness region (these models are typically viable over the whole
  # supplied range). Bounds are already on the (empty) community.
  community
}

explicit_community_check_for_inviable_strategies <- function(community) {
  runner <- community_make_demography_runner(community)
  op <- runner(community$birth_rate)
  eps <- community$demography_control$equilibrium_extinct_birth_rate
  drop <- op < eps
  attr(op, "drop") <- drop
  op
}

# ---- concrete toy models ----------------------------------------------------

##' DD99: Dieckmann & Doebeli 1999 competition model (Nature 400:354-357).
##'
##' Continuous-time logistic competition for a Gaussian resource: carrying
##' capacity \code{K(x) = K0 exp(-(x-x0)^2/(2 sigma_K^2))} and competition kernel
##' \code{C(d) = exp(-d^2/(2 sigma_C^2))}. Invasion fitness
##' \code{s(y) = r (1 - sum_i N_i C(y-x_i)/K(y))}. The singular strategy is
##' \code{x* = x0}; it is an evolutionary \emph{branching point} iff
##' \code{sigma_C < sigma_K} and an ESS iff \code{sigma_C > sigma_K} -- the
##' classic disruptive-selection oracle.
##'
##' @title DD99 (Dieckmann & Doebeli 1999) harness
##' @param r intrinsic growth rate
##' @param K0 maximum carrying capacity
##' @param x0 trait optimising the resource (singular strategy)
##' @param sigma_K width of the resource/carrying-capacity kernel
##' @param sigma_C width of the competition kernel
##' @param trait_name name of the evolving trait
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_dd99 <- function(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4,
                         trait_name = "x") {
  pars <- list(r = r, K0 = K0, x0 = x0, sigma_K = sigma_K, sigma_C = sigma_C)
  harness_explicit(
    fitness     = dd99_fitness,
    equilibrium = function(x_res, pars) dd99_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_name,
    label       = "dd99"
  )
}

##' DD99 in two (or more) trait dimensions (cf. Ito & Dieckmann 2007).
##'
##' The multi-trait extension of \code{\link{harness_dd99}}: the resource and
##' competition kernels are products of per-dimension Gaussians, so each trait
##' dimension `d` has its own optimum \code{x0[d]} and widths
##' \code{sigma_K[d]}, \code{sigma_C[d]}. The singular strategy is \code{x* = x0}
##' and dimension `d` is a branching direction iff
##' \code{sigma_C[d] < sigma_K[d]}. Use \code{trait_scale = "linear"} in
##' \code{community_start()} since the optima sit at 0.
##'
##' @title DD99 multi-trait harness
##' @param r intrinsic growth rate
##' @param K0 maximum carrying capacity
##' @param x0 length-k vector of per-dimension optima (the singular strategy)
##' @param sigma_K length-k vector of resource-kernel widths
##' @param sigma_C length-k vector of competition-kernel widths
##' @param trait_names length-k character vector naming the traits
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_dd99_nd <- function(r = 1, K0 = 500, x0 = c(0, 0),
                            sigma_K = c(1, 1), sigma_C = c(0.4, 0.4),
                            trait_names = c("x1", "x2")) {
  k <- length(trait_names)
  if (length(x0) != k || length(sigma_K) != k || length(sigma_C) != k) {
    stop("x0, sigma_K, sigma_C must each have length length(trait_names)")
  }
  pars <- list(r = r, K0 = K0, x0 = x0, sigma_K = sigma_K, sigma_C = sigma_C)
  harness_explicit(
    fitness     = dd99_nd_fitness,
    equilibrium = function(x_res, pars) dd99_nd_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_names,
    label       = "dd99_nd"
  )
}

##' GK98: Geritz, Kisdi, Meszena & Metz 1998 soft-selection model (Evol. Ecol.
##' 12:35-57; the worked Levene example).
##'
##' \code{m} patches with optima \code{mu} and capacities \code{K}; Gaussian
##' within-patch survival of width \code{sigma}. Invasion fitness reduces, for a
##' single resident, to \code{S(y) = log sum_j c_j f_j(y)/f_j(x)} with
##' \code{c_j = K_j/sum K}. For the symmetric three-patch default
##' (\code{mu = (-d, 0, d)}, equal \code{K}) the singular strategy is
##' \code{x* = 0}: convergence stable for all \code{d/sigma}, and an evolutionary
##' branching point iff \code{d/sigma > sqrt(3/2) ~= 1.2247} (a CSS below that).
##'
##' @title GK98 (Geritz et al. 1998 soft-selection) harness
##' @param d patch-optimum spacing for the symmetric default (optima -d, 0, d)
##' @param sigma within-patch survival width
##' @param mu patch optima (overrides the symmetric default built from \code{d})
##' @param K patch capacities (defaults to equal capacities)
##' @param trait_name name of the evolving trait
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_gk98 <- function(d = 1.5, sigma = 1, mu = c(-d, 0, d),
                         K = rep(1, length(mu)), trait_name = "x") {
  pars <- list(sigma = sigma, mu = mu, K = K)
  harness_explicit(
    fitness     = gk98_fitness,
    equilibrium = function(x_res, pars) gk98_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_name,
    label       = "gk98"
  )
}

##' GM99: Geritz, van der Meijden & Metz 1999 seed-size safe-site model.
##'
##' Plants compete for safe sites; seeds arrive by a Poisson process and
##' seedlings undergo size-asymmetric lottery competition. Pre-competitive
##' survival \code{s(x) = max(0, 1 - 2 exp(-beta x))}, competitive ability
##' \code{c(x) = exp(alpha x)}, invasion fitness (a lifetime reproductive ratio)
##' \code{W = (R/x) s(x) g(x_mut, x, N)} (mutant trait \code{x_mut} against
##' resident \code{x}) averaged over the Poisson number of
##' competitors. Only the products \code{alpha*R} and \code{beta*R} matter. There
##' is no closed-form singular strategy; size-asymmetric competition
##' (large \code{alpha}) drives evolutionary branching in seed size (see
##' Fig. 5 of the paper: e.g. \code{alpha*R = 4.5} vs \code{7.0} at \code{beta*R = 15}).
##'
##' This is the Geritz \emph{1999} seed-size model (the one in the original MATLAB code),
##' distinct from the GK98 1998 soft-selection model (\code{\link{harness_gk98}}).
##'
##' @title GM99 (Geritz et al. 1999 seed-size) harness
##' @param alpha competitive asymmetry (larger = more size-asymmetric)
##' @param beta habitat parameter controlling pre-competitive survival
##' @param R resources per safe site (only alpha*R and beta*R matter)
##' @param trait_name name of the evolving trait (seed size)
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_gm99 <- function(alpha = 6, beta = 25, R = 1, trait_name = "x") {
  pars <- list(R = R, alpha = alpha, beta = beta)
  harness_explicit(
    fitness     = gm99_fitness,
    equilibrium = function(x_res, pars) gm99_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_name,
    label       = "gm99"
  )
}

##' JJ12: migratory-bird arrival-time model (Johansson & Jonzen 2012; the
##' simplified analytic form of Brannstrom, Johansson & von Festenberg 2013,
##' Games 4:304-328, section 4).
##'
##' Trait = arrival time. Early arrival raises competitive ability
##' \code{C(x) = exp(-a x)} for a fixed number \code{K} of territories;
##' reproduction \code{R(x) = R0 exp(-(x - x_opt)^2 / (2 sigma^2))} peaks at the
##' seasonal optimum; survival \code{p}. The model has a single continuously
##' stable strategy (no branching) with the closed form
##' \deqn{x^* = x_{opt} - a\,\sigma^2,}
##' which arrives \emph{before} the population optimum (a tragedy of the commons).
##' This exact answer makes the model a benchmark for numerical
##' selection-gradient / singular-strategy solvers.
##'
##' @title JJ12 (Johansson & Jonzen 2012 bird arrival) harness
##' @param a competitive advantage of early arrival (>= 0)
##' @param x_opt seasonal optimum arrival time
##' @param sigma width of the benign season
##' @param R0 maximum reproductive output
##' @param K number of territories
##' @param p year-to-year survival, in (0, 1)
##' @param trait_name name of the evolving trait
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_jj12 <- function(a = 0.1, x_opt = 0, sigma = 1, R0 = 1, K = 1, p = 0.5,
                         trait_name = "x") {
  pars <- list(a = a, x_opt = x_opt, sigma = sigma, R0 = R0, K = K, p = p)
  harness_explicit(
    fitness     = jj12_fitness,
    equilibrium = function(x_res, pars) jj12_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_name,
    label       = "jj12"
  )
}

##' @export
print.harness <- function(x, ...) {
  cat(sprintf("<harness: %s>\n", paste(class(x)[-length(class(x))], collapse = ", ")))
  if (inherits(x, "harness_plant")) {
    cat(sprintf("  plant model: %s (plant %s)\n", x$model, x$version))
  }
  if (!is.null(x$pars)) {
    flat <- vapply(x$pars, function(v) paste(format(v), collapse = ","), character(1))
    cat("  parameters:", paste(names(flat), flat, sep = "=", collapse = ", "), "\n")
  }
  invisible(x)
}
