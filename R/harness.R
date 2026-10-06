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
# model state; its class is used for printing/dispatch. A harness that can
# differentiate its own fitness says so in `$provides` and sets the matching
# closures in community$fitness_derivatives when it builds the fitness
# function; harness_fd() supplies finite differences for the rest (see
# R/derivatives.R). Two families ship:
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

##' Which derivatives a harness supplies itself.
##'
##' @title Query a harness's derivatives
##' @param harness A \code{harness}.
##' @param what A derivative name: \code{"fitness_gradient"} or
##' \code{"fitness_hessian"}. With \code{NULL}, the names of all it provides.
##' @return A logical, or with \code{what = NULL} a character vector.
##' @author Daniel Falster
##' @export
harness_provides <- function(harness, what = NULL) {
  provides <- harness$provides
  if (is.null(provides)) provides <- character(0)
  if (is.null(what)) {
    return(provides)
  }
  what <- match.arg(what, fitness_derivative_names)
  what %in% provides
}

##' Reach an explicit harness's equilibrium by iterating its demography.
##'
##' The reference models return their resident equilibrium in closed form (or
##' from their own internal solver), so the package's equilibrium solvers never
##' iterate on them and nothing about warm starts, seeding or convergence is
##' exercised. This wrapper replaces that shortcut with the model's own
##' generation-to-generation demography,
##' \deqn{n_{t+1} = n_t \exp(s(x; x, n_t)),}
##' where \eqn{s} is the model's invasion fitness of each resident against the
##' current community --- for GK98 this is its soft-selection recursion (Eq. 19),
##' for JJ12 the territory recursion, for GM99 the seeds-per-seed return; for
##' DD99, whose fitness is a per-capita rate, it is a unit time step.
##'
##' The wrapper supplies only that one-step map, as the demography runner,
##' because that is all a model is asked to know: what its residents produce
##' next generation. Finding the fixed point of the map is an algorithm, and
##' it belongs to regnans, chosen per community, so that the same solver code
##' runs against the plant SCM and against a reference model --- a solver
##' verified here is the one plant gets. \emph{Which} solver iterates the map
##' to the fixed point is therefore the community's choice, exactly as for the
##' plant model:
##' \code{community$demography_control$equilibrium_solver_name}, set through
##' \code{\link{demographic_step_control}} --- plain fixed-point iteration
##' (\code{"equilibrium_iteration"}, the default), root-finding on
##' \eqn{n - f(n)} with \code{"equilibrium_solve_nleqslv"} or
##' \code{"equilibrium_solve_dfsane"}, or the \code{"equilibrium_hybrid"}
##' combination --- together with its tolerance and step limit. The starting
##' density is the resident's \code{birth_rate}. Use it to test and time the
##' solvers on models whose answers are known.
##'
##' @title Iterate an explicit harness's demography to equilibrium
##' @param harness A harness from \code{\link{harness_explicit}} or one of the
##' shipped reference models.
##' @return The harness, with its demography runner replaced by the model's
##' one-generation map.
##' @examples
##' h <- harness_iterate_demography(harness_gm99(alpha = 7, beta = 15))
##' comm <- community_start(bounds(x = c(0.08, 0.9)), trait_scale = "log", harness = h,
##'                         demography_control = demographic_step_control(
##'                           list(equilibrium_solver_name = "equilibrium_solve_nleqslv",
##'                                equilibrium_eps = 1e-8))) |>
##'   community_add(trait_matrix(0.3, "x"), birth_rate = 1) |>
##'   community_demography()
##' comm$birth_rate                       # the same equilibrium the closed form gives
##' NROW(attr(comm, "progress"))          # how many one-generation steps it cost
##' @author Daniel Falster
##' @export
harness_iterate_demography <- function(harness) {
  if (!inherits(harness, "harness_explicit")) {
    stop("harness_iterate_demography() applies to explicit (reference-model) harnesses")
  }
  harness$mode <- "iterate_demography"
  harness$fns$make_demography_runner <- explicit_community_make_dynamics_runner
  harness
}

## The dynamics map n -> n exp(s(x; x, n)) as a demography runner, recording
## the same state the equilibrium runner does so the cleanup is shared.
explicit_community_make_dynamics_runner <- function(community) {
  h <- community$harness
  x_res <- explicit_resident_traits(community)
  last_offspring_production <- NULL
  history <- list()
  function(birth_rates) {
    n <- as.numeric(birth_rates)
    s <- as.numeric(h$fitness(x_res, x_res, n))
    out <- n * exp(s)
    out[!is.finite(out)] <- 0
    last_offspring_production <<- out
    history[[length(history) + 1L]] <<- list(`in` = birth_rates, out = out)
    out
  }
}

##' Supply finite-difference derivatives for whatever a harness lacks.
##'
##' The methods in this package are written against exact derivatives of
##' invasion fitness (\code{\link{community_fitness_gradient}},
##' \code{\link{community_fitness_hessian}}). A model that cannot differentiate
##' itself still works, at the cost and accuracy of finite differences, once its
##' harness has been through this function: it wraps the harness's
##' fitness-function connector so that any derivative not provided by the model
##' is computed by the finite-difference stencils in \code{R/util_gradient.R},
##' with the settings in \code{\link{derivative_control}}. Derivatives the
##' model does provide are left alone, and the source of each is recorded in
##' \code{community$fitness_derivatives$source}.
##'
##' \code{\link{community_start}} applies this to every harness, so it is only
##' needed when building a community by hand. Applying it twice is harmless.
##'
##' @title Finite-difference derivatives for a harness
##' @param harness A \code{harness}.
##' @return The harness, with its fitness-function connector wrapped.
##' @author Daniel Falster
##' @export
harness_fd <- function(harness) {
  if (!inherits(harness, "harness")) {
    stop("harness_fd() needs a harness object")
  }
  if (isTRUE(harness$fd)) {
    return(harness)
  }
  model_update <- harness$fns$update_fitness_function
  harness$fns$update_fitness_function <- function(community) {
    fd_fill_fitness_derivatives(model_update(community))
  }
  harness$fd <- TRUE
  harness
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
##' (both typically backed by C++), plus optional derivatives of fitness in the
##' mutant direction. Concrete instances --- \code{harness_dd99},
##' \code{harness_gk98}, \code{harness_gm99}, \code{harness_jj12} --- supply
##' all of them. "Explicit" refers to the mechanism (the fitness/equilibrium are
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
##' @param fitness_gradient optional function(x_mut, x_res, n_res, pars) -> the
##'   gradient of invasion fitness with respect to the mutant trait(s), one row
##'   per mutant and one column per trait.
##' @param fitness_hessian optional function(x_mut, x_res, n_res, pars) -> the
##'   k by k Hessian of invasion fitness with respect to the mutant trait(s)
##'   for a single mutant.
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_explicit <- function(fitness, equilibrium, pars, trait_names,
                             label = "explicit",
                             fitness_gradient = NULL, fitness_hessian = NULL) {
  # close pars into the model primitives so connectors call them with (x, ...)
  with_pars <- function(f) {
    if (is.null(f)) NULL else function(x_mut, x_res, n_res) f(x_mut, x_res, n_res, pars)
  }
  h <- list(
    pars        = pars,
    trait_names = trait_names,
    label       = label,
    fitness     = with_pars(fitness),
    equilibrium = function(x_res) equilibrium(x_res, pars),
    fitness_gradient = with_pars(fitness_gradient),
    fitness_hessian  = with_pars(fitness_hessian),
    provides = fitness_derivative_names[c(!is.null(fitness_gradient),
                                          !is.null(fitness_hessian))],
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

  ## Mutants in the shape the model primitives expect: a plain vector for a
  ## one-trait model, rows of a k-column matrix otherwise (a bare vector is
  ## one mutant).
  model_points <- function(x) {
    if (k == 1L) {
      if (is.matrix(x)) x <- x[, 1]
      as.numeric(x)
    } else if (is.matrix(x)) {
      x
    } else {
      matrix(x, nrow = 1L)
    }
  }

  community$fitness_function <- function(x) h$fitness(model_points(x), x_res, n_res)

  der <- list()
  if (is.function(h$fitness_gradient)) {
    der$gradient <- function(y, control) h$fitness_gradient(model_points(y), x_res, n_res)
  }
  if (is.function(h$fitness_hessian)) {
    der$hessian <- function(y, control) h$fitness_hessian(model_points(y), x_res, n_res)
  }
  community$fitness_derivatives <- der
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
    label       = "dd99",
    fitness_gradient = dd99_fitness_gradient,
    fitness_hessian  = dd99_fitness_hessian
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
    label       = "dd99_nd",
    fitness_gradient = dd99_nd_fitness_gradient,
    fitness_hessian  = dd99_nd_fitness_hessian
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
    label       = "gk98",
    fitness_gradient = gk98_fitness_gradient,
    fitness_hessian  = gk98_fitness_hessian
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
    label       = "gm99",
    fitness_gradient = gm99_fitness_gradient,
    fitness_hessian  = gm99_fitness_hessian
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
    label       = "jj12",
    fitness_gradient = jj12_fitness_gradient,
    fitness_hessian  = jj12_fitness_hessian
  )
}

##' @export
print.harness <- function(x, ...) {
  cat(sprintf("<harness: %s>\n", paste(class(x)[-length(class(x))], collapse = ", ")))
  if (inherits(x, "harness_plant")) {
    cat(sprintf("  plant model: %s (plant %s)\n", x$model, x$version))
  }
  if (identical(x$mode, "iterate_demography")) {
    cat("  equilibrium: iterated from the model's own demography\n")
  }
  if (!is.null(x$pars)) {
    flat <- vapply(x$pars, function(v) paste(format(v), collapse = ","), character(1))
    cat("  parameters:", paste(names(flat), flat, sep = "=", collapse = ", "), "\n")
  }
  provides <- harness_provides(x)
  cat("  derivatives:",
      if (length(provides) > 0L) paste(provides, collapse = ", ") else "none",
      if (isTRUE(x$fd)) "(finite differences for the rest)" else "", "\n")
  invisible(x)
}
