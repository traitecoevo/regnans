##' Controls how invasion-fitness landscapes are sampled.
##'
##' Returns a list of options read by \code{\link{community_fitness_landscape}}
##' and its methods. Passing a list via \code{control} overrides the defaults;
##' unknown names are an error. This is the \code{fitness_control} slot of a
##' \code{community} (see \code{\link{community_start}}), and is the fitness
##' analogue of \code{\link{demographic_step_control}}.
##'
##' Options:
##' \itemize{
##'   \item \code{method} --- \code{"grid"} (the default: evenly spaced points on
##'     the community's trait scale, augmented with the residents) or
##'     \code{"bayesopt"} (Gaussian-process surrogate via \pkg{mlr3mbo}).
##'   \item \code{n_evals} --- number of fitness evaluations.
##'   \item \code{n_init} --- size of the initial design for \code{"bayesopt"}.
##' }
##'
##' @title Options controlling fitness-landscape construction
##' @param control A list of values to modify from the defaults.
##' @return A list with elements \code{method}, \code{n_evals}, \code{n_init}.
##' @author Daniel Falster
##' @export
fitness_landscape_control <- function(control = NULL) {
  defaults <- list(
    method  = "grid",
    n_evals = 51,
    n_init  = 20
  )

  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown fitness control parameters ", paste(extra, collapse = ", "))
  }
  modifyList(defaults, control)
}

## Fitness-landscape options for a community, tolerating communities built
## before fitness_control had a default (or by hand, without going through
## community_start()).
community_fitness_control <- function(community) {
  fitness_landscape_control(community$fitness_control)
}

##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param community A community object
##' @param method used to construct landscape; defaults to
##' \code{community$fitness_control$method} and falls back to \code{"grid"}
##' (see \code{\link{fitness_landscape_control}}).
##' @param ... additional arguments passed to the chosen landscape method
##' (e.g. \code{bounds}, \code{n_evals}).
##' @author Daniel Falster
##' @rdname community_fitness_landscape
##' @export
community_fitness_landscape <- function(community, method = NULL, ...) {

  ## fitness_control may be NULL on a hand-built community; normalise so
  ## method/n_evals/n_init always resolve.
  community$fitness_control <- community_fitness_control(community)
  if (is.null(method)) {
    method <- community$fitness_control$method
  }

  if (is.null(community$fitness_function)) {
    community <- community %>% community_demography()
  }
  plant_log_assembler(sprintf(
    "Calulcating fitness landscape for %d strategy communtiy using %s", nrow(community$traits), method))  
  
  if(method == "grid") {
    community_fitness_landscape_grid(community, ...)
  } else if(method == "bayesopt") {
    community_fitness_landscape_bayesopt(community, ...)
  } else {
    stop("Unknown fitness landscape method")
  }
}

community_fitness_landscape_grid <- function(community, bounds = community$bounds, n_evals = community$fitness_control$n_evals) {

  ## Space the grid on the community's trait scale (log for positive plant
  ## traits, linear for traits spanning zero). seq_log_range() is the log case.
  tf <- community_trait_transform(community)
  lo <- tf$fwd(bounds[1, 1])
  hi <- tf$fwd(bounds[1, 2])
  x <- tf$inv(seq(lo, hi, length.out = n_evals))

  # add residents - also points offset from resident to capture local gradient
  res <- community$traits[, 1]
  if (length(res) > 0L) {
    eps <- 0.005 * (hi - lo)
    x <- c(x, tf$inv(tf$fwd(res) - eps), res, tf$inv(tf$fwd(res) + eps))
  }
  x <- sort(unique(x))

  y <- community$fitness_function(x)

  community$fitness_points <-
    dplyr::tibble(x = x, fitness = y) %>%
    dplyr::mutate(resident = ifelse(x %in% community$traits, TRUE, FALSE))

  names(community$fitness_points)[1] <- community$trait_names

  community
}

# Construct a fitness landscape using Bayesian optimisation.
# Uses the mlr3mbo package. Internal helper (not exported, no Rd).

community_fitness_landscape_bayesopt <- function(community, bounds = community$bounds, n_evals = community$fitness_control$n_evals, n_init = community$fitness_control$n_init) {
 
  set.seed(1)
  
  obfun <- bbotk::ObjectiveRFun$new(
    fun = function(xs) list(community$fitness_function(exp(xs$x))),
    domain = paradox::ps(x = paradox::p_dbl(lower = log(community$bounds[1]), upper = log(community$bounds[2]))),
    codomain = paradox::ps(y = paradox::p_dbl(tags = "maximize"))
  )

  optimizer <- bbotk::opt("mbo",
    loop_function = mlr3mbo::bayesopt_ego,
    surrogate = fitness_surrogate_start(),
    acq_function = mlr3mbo::acqf("ei"),
    acq_optimizer = mlr3mbo::acqo(
      bbotk::opt("nloptr", algorithm = "NLOPT_GN_ORIG_DIRECT"),
      terminator = bbotk::trm("stagnation", iters = 100, threshold = 1e-5)
    )
  )

  instance <- bbotk::OptimInstanceSingleCrit$new(
    objective = obfun,
    terminator = bbotk::trm("evals", n_evals = n_evals)
  )

  # Initial data -- 
  # space n_evals and add residents
  x <- sort(unique(c(seq_log_range(bounds, min(n_init, n_evals)), 
    0.995 * community$traits[, 1], 
    community$traits[, 1],
    1.005*community$traits[,1])))

  initial_design <- data.table::data.table(x = log(x))
  instance$eval_batch(initial_design)

  # run optimisation
  optimizer$optimize(instance)

  # Store points
  community$fitness_surrogate_archive <- instance$archive

  community$fitness_points <-
    instance$archive$data %>% 
    dplyr::as_tibble() %>%
    dplyr::mutate(x = exp(x)) %>%  #back transform x
    dplyr::select(x = x, fitness = y, batch_nr) %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(resident = ifelse(x %in% community$traits, TRUE, FALSE))
  names(community$fitness_points)[1] <- community$trait_names

  # Store surrogate
  community <- community_fitness_surrogate_create(community)

  #community$fitness_surrogate_function(0.01)
  #community$fitness_function(0.01)
  
  # # make predictions
  # xdt <- data.table::data.table(x = seq_range(log(community$bounds), length.out = 101))
  # surrogate_pred <- surrogate$predict(xdt)
  # pred <- bind_cols(xdt %>% as_tibble(), surrogate_pred %>% as_tibble()) %>%
  #   rename(y = mean)

  # # plot true function in black
  # # surrogate prediction (mean +- se in grey)
  # # known optimum in darkred
  # # found optimum in darkgreen
  # ggplot(aes(x, y), data = pred) +
  #   geom_point(data = instance$archive$data %>% as_tibble()) +
  #   geom_line(col = "red") +
  #   geom_ribbon(aes(x = x, ymin = y - se, ymax = y + se), fill = "grey", alpha = 0.2) +
  #   theme_minimal()
  community
}

fitness_surrogate_start <- function(archive = NULL) {
  requireNamespace("mlr3learners", quietly = TRUE) # registers the "regr.km" learner
  mlr3mbo::srlrn(mlr3::lrn("regr.km", covtype = "matern3_2", control = list(trace = FALSE)), archive = archive)
}

community_fitness_surrogate_create <- function(community, 
  archive = community$fitness_surrogate_archive) {

  requireNamespace("mlr3learners", quietly = TRUE) # registers the "regr.km" learner
  
  community$fitness_surrogate_object <- fitness_surrogate_start(archive)

  if(!is.null(archive))
    community$fitness_surrogate_object$update()
    
  community$fitness_surrogate_function <- function(x, se = FALSE) {
    xdt <- data.table::data.table(x = log(x))

    surrogate_pred <-
      community$fitness_surrogate_object$predict(xdt) %>%
      dplyr::as_tibble() %>%
      dplyr::rename(y = mean)

    if (!se) {
      surrogate_pred <- surrogate_pred[["y"]]
    }

    surrogate_pred
  }
  
  community
}
