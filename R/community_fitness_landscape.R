##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param community A community object
##' @param method used to construct landscape
##' @param ... additional arguments passed to the chosen landscape method
##' (e.g. \code{bounds}, \code{n_evals}).
##' @author Daniel Falster
##' @rdname community_fitness_landscape
##' @export
community_fitness_landscape <- function(community, method = community$fitness_control$method, ...) {

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
