#' Selection gradient of the residents
#'
#' The derivative of invasion fitness with respect to the mutant trait,
#' evaluated at each resident's own traits: the direction in which selection
#' pushes that resident. Solve the community to demographic equilibrium first
#' (\code{\link{community_demography}}). Derivative settings come from
#' \code{\link{derivative_control}}.
#' @param community A \code{community} solved to demographic equilibrium.
#' @author Daniel Falster
#' @export
#' @return The community with \code{selection_gradient} (a vector for a single
#' resident, otherwise a residents-by-traits matrix) and \code{resident_fitness}
#' set.
community_selection_gradient <- function(community) {

  msg <- sprintf("Calculating selection gradient for [%s] = [%s]",
    paste(community$trait_names, collapse = ", "), 
    paste(community$traits, collapse = ", ")
  )
  plant_log_assembler(msg)

  g <- community_fitness_gradient(community)
  ret <- if (nrow(g) == 1L) as.vector(g[1, ]) else g

  msg <- sprintf("Solved! Selection gradient for [%s] = [%s] is [%s]", 
    paste(community$trait_names, collapse = ", "), 
    paste(community$traits, collapse = ", "),
    paste(ret, collapse = ", ")
  )
  plant_log_assembler(msg)

  community[["resident_fitness"]] <- attr(g, "value")
  community[["selection_gradient"]] <- ret
  community
}
