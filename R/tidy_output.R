#' Tidies output from assembler_run into tidy dataframe.
#'
#' @param obj output of an assembly run
#'
#' @return A tibble of the assembly, one row per resident per step. The
#'   community's trait scale is kept in \code{attr(, "trait_scale")} for
#'   \code{\link{plot_community}}.
#' @export
#'
tidy_assembly <- function(obj){
  out <- obj$history %>%
    purrr::map(tidy_community) %>%
    dplyr::bind_rows(.id = "step") %>%
    mutate(strategy_id = strategy_id %>% as.factor() %>% as.character())
  attr(out, "trait_scale") <- community_trait_transform(obj$community)$scale
  out
}

tidy_community <- function(community){
  
  out <- community$traits %>% dplyr::as_tibble()
  
  traits <- names(out)

  out %>%
    dplyr::mutate(
      resident = TRUE,
      births = community$birth_rate
    ) %>%
    tidyr::unite("strategy_id", any_of(traits), remove = FALSE) %>%
    tidyr::nest(data = dplyr::any_of(traits)) %>%
    dplyr::rename("traits" = data) %>%
    dplyr::mutate(fitness_landscape = list(community$fitness_points)) %>%
    relocate(strategy_id, .before = resident)
}

