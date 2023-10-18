#' Render the manuscript rmd
#'
#' @name render_rmd
#' @description A
#' @param modelExtractsTarget table of R2 and effects for quick reference in the rmd
#' @param effectPlotsTarget ggplot objects for plots showing the effects for quick reference in the rmd
#' @return Nothing, PDF (or output) will be saved to a folder.
#'
#' @export
render_rmd <- function(modelExtracts, effectPlots,
                       areaSpecCurve,
                       twoStepSpecCurve,
                       poisSpecCurve,
                       ssfSpecCurve){
  
  rmarkdown::render(input = here::here("notebook", "manuscript",
                                       "multiverseHabitatPopulationManuscript.Rmd"),
                    output_file = here::here("notebook", "manuscript",
                                             "multiverseHabitatPopulationManuscript.pdf"))
}
