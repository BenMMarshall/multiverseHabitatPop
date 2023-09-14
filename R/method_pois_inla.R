#' Calculate population selection using a Poisson model
#'
#' @name method_pois_inla
#' @description A
#' @param availUseData comes from area_based_extraction
#' @param sampleGroups list if IDs for the sample groups
#' @param optionsList options list that will be used to loop through options
#' @return Population estimates for area based methods.
#'
#' @export
method_pois_inla <- function(availUseData, sampleGroups, optionsList){
  
  optionsForm <- optionsList$MethodPois_mf
  optionsApoints <- optionsList$MethodPois_as
  
  # listLength <- length(optionsType)*
  #   length(optionsMethod)*
  #   length(optionsContour)*
  #   length(optionsAPoints)*
  #   length(optionsSPSamp)*
  #   length(optionsTest)
  # 
  # companaResultsList <- vector("list", length = listLength)
  # i <- 0
  # for(typ in optionsType){
  #   for(met in optionsMethod){
  #     for(con in optionsContour){
  #       for(aPo in optionsAPoints){
  #         for(spS in optionsSPSamp){
  #           
  #           for(sampID in names(sampleGroups)){
  #             
  #             IDs <- optionsList_samples[[sampID]]
  #             
  #             IDs <- paste0("simData_i", sprintf("%03d", IDs))
  #             
  #             use <- availUseData %>% 
  #               filter(type == typ,
  #                      method == met,
  #                      contour == con,
  #                      aPoints == aPo,
  #                      spSamp == spS) %>% 
  #               filter(id %in% IDs) %>% 
  #               select(used_c0, used_c2)
  #             
  #             avail <- availUseData %>% 
  #               filter(type == typ,
  #                      method == met,
  #                      contour == con,
  #                      aPoints == aPo,
  #                      spSamp == spS) %>% 
  #               filter(id %in% IDs) %>% 
  #               select(avail_c0, avail_c2)
  #             
  #             names(use) <- c("c0", "c2")
  #             names(avail) <- c("c0", "c2")
  #             
  #             for(tes in optionsTest){
  #               
  #               companaOUT <- compana(used = use, avail = avail,
  #                                     test = tes)
  #               
  #               companaResultsDF <- data.frame(
  #                 sampleID = sampID,
  #                 sampleSize = length(IDs),
  #                 type = typ,
  #                 areaMethod = met,
  #                 contour = con,
  #                 aPoints = aPo,
  #                 spSamp = spS,
  #                 test = tes,
  #                 companaLambda = companaOUT$test["Lambda"],
  #                 companaP = companaOUT$test["P"]
  #               )
  #               
  #               i <- i+1
  #               
  #               companaResultsList[[i]] <- companaResultsDF
  #               
  #             } # tes
  #             
  #           } # samp
  #           
  #         }
  #       }
  #     }
  #   }
  # }
  # 
  # return(do.call(rbind, companaResultsList))
  
}