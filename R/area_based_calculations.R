#' Calculate the population results for the area based methods
#'
#' @name area_based_calculations
#' @description A
#' @param availUseData comes from area_based_extraction
#' @param sampleGroups list if IDs for the sample groups
#' @param optionsList options list that will be used to loop through options
#' @return Population estimates for area based methods.
#'
#' @export
area_based_calculations <- function(availUseData, sampleGroups, optionsList){
  
  optionsType <- unique(availUseData$type)
  optionsMethod <- unique(availUseData$method)
  optionsContour <- unique(availUseData$contour)
  optionsAPoints <- unique(availUseData$availablePoints)
  optionsSPSamp <- unique(availUseData$samplingPattern)
  
  # optionsList$areaBasedMethod
  optionsTest <- optionsList$areaBasedTest
  
  listLength <- length(optionsType)*
    length(optionsMethod)*
    length(optionsContour)*
    length(optionsAPoints)*
    length(optionsSPSamp)*
    length(optionsTest) *
    length(names(sampleGroups))
  
  companaResultsList <- vector("list", length = listLength)
  i <- 0
  for(typ in optionsType){
    for(met in optionsMethod){
      for(con in optionsContour){
        for(aPo in optionsAPoints){
          for(spS in optionsSPSamp){
            
            for(sampID in names(sampleGroups)){
              
              IDs <- optionsList_samples[[sampID]]

              IDs <- paste0(stringr::str_extract(availUseData$id[1], "^.*_i"),
                            sprintf("%03d", IDs))
              
              use <- availUseData %>% 
                dplyr::filter(type == typ,
                       method == met,
                       contour == con,
                       availablePoints == aPo,
                       samplingPattern == spS) %>% 
                dplyr::filter(id %in% IDs) %>% 
                dplyr::select(used_c0, used_c2)
              
              avail <- availUseData %>% 
                dplyr::filter(type == typ,
                       method == met,
                       contour == con,
                       availablePoints == aPo,
                       samplingPattern == spS) %>% 
                dplyr::filter(id %in% IDs) %>% 
                dplyr::select(avail_c0, avail_c2)
              
              names(use) <- c("c0", "c2")
              names(avail) <- c("c0", "c2")
              
              for(tes in optionsTest){
                
                companaOUT <- compana(used = use, avail = avail,
                                      test = tes)
                
                companaResultsDF <- data.frame(
                  sampleID = sampID,
                  sampleSize = length(IDs),
                  trackFreq = availUseData$trackFreq[1],
                  trackDura = availUseData$trackDura[1],
                  type = typ,
                  areaMethod = met,
                  contour = con,
                  availablePoints = aPo,
                  samplingPattern = spS,
                  test = tes,
                  companaLambda = companaOUT$test["Lambda"],
                  companaP = companaOUT$test["P"]
                )
                
                i <- i+1
                
                companaResultsList[[i]] <- companaResultsDF
                
              } # tes
              
            } # samp
            
          }
        }
      }
    }
  }
  
  return(do.call(rbind, companaResultsList))
  
}