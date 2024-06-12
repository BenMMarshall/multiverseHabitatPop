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
area_based_calculations <- function(availUseData, sampleGroups, optionsList, optionsListArea){
  # availUseData <- areaBasedAvailUse_60_12
  optionsType <- unique(availUseData$type)
  optionsMethod <- unique(availUseData$method)
  optionsContour <- unique(availUseData$contour)
  optionsAPoints <- unique(availUseData$availablePoints)
  optionsSPSamp <- unique(availUseData$samplingPattern)
  optionsLand <- unique(availUseData$classLandscape)
  
  # optionsList$areaBasedMethod
  optionsTest <- optionsListArea$areaBasedTest
  
  listLength <- length(optionsType)*
    length(optionsMethod)*
    length(optionsLand)*
    length(optionsContour)*
    length(optionsAPoints)*
    length(optionsSPSamp)*
    length(optionsTest)*
    length(names(sampleGroups))
  
  companaResultsList <- vector("list", length = listLength)
  i <- 0
  for(typ in optionsType){
    # typ <- optionsType[1]
    for(met in optionsMethod){
      # met <- optionsMethod[1]
      for(con in optionsContour){
        # con <- optionsContour[1]
        for(aPo in optionsAPoints){
          # aPo <- optionsAPoints[1]
          for(spS in optionsSPSamp){
            # spS <- optionsSPSamp[1]
            for(land in optionsLand){
              # land <- optionsLand[1]
              for(sampID in names(sampleGroups)){
                # sampID <- names(sampleGroups)[1]
                IDs <- optionsList_samples[[sampID]]
                
                IDs <- paste0(stringr::str_extract(availUseData$id[1], "^.*_i"),
                              sprintf("%03d", IDs))
                
                use <- availUseData %>% 
                  dplyr::filter(type == typ,
                                classLandscape == land,
                                method == met,
                                contour == con,
                                availablePoints == aPo,
                                samplingPattern == spS) %>% 
                  dplyr::filter(id %in% IDs) %>% 
                  dplyr::select(used_c0, used_c2)
                
                avail <- availUseData %>% 
                  dplyr::filter(type == typ,
                                classLandscape == land,
                                method == met,
                                contour == con,
                                availablePoints == aPo,
                                samplingPattern == spS) %>% 
                  dplyr::filter(id %in% IDs) %>% 
                  dplyr::select(avail_c0, avail_c2)
                
                names(use) <- c("c0", "c2")
                names(avail) <- c("c0", "c2")
                
                for(tes in optionsTest){
                  # tes <- optionsTest[1]
                  companaOUT <- compana(used = use, avail = avail,
                                        test = tes)
                  
                  # names(companaOUT$rank[1])
                  # companaOUT$rm
                  # companaOUT$rm
                  # companaOUT$rmv
                  
                  companaResultsDF <- data.frame(
                    sampleID = sampID,
                    species = availUseData$species[1],
                    sampleSize = length(IDs),
                    trackFreq = availUseData$trackFreq[1],
                    trackDura = availUseData$trackDura[1],
                    classLandscape = land,
                    type = typ,
                    areaMethod = met,
                    contour = con,
                    availablePoints = aPo,
                    samplingPattern = spS,
                    test = tes,
                    companaHabDiff = companaOUT$rmv[2,1],
                    companaLambda = companaOUT$test["Lambda"],
                    companaP = companaOUT$test["P"]
                  )
                  
                  i <- i+1
                  
                  companaResultsList[[i]] <- companaResultsDF
                  
                } # tes
                
              } # samp
            } #classLandscape
          }
        }
      }
    }
  }
  
  return(do.call(rbind, companaResultsList))
  
}