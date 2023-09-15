#' Sample individual SSF models
#'
#' @name sample_ssf_results
#' @description create a population estimate
#' @param ssfResults all sff models and corresponding options data.frame
#' @param sampleGroups list if IDs for the sample groups
#' @return Population estimates for SSF based methods.
#'
#' @export
sample_ssf_results <- function(ssfResults, sampleGroups){
  
  ssfSampledList <- vector("list", length = length(names(ssfResults)) *
                             length(names(sampleGroups)))
  names(ssfSampledList) <- do.call(paste0, expand.grid(names(ssfResults),
                                                       paste0("_", names(sampleGroups))))
  for(regime in names(ssfResults)){
    # regime <- "ssfOUT_15_1"
    ssfRegimeResults <- ssfResults[[regime]]
    
    for(sampID in names(sampleGroups)){
      # sampID <- names(optionsList_samples)[1]
      
      IDs <- sampleGroups[[sampID]]
      IDs <- paste0("simData_i", sprintf("%03d", IDs))
      
      # bit convoluted but this can get the ID (ie species)
      # IDs <- paste0(stringr::str_extract(ssfRegimeResults$simData_i001[[1]]$options$id[1],
      #                                    "^.*_i"),
      #               sprintf("%03d", IDs))
      
      ssfSampleList <- ssfRegimeResults[names(ssfRegimeResults) %in% IDs]
      
      sampleModelList <- lapply(ssfSampleList, function(x){
        x[[1]]$model
      })
      sampleEstList <- lapply(ssfSampleList, function(x){
        x[[1]]$options
      })
      sampleEst <- do.call(rbind, sampleEstList)
      
      # calculate CI surrounding the naive mean
      nEst <- length(sampleEst$Estimate)
      meanEst <- mean(sampleEst$Estimate)
      sdEst <- sd(sampleEst$Estimate)
      marginEst <- qt(0.975, df = nEst - 1) * sdEst/sqrt(nEst)
      lowEst <- meanEst - marginEst
      highEst <- meanEst + marginEst
      
      modelAvgOUT <- MuMIn::model.avg(sampleModelList)
      modelAvgOUTCI <- confint(modelAvgOUT)
      
      ssfResults$ssfOUT_15_1$simData_i001[[1]]$options$trackFreq
      
      modelAvgData <- data.frame(
        sampleID = sampID,
        sampleSize = length(IDs),
        trackFreq = sampleEst$trackFreq[1],
        trackDura = sampleEst$trackDura[1],
        method = sampleEst$analysis[1],
        modelFormula = sampleEst$modelForm[1],
        stepDist = sampleEst$stepDist[1],
        turnDist = sampleEst$turnDist[1],
        availableSteps = sampleEst$availablePerStep[1],
        averagingMethod = "Model average",
        modelAvg = modelAvgOUT$coefficients[1,1],
        modelAvgSE = summary(modelAvgOUT)$coefmat.full[1,2],
        modelAvgLower = modelAvgOUTCI[1,1],
        modelAvgUpper = modelAvgOUTCI[1,2]
      )
      
      naiveAvgData <- data.frame(
        sampleID = sampID,
        sampleSize = length(IDs),
        trackFreq = sampleEst$trackFreq[1],
        trackDura = sampleEst$trackDura[1],
        method = sampleEst$analysis[1],
        modelFormula = sampleEst$modelForm[1],
        stepDist = sampleEst$stepDist[1],
        turnDist = sampleEst$turnDist[1],
        availableSteps = sampleEst$availablePerStep[1],
        averagingMethod = "Naive average",
        modelAvg = mean(sampleEst$Estimate),
        modelAvgSE = sd(sampleEst$Estimate)/sqrt(length(sampleEst$Estimate)),
        modelAvgLower = lowEst,
        modelAvgUpper = highEst
      )
      
      ssfSampleOUT <- rbind(modelAvgData, naiveAvgData)
      
      ssfSampledList[[paste0(regime, "_", sampID)]] <- ssfSampleOUT
      
    }
    
    # for(indiID in names(ssfRegimeResults)){
    #   indiID <- "simData_i001"
    #   
    #   ssfRegimeResults[[indiID]]
    #   
    # }
    
    
  }
  return(do.call(rbind, ssfSampledList))
  
}