#' Sample individual SSF models
#'
#' @name sample_ssf_results
#' @description create a population estimate
#' @param ssfResults all sff models and corresponding options data.frame
#' @param sampleGroups list if IDs for the sample groups
#' @return Population estimates for SSF based methods.
#'
#' @export
sample_ssf_results <- function(allIndividualData, sampleGroups, optionsList){
  
  optionsForm <- optionsList$MethodSSF_mf
  optionsLand <- optionsList$MethodSSF_land
  optionsASteps <- optionsList$MethodSSF_as
  optionsStepD <- optionsList$MethodSSF_sd
  optionsTurnD <- optionsList$MethodSSF_td
  
  listLength <- length(optionsForm) *
    length(optionsLand)*
    length(optionsASteps) *
    length(optionsStepD) *
    length(optionsTurnD) *
    length(names(sampleGroups)) 
  # length(names(ssfResults))
  
  ssfSampledList <- vector("list", length = listLength)
  i <- 0
  # for(regime in names(ssfResults)){
  #   # regime <- "ssfOUT_15_1"
  #   ssfRegimeResults <- ssfResults[[regime]]
  # ssfRegimeResults <- ssfResults
  
  for(sampID in names(sampleGroups)){
    # sampID <- names(optionsList_samples)[1]
    
    print(sampID)
    IDs <- sampleGroups[[sampID]]
    IDs <- paste0("simData_i", sprintf("%03d", IDs))
    
    # bit convoluted but this can get the ID (ie species)
    # IDs <- paste0(stringr::str_extract(ssfRegimeResults$simData_i001[[1]]$options$id[1],
    #                                    "^.*_i"),
    #               sprintf("%03d", IDs))
    
    #### THIS NEEDDS TO GET SAMPLED INDIVIDUALS FROM sampDuraFreqData_xx_x data
    sampledIndividualData <- allIndividualData[names(allIndividualData) %in% IDs]
    # landscapeList <- allIndividualData$landscape
    # landscape <- landscapeList[names(landscapeList) %in% c("classRaster", "classRasterScram")]
    
    # mf <- "mf.ss"
    # td <- 15
    # tf <- 1
    # step <- "gamma"
    # turn <- "vonmises"
    # as <- 2
    
    for(mf in optionsForm){
      # mf <- optionsForm[1]
      
      # trackDura and trackFreq loops not stricktly required because regime is
      # top level list, but helpful for subsetting and carrying forward info
      for(td in unique(sampledIndividualData[[1]]$trackDura)){
        # td<- unique(sampledIndividualData[[1]])[1]
        for(tf in unique(sampledIndividualData[[1]]$trackFreq)){
          # tf <- unique(sampledIndividualData[[1]]$trackFreq)[1]
          for(step in optionsStepD){
            # step <- optionsStepD[1]
            for(turn in optionsTurnD){
              # turn <- optionsTurnD[1]
              for(as in optionsASteps){
                # as <- optionsASteps[1]
                
                for(land in optionsLand){
                  # land <- unique(names(landscape))[1]
                  optionsList_sffTemp <- list(
                    Method_method = c("ssf"),
                    Method_land = land,
                    MethodSSF_as = as,
                    MethodSSF_mf = mf,
                    MethodSSF_sd = step,
                    MethodSSF_td = turn
                  )
                  
                  sampledIndividualData$landscape <- allIndividualData$landscape
                  # we need to pass the landscape data into the function with the sampledIndidata
                  ssfResultsIndividuals <- wrapper_indi_ssf(allIndividualData = sampledIndividualData,
                                                            optionsList = optionsList_sffTemp)
                  
                  sampleModelList <- lapply(ssfResultsIndividuals, function(x){
                    x[[1]]$model
                  })
                  sampleEstList_sample <- lapply(ssfResultsIndividuals, function(x){
                    x[[1]]$options
                  })
                  sampleEst_sample <- do.call(rbind, sampleEstList_sample)
                  
                  print(all(sampleEst_sample$modelFormula == mf,
                            sampleEst_sample$trackFreq == tf,
                            sampleEst_sample$trackDura == td,
                            sampleEst_sample$classLandscape == land,
                            sampleEst_sample$stepDist == step,
                            sampleEst_sample$turnDist == turn,
                            sampleEst_sample$availablePerStep == as))
                  print(paste(tf,
                              td,
                              land,
                              mf,
                              step,
                              turn,
                              as, sep = " - "))
                  
                  # calculate CI surrounding the naive mean
                  nEst <- length(sampleEst_sample$Estimate)
                  meanEst <- mean(sampleEst_sample$Estimate)
                  sdEst <- sd(sampleEst_sample$Estimate)
                  marginEst <- qt(0.975, df = nEst - 1) * sdEst/sqrt(nEst)
                  lowEst <- meanEst - marginEst
                  highEst <- meanEst + marginEst
                  
                  modelAvgOUT <- MuMIn::model.avg(sampleModelList)
                  # need an try here as confint can fail, in those cases return NA
                  modelAvgOUTCI <- try(
                    confint(modelAvgOUT)
                  )
                  if(class(modelAvgOUTCI)[1] == "try-error"){
                    modelAvgOUTCI <- matrix(NA, 2, 2)
                  }
                  
                  modelAvgData <- data.frame(
                    sampleID = sampID,
                    sampleSize = length(IDs),
                    trackFreq = sampleEst_sample$trackFreq[1],
                    trackDura = sampleEst_sample$trackDura[1],
                    analysis = sampleEst_sample$analysis[1],
                    classLandscape = land,
                    modelFormula = sampleEst_sample$modelFormula[1],
                    stepDist = sampleEst_sample$stepDist[1],
                    turnDist = sampleEst_sample$turnDist[1],
                    availablePerStep = sampleEst_sample$availablePerStep[1],
                    averagingMethod = "Model average",
                    modelAvg = modelAvgOUT$coefficients[1,1],
                    modelAvgSE = summary(modelAvgOUT)$coefmat.full[1,2],
                    modelAvgLower = modelAvgOUTCI[1,1],
                    modelAvgUpper = modelAvgOUTCI[1,2]
                  )
                  
                  naiveAvgData <- data.frame(
                    sampleID = sampID,
                    sampleSize = length(IDs),
                    trackFreq = sampleEst_sample$trackFreq[1],
                    trackDura = sampleEst_sample$trackDura[1],
                    analysis = sampleEst_sample$analysis[1],
                    classLandscape = land,
                    modelFormula = sampleEst_sample$modelFormula[1],
                    stepDist = sampleEst_sample$stepDist[1],
                    turnDist = sampleEst_sample$turnDist[1],
                    availablePerStep = sampleEst_sample$availablePerStep[1],
                    averagingMethod = "Naive average",
                    modelAvg = mean(sampleEst_sample$Estimate),
                    modelAvgSE = sd(sampleEst_sample$Estimate)/sqrt(length(sampleEst_sample$Estimate)),
                    modelAvgLower = lowEst,
                    modelAvgUpper = highEst
                  )
                  
                  ssfSampleOUT <- rbind(modelAvgData, naiveAvgData)
                  
                  i <- i+1
                  # ssfSampledList[[paste0(regime, "_", sampID)]] <- ssfSampleOUT
                  ssfSampledList[[i]] <- ssfSampleOUT
                  
                } # landscape
              }
            }
          }
        }
      }
    }
    
  } # sample loop
  
  # for(indiID in names(ssfRegimeResults)){
  #   indiID <- "simData_i001"
  #   
  #   ssfRegimeResults[[indiID]]
  #   
  # }
  
  # } # regime loop
  return(do.call(rbind, ssfSampledList))
  
}