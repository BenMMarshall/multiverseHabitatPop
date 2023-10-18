#' Sample individual SSF models
#'
#' @name sample_ssf_results
#' @description create a population estimate
#' @param ssfResults all sff models and corresponding options data.frame
#' @param sampleGroups list if IDs for the sample groups
#' @return Population estimates for SSF based methods.
#'
#' @export
sample_ssf_results <- function(ssfResults, sampleGroups, optionsList){
  # sampleGroups <- optionsList_samples
  # ssfSampledList <- vector("list", length = length(names(ssfResults)) *
  #                            length(names(sampleGroups)))
  # names(ssfSampledList) <- do.call(paste0, expand.grid(names(ssfResults),
  #                                                      paste0("_", names(sampleGroups))))
  
  optionsForm <- optionsList$MethodSSF_mf
  optionsASteps <- optionsList$MethodSSF_as
  optionsStepD <- optionsList$MethodSSF_sd
  optionsTurnD <- optionsList$MethodSSF_td
  
  listLength <- length(optionsForm) *
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
    ssfRegimeResults <- ssfResults
    
    for(sampID in names(sampleGroups)){
      # sampID <- names(optionsList_samples)[1]
      
      IDs <- sampleGroups[[sampID]]
      IDs <- paste0("simData_i", sprintf("%03d", IDs))
      
      # bit convoluted but this can get the ID (ie species)
      # IDs <- paste0(stringr::str_extract(ssfRegimeResults$simData_i001[[1]]$options$id[1],
      #                                    "^.*_i"),
      #               sprintf("%03d", IDs))
      
      ssfSampleListStart <- ssfRegimeResults[names(ssfRegimeResults) %in% IDs]
      
      # length(ssfSampleList$simData_i001)
      
      # sampleModelList <- lapply(ssfSampleListStart, function(x){
      #   x[[1]]$model
      # })
      # sampleEstList <- lapply(ssfSampleListStart, function(x){
      #   x[[1]]$options
      # })
      # sampleEst <- do.call(rbind, sampleEstList)
      
      sampleEstList_all <- lapply(ssfSampleListStart, function(x){
        zList <- lapply(x, function(y){
          y$options
        })
        z <- do.call(rbind, zList)
      })
      sampleEst_all <- do.call(rbind, sampleEstList_all)
      
      # mf <- "mf.ss"
      # td <- 15
      # tf <- 1
      # step <- "gamma"
      # turn <- "vonmises"
      # as <- 2
      
      
      for(mf in unique(sampleEst_all$modelFormula)){
        # mf <- unique(sampleEst$modelFormula)[1]
        
        # trackDura and trackFreq loops not stricktly required because regime is
        # top level list, but helpful for subsetting and carrying forward info
        for(td in unique(sampleEst_all$trackDura)){
          # td<- unique(sampleEst$trackDura)[1]
          for(tf in unique(sampleEst_all$trackFreq)){
            # tf <- unique(sampleEst$trackFreq)[1]
            for(step in unique(sampleEst_all$stepDist)){
              # step <- unique(sampleEst$stepDist)[1]
              for(turn in unique(sampleEst_all$turnDist)){
                # turn <- unique(sampleEst$turnDist)[1]
                for(as in unique(sampleEst_all$availablePerStep)){
                  # as <- unique(sampleEst$availablePerStep)[1]
                  
                  currentOptions <- sampleEst_all %>% 
                    dplyr::filter(
                      modelFormula == mf,
                      trackDura == td,
                      trackFreq == tf,
                      stepDist == step,
                      turnDist == turn,
                      availablePerStep == as
                    )
                  
                  optionsIndex <- as.numeric(str_extract(row.names(currentOptions),
                                                         "(?<=\\.).{1,}$"))[1]
                  
                  print(
                    paste0(optionsIndex, " --- ",
                           mf, " - ", td, " - ", tf, " - ", step, " - ", turn, " - ", as)
                  )
                  
                  sampleModelList <- lapply(ssfSampleListStart, function(x){
                    x[[optionsIndex]]$model
                  })
                  sampleEstList_sample <- lapply(ssfSampleListStart, function(x){
                    x[[optionsIndex]]$options
                  })
                  sampleEst_sample <- do.call(rbind, sampleEstList_sample)
                  
                  print(all(sampleEst_sample$modelFormula == mf,
                            sampleEst_sample$trackFreq == tf,
                            sampleEst_sample$trackDura == td,
                            sampleEst_sample$stepDist == step,
                            sampleEst_sample$turnDist == turn,
                            sampleEst_sample$availablePerStep == as))
                  
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
                  
                  # ssfResults$ssfOUT_15_1$simData_i001[[1]]$options$trackFreq
                  
                  modelAvgData <- data.frame(
                    sampleID = sampID,
                    sampleSize = length(IDs),
                    trackFreq = sampleEst_sample$trackFreq[1],
                    trackDura = sampleEst_sample$trackDura[1],
                    analysis = sampleEst_sample$analysis[1],
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