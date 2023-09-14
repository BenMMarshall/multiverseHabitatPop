
targets::tar_load("sampDuraFreqData_15_1")

mod <- method_indi_ssf(
  allIndividualData = sampDuraFreqData_15_1,
  optionsList = optionsList_sff
)

movementData = sampDuraFreqData_15_1$simData_i002$locations
landscape = sampDuraFreqData_15_1$landscape
# below can all be programmed as single values as the
# targets workflow will be used to feed multiple values
# in
methodForm = "mf.is"
stepDist = "gamma"
turnDist = "vonmises"
availableSteps = 10

mean(c(ssfOUT$model$coefficients[1], ssfOUT2$model$coefficients[1]))
MuMIn::model.avg(ssfOUT$model, ssfOUT2$model)

library(dplyr)

targets::tar_load("ssfResults")

ssfResults$ssfOUT_15_1$simData_i001[[1]]$options
ssfResults$ssfOUT_15_1$simData_i001[[1]]$model


optionsList_samples <- list(c(1,2))

names(optionsList_samples) <- paste0("samp", 1:length(optionsList_samples))

sampID <- names(optionsList_samples)[1]


ssfSampledList <- vector("list", length = length(names(ssfResults)) *
                           length(names(optionsList_samples)))
names(ssfSampledList) <- do.call(paste0, expand.grid(names(ssfResults),
                                                     paste0("_", names(optionsList_samples))))
for(regime in names(ssfResults)){
  # regime <- "ssfOUT_15_1"
  ssfRegimeResults <- ssfResults[[regime]]
  
  for(sampID in names(optionsList_samples)){
    # sampID <- names(optionsList_samples)[1]
    
    IDs <- optionsList_samples[[sampID]]
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
    
    ssfSampleOUT <- data.frame(
      sampleID = sampID,
      sampleSize = length(IDs),
      method = sampleEst$analysis[1],
      modelFormula = sampleEst$modelForm[1],
      stepDist = sampleEst$stepDist[1],
      turnDist = sampleEst$turnDist[1],
      availableSteps = sampleEst$availablePerStep[1],
      naiveAvg = mean(sampleEst$Estimate),
      naiveAvgSE = sd(sampleEst$Estimate)/sqrt(length(sampleEst$Estimate)),
      naiveAvg025Lower = lowEst,
      naiveAvg975Upper = highEst,
      modelAvg = modelAvgOUT$coefficients[1,1],
      modelAvgSE = summary(modelAvgOUT)$coefmat.full[1,2],
      modelAvgLower = modelAvgOUTCI[1,1],
      modelAvgUpper = modelAvgOUTCI[1,2]
    )
    
  }
  
  # for(indiID in names(ssfRegimeResults)){
  #   indiID <- "simData_i001"
  #   
  #   ssfRegimeResults[[indiID]]
  #   
  # }
  
  ssfSampledList[[paste0(regime, "_", sampID)]] <- ssfSampleOUT
  
}
return(ssfSampledList)


optionsTrackFreq <- unique(ssfResults$trackFreq)
optionsTrackDura <- unique(ssfResults$trackDura)
optionsASteps <- unique(ssfResults$availablePerStep)
optionsStepD <- unique(ssfResults$stepDist)
optionsTurnD <- unique(ssfResults$turnDist)
optionsForm <- unique(ssfResults$modelForm)

listLength <- length(optionsTrackFreq) *
  length(optionsTrackDura) *
  length(optionsASteps) *
  length(optionsStepD) *
  length(optionsTurnD) *
  length(optionsForm) *
  length(optionsList_samples)

ssfPopulationList <- vector("list", length = listLength)
i <- 0
for(tf in optionsTrackFreq){
  
  for(td in optionsTrackDura){
    
    for(as in optionsASteps){
      
      for(sdist in optionsStepD){
        
        for(tdist in optionsTurnD){
          
          for(form in optionsForm){
            
            ssfOption <- ssfResults %>% 
              filter(
                # trackFreq == 1,
                # trackDura == 15,
                # availablePerStep == 10,
                # turnDist == "vonmises",
                # stepDist == "gamma",
                # modelForm == "mf.is"
              )
            
            for(sampID in names(optionsList_samples)){
              
              IDs <- optionsList_samples[[sampID]]
              
              IDs <- paste0(stringr::str_extract(ssfOption$id[1], "^.*_i"),
                            sprintf("%03d", IDs))
              
              ssfSample <- ssfOption %>% 
                filter(id %in% IDs)
              
              mean(ssfSample$Estimate)
              
              ssfPopulationOUT
              
              i <- i+1
              ssfPopulationList[[i]] <- ssfPopulationOUT
              
            }
          }
        }
      }
    }
  }
}



