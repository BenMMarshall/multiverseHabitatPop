
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

optionsList_samples <- list(c(1,2))

names(optionsList_samples) <- paste0("samp", 1:length(optionsList_samples))

sampID <- names(optionsList_samples)[1]

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



