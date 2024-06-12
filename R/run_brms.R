#' Run all brm models
#'
#' @name run_brms
#' @description A
#' @param resultsData output tar_combine resulting in ssfResults or areaResults
#' @return a
#'
#' @export
run_brms <- function(resultsData,
                     iter = 4000,
                     warmup = 750,
                     thin = 4){
  # resultsData <- poisResults
  if("companaHabDiff" %in% names(resultsData)){
    # resultsData <- areaBasedResults
    modelData <- resultsData %>% 
      # track freq translated into the more intuitive track per hour
      dplyr::mutate(trackFreq = round(1/as.numeric(trackFreq), digits = 2)) %>% 
      dplyr::group_by(classLandscape) %>% 
      dplyr::mutate(medEst = median(companaHabDiff, na.rm = TRUE),
                    absDeltaEst = abs(companaHabDiff - medEst)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
        trackFreqScaled = (trackFreq - mean(trackFreq))/sd(trackFreq),
        trackDuraScaled = (trackDura - mean(trackDura))/sd(trackDura),
        contourScaled = (contour - mean(contour))/sd(contour),
        availablePointsScaled  = (availablePoints - mean(availablePoints))/sd(availablePoints)
      )
    
    formAbsDelta <- brms::bf(
      absDeltaEst ~ 1 + sampleSizeScaled + trackFreqScaled + trackDuraScaled +
        type + areaMethod + contourScaled +
        availablePointsScaled + samplingPattern + test +
        (1|sampleID) + (classLandscape|classLandscape) + (species|species)
    )
    
    brmpriors <- c(
      # data quantity decreases deviation from median effect
      brms::set_prior("cauchy(-0.1, 3)", coef = "sampleSizeScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackFreqScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackDuraScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "availablePointsScaled"),
      # other kept as weak positive priors
      brms::set_prior("cauchy(0.1, 3)", coef = "contourScaled"),
      brms::set_prior("cauchy(0.1, 3)", coef = "samplingPatternst"),
      brms::set_prior("cauchy(0.1, 3)", coef = "testrandomisation"),
      brms::set_prior("cauchy(0.1, 3)", coef = "areaMethodMCP"),
      brms::set_prior("cauchy(0.1, 3)", coef = "typeIII")
    )
    
    modelSave <- here::here("notebook", "modelOutput",
                            "absDelta_areaBased.txt")
    modelFile <- here::here("notebook", "modelOutput",
                            "absDelta_areaBased")
    
  } else if(resultsData$analysis[1] == "ssf"){
    
    modelData <- resultsData %>% 
      # track freq translated into the more intuitive track per hour
      dplyr::mutate(trackFreq = round(1/as.numeric(trackFreq), digits = 2)) %>% 
      dplyr::group_by(classLandscape) %>% 
      dplyr::mutate(medEst = median(modelAvg, na.rm = TRUE),
                    absDeltaEst = abs(modelAvg - medEst)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
        trackFreqScaled = (trackFreq - mean(trackFreq))/sd(trackFreq),
        trackDuraScaled = (trackDura - mean(trackDura))/sd(trackDura),
        availablePerStepScaled  = (availablePerStep - mean(availablePerStep))/sd(availablePerStep)
      )
    
    formAbsDelta <- brms::bf(
      absDeltaEst ~ 1 + sampleSizeScaled + trackFreqScaled + trackDuraScaled +
        modelFormula + stepDist + turnDist + availablePerStepScaled + averagingMethod +
        (1|sampleID) + (classLandscape|classLandscape) + (species|species)
    )
    
    
    # print(brms::get_prior(formAbsDelta,
    #                       data = modelData))
    
    brmpriors <- c(
      # data quantity decreases deviation from median effect
      brms::set_prior("cauchy(-0.1, 3)", coef = "sampleSizeScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackFreqScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackDuraScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "availablePerStepScaled"),
      # other kept as weak positive priors
      brms::set_prior("cauchy(0.1, 3)", coef = "modelFormulamf.ss"),
      brms::set_prior("cauchy(0.1, 3)", coef = "stepDistgamma"),
      brms::set_prior("cauchy(0.1, 3)", coef = "turnDistvonmises"),
      brms::set_prior("cauchy(0.1, 3)", coef = "averagingMethodNaiveaverage")
    )
    
    modelSave <- here::here("notebook", "modelOutput",
                            "absDelta_ssf.txt")
    modelFile <- here::here("notebook", "modelOutput",
                            "absDelta_ssf")
    
  } else if(resultsData$analysis[1] == "Poisson"){
    # resultsData <- poisResults
    prefDiffDF <- resultsData %>% 
      # filter(term == "layerc2") %>% 
      mutate(key = paste0(sampleID, trackFreq, trackDura, modelFormula, availablePerStep, stepDist, classLandscape,
                          turnDist, species)) %>% 
      group_by(key) %>% 
      summarise(prefDiff = diff(mean))
    
    resultsData <- resultsData %>% 
      filter(term == "layerc2") %>% 
      mutate(key = paste0(sampleID, trackFreq, trackDura, modelFormula, availablePerStep, stepDist, classLandscape,
                          turnDist, species)) %>% 
      left_join(prefDiffDF)
    
    modelData <- resultsData %>% 
      # track freq translated into the more intuitive track per hour
      dplyr::mutate(trackFreq = round(1/as.numeric(trackFreq), digits = 2)) %>% 
      dplyr::group_by(classLandscape) %>% 
      dplyr::mutate(medEst = median(prefDiff, na.rm = TRUE),
                    absDeltaEst = abs(prefDiff - medEst)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
        trackFreqScaled = (trackFreq-mean(trackFreq))/sd(trackFreq),
        trackDuraScaled = (trackDura-mean(trackDura))/sd(trackDura),
        availablePerStepScaled  = (availablePerStep - mean(availablePerStep))/sd(availablePerStep)
      )
    
    formAbsDelta <- brms::bf(
      absDeltaEst ~ 1 + sampleSizeScaled + trackFreqScaled + trackDuraScaled +
        modelFormula + stepDist + turnDist + availablePerStepScaled +
        (1|sampleID) + (classLandscape|classLandscape) + (species|species)
    )
    
    brmpriors <- c(
      # data quantity decreases deviation from median effect
      brms::set_prior("cauchy(-0.1, 3)", coef = "sampleSizeScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackFreqScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackDuraScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "availablePerStepScaled"),
      # other kept as weak positive priors
      brms::set_prior("cauchy(0.1, 3)", coef = "modelFormulamf.ss"),
      brms::set_prior("cauchy(0.1, 3)", coef = "stepDistgamma"),
      brms::set_prior("cauchy(0.1, 3)", coef = "turnDistvonmises")
    )
    
    modelSave <- here::here("notebook", "modelOutput",
                            "absDelta_pois.txt")
    modelFile <- here::here("notebook", "modelOutput",
                            "absDelta_pois")
    
  } else if(resultsData$analysis[1] == "TwoStep"){
    
    modelData <- resultsData %>% 
      # track freq translated into the more intuitive track per hour
      dplyr::mutate(trackFreq = round(1/as.numeric(trackFreq), digits = 2)) %>% 
      dplyr::group_by(classLandscape) %>% 
      dplyr::mutate(medEst = median(twoStepBeta, na.rm = TRUE),
                    absDeltaEst = abs(twoStepBeta - medEst)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(
        sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
        trackFreqScaled = (trackFreq-mean(trackFreq))/sd(trackFreq),
        trackDuraScaled = (trackDura-mean(trackDura))/sd(trackDura),
        availablePerStepScaled  = (availablePerStep - mean(availablePerStep))/sd(availablePerStep)
      )
    
    formAbsDelta <- brms::bf(
      absDeltaEst ~ 1 + sampleSizeScaled + trackFreqScaled + trackDuraScaled +
        modelFormula + stepDist + turnDist + availablePerStepScaled +
        (1|sampleID) + (classLandscape|classLandscape) + (species|species)
    )
    
    brmpriors <- c(
      # data quantity decreases deviation from median effect
      brms::set_prior("cauchy(-0.1, 3)", coef = "sampleSizeScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackFreqScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "trackDuraScaled"),
      brms::set_prior("cauchy(-0.1, 3)", coef = "availablePerStepScaled"),
      # other kept as weak positive priors
      brms::set_prior("cauchy(0.1, 3)", coef = "modelFormulamf.ss"),
      brms::set_prior("cauchy(0.1, 3)", coef = "stepDistgamma"),
      brms::set_prior("cauchy(0.1, 3)", coef = "turnDistvonmises")
    )
    
    modelSave <- here::here("notebook", "modelOutput",
                            "absDelta_twoStep.txt")
    modelFile <- here::here("notebook", "modelOutput",
                            "absDelta_twoStep")
  }
  
  modOUT_absDelta <- brms::brm(formula = formAbsDelta,
                               data = modelData,
                               family = gaussian,
                               prior = brmpriors,
                               warmup = warmup, iter = iter, chains = 4,
                               cores = 4,
                               thin = thin,
                               # control = list(adapt_delta = 0.90,
                               #                max_treedepth = 15),
                               seed = 1,
                               save_pars = brms::save_pars(all = TRUE),
                               save_model = modelSave,
                               file = modelFile
  )
  
  return(modOUT_absDelta)
  
}