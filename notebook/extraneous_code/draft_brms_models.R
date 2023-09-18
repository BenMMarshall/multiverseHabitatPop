
library(dplyr)
library(brms)

targets::tar_load("areaBasedResults")
targets::tar_load("ssfSampled")
targets::tar_load("poisResults")

names(areaBasedResults)
names(ssfSampled)
names(poisResults)

modelData <- areaBasedResults %>% 
  # track freq translated into the more intuitive track per hour
  dplyr::mutate(trackFreq = round(1/as.numeric(trackFreq), digits = 2)) %>% 
  dplyr::mutate(medEst = median(companaLambda),
                absDeltaEst = companaLambda - medEst) %>% 
  dplyr::mutate(
    sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
    trackFreqScaled = (trackFreq - mean(trackFreq))/sd(trackFreq),
    trackDuraScaled = (trackDura - mean(trackDura))/sd(trackDura),
    contourScaled = (contour - mean(contour))/sd(contour),
    availablePointsScaled  = (availablePoints - mean(availablePoints))/sd(availablePoints)
  )

modelData <- ssfSampled %>% 
  # track freq translated into the more intuitive track per hour
  dplyr::mutate(trackFreq = round(1/as.numeric(trackFreq), digits = 2)) %>% 
  dplyr::mutate(medEst = median(modelAvg),
                absDeltaEst = modelAvg - medEst) %>% 
  dplyr::mutate(
    sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
    trackFreqScaled = (trackFreq - mean(trackFreq))/sd(trackFreq),
    trackDuraScaled = (trackDura - mean(trackDura))/sd(trackDura),
    availablePerStepScaled  = (availablePerStep - mean(availablePerStep))/sd(availablePerStep)
  )

modelData <- poisResults %>% 
  # track freq translated into the more intuitive track per hour
  dplyr::mutate(trackFreq = round(1/as.numeric(trackFreq), digits = 2)) %>% 
  dplyr::mutate(medEst = median(mean),
                absDeltaEst = mean - medEst) %>% 
  dplyr::mutate(
    sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
    trackFreqScaled = (trackFreq-mean(trackFreq))/sd(trackFreq),
    trackDuraScaled = (trackDura-mean(trackDura))/sd(trackDura),
    availablePerStepScaled  = (availablePerStep - mean(availablePerStep))/sd(availablePerStep)
  )

modelData$availablePerStep

formAbsDelta_areaBased <- brms::bf(
  absDeltaEst ~ 1 + sampleSizeScaled + trackFreqScaled + trackDuraScaled +
    type + areaMethod + contourScaled +
    availablePointsScaled + samplingPattern + test +
    (1|sampleID)
)

formAbsDelta_ssf <- brms::bf(
  absDeltaEst ~ 1 + sampleSizeScaled + trackFreqScaled + trackDuraScaled +
    modelFormula + stepDist + turnDist + availablePerStep + averagingMethod +
    (1|sampleID)
)

formAbsDelta_pois <- brms::bf(
  absDeltaEst ~ 1 + sampleSizeScaled + trackFreqScaled + trackDuraScaled +
    modelFormula + stepDist + turnDist + availablePerStepScaled +
    (1|sampleID)
)

brms::get_prior(formAbsDelta_pois,
                data = modelData)

# library(ggplot2)
# ggplot() +
#   geom_density(data = data.frame(x = rcauchy(20000, location = 0.1, scale = 3)),
#                aes(x = x))

brmprior_areaBased <- c(
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

brmprior_ssf <- c(
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

brmprior_pois <- c(
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
