
# can we see the trade-off between sample size and tracking regime, ie can be
# guide what would be worth inversting in

library(here)
library(dplyr)
library(ggplot2)

resultsData <- read.csv(here("data", "poisEstimateOutputs.csv.gz"))

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
  ) %>% 
  mutate(dataQuant = sampleSize * (trackFreq*24) * trackDura)

modelData %>% 
  ggplot() +
  geom_point(aes(x = dataQuant, y = prefDiff,
             colour = sampleSize), position = position_jitter(),
             alpha = 0.2) +
  facet_grid(cols = vars(classLandscape)) +
  scale_x_log10()

modelData %>% 
  ggplot() +
  geom_point(aes(x = dataQuant, y = prefDiff,
             colour = trackFreq)) +
  facet_grid(cols = vars(classLandscape)) +
  scale_x_log10()

modelData %>% 
  filter(trackFreq == 1) %>% 
  filter((sampleSize == 5 & trackDura == 30) | (sampleSize == 10 & trackDura == 15)) %>% 
  ggplot() +
  geom_point(aes(x = sampleSize, y = prefDiff,
                 colour = dataQuant, shape = as.factor(trackDura)), position = position_jitter(),
             alpha = 1) +
  facet_grid(cols = vars(classLandscape)) +
  scale_x_log10()

# those over 480 could be used for comparisons???
length(table(modelData$dataQuant))

