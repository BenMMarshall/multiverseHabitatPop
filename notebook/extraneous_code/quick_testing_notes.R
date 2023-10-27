library(dplyr)
lapply(c("tibble",
         "qs",
         "dplyr",
         "stringr",
         "abmAnimalMovement",
         "INLA",
         "TwoStepCLogit",
         "MuMIn",
         "adehabitatHS",
         "amt",
         "here",
         "ggplot2",
         "patchwork",
         "bayesplot",
         "tidybayes"
), require, character.only = TRUE)

targets::tar_load("modelExtracts")

list2env(modelExtracts, envir = environment())
betasOutputs
unique(betasOutputs$.variable)

palette <- multiverseHabitat::get_palette()

modelPalette <- unname(palette[1:4])

library(ggtext)

names(modelPalette) <- c(
  "<b style='color:#AD6DED'>Two-step</b>",
  "<b style='color:#7D26D4'>Poisson</b>",
  "<b style='color:#4F0E99'>Step Selection</b>",
  "<b style='color:#E87D13'>Area Based</b>")

betasOutputsPlotData <- betasOutputs %>% 
  mutate(facetSplit = factor(case_when(
    .variable %in% c("b_Intercept",
                     "b_sampleSizeScaled",
                     "b_trackDuraScaled",
                     "b_trackFreqScaled") ~  "Sampling Choices",
    .variable %in% c("b_availablePerStepScaled",
                     "b_stepDistgamma",
                     "b_turnDistvonmises",
                     "b_modelFormulamf.ss") ~  "Step Generation Choices",
    .variable %in% c("b_averagingMethodNaiveaverage") ~  "<b style='color:#4F0E99'>SSF Only</b>",
    .variable %in% c("b_areaMethodMCP",
                     "b_contourScaled",
                     "b_typeIII",
                     "b_availablePointsScaled",
                     "b_samplingPatternst",
                     "b_testrandomisation") ~  "<b style='color:#E87D13'>Area Based Only</b>",
  ), levels = c(
    "Sampling Choices",
    "Step Generation Choices",
    "<b style='color:#4F0E99'>SSF Only</b>",
    "<b style='color:#E87D13'>Area Based Only</b>")
  )) %>% 
  mutate(
    .variable = factor(case_when(
      .variable == "b_trackDuraScaled" ~ "\u03B2 Tracking Duration",
      .variable == "b_trackFreqScaled" ~ "\u03B2 Tracking Frequency",
      .variable == "b_sampleSizeScaled" ~ "\u03B2 Sample Size",
      .variable == "b_modelFormulamf.ss" ~ "\u03B2 Model Formula: Not Integrated",
      .variable == "b_stepDistgamma" ~ "\u03B2 Step Distribution: Gamma",
      .variable == "b_turnDistvonmises" ~ "\u03B2 Turn Distribution: Von Mises",
      .variable == "b_availablePerStepScaled" ~ "\u03B2 Available Points Per Step",
      .variable == "b_averagingMethodNaiveaverage" ~ "\u03B2 Averaging Method: Naive",
      .variable == "b_areaMethodMCP" ~ "\u03B2 Available Area: MCP",
      .variable == "b_samplingPatternst" ~ "\u03B2 Sampling Pattern: Stratified",
      .variable == "b_typeIII" ~ "\u03B2 Desigen Type: III",
      .variable == "b_contourScaled" ~ "\u03B2 Available Area Contour",
      .variable == "b_availablePointsScaled" ~ "\u03B2 Available Points Multipiler",
      .variable == "b_testrandomisation" ~ "\u03B2 Compana Test: Randomisation"
    ),
    levels = rev(c(
      "\u03B2 Sample Size",
      "\u03B2 Tracking Duration",
      "\u03B2 Tracking Frequency",
      
      "\u03B2 Available Points Per Step",
      "\u03B2 Step Distribution: Gamma",
      "\u03B2 Turn Distribution: Von Mises",
      "\u03B2 Model Formula: Not Integrated",
      
      "\u03B2 Averaging Method: Naive",
      
      "\u03B2 Available Area: MCP",
      "\u03B2 Available Area Contour",
      "\u03B2 Desigen Type: III",
      "\u03B2 Available Points Multipiler",
      "\u03B2 Sampling Pattern: Stratified",
      "\u03B2 Compana Test: Randomisation"
    ))
    )) %>% 
  mutate(model = factor(case_when(
    model == "areaBrms" ~ "<b style='color:#E87D13'>Area Based</b>",
    model == "ssfBrms" ~ "<b style='color:#4F0E99'>Step Selection</b>",
    model == "poisBrms" ~ "<b style='color:#7D26D4'>Poisson</b>",
    model == "twoStepBrms" ~ "<b style='color:#AD6DED'>Two-step</b>"),
    levels = c(
      "<b style='color:#E87D13'>Area Based</b>",
      "<b style='color:#4F0E99'>Step Selection</b>",
      "<b style='color:#7D26D4'>Poisson</b>",
      "<b style='color:#AD6DED'>Two-step</b>"
    ))
  ) %>% 
  filter(!.variable == "b_Intercept")


gradLimits <- range(c(betasOutputsPlotData$.lower, betasOutputsPlotData$.upper))
labelLocation <- data.frame(gradLimits)
labelLocationText <- data.frame(gradIndent = gradLimits + c(0.15, -0.15))
labelText <- c("Closer to median\npreference estimate",
               "Farther from median\npreference estimate")
arrowAdj <- c(0.015, -0.015)
facetSplit <- factor(rep("<b style='color:#E87D13'>Area Based Only</b>", 2), levels = c(
  "Sampling Choices",
  "Step Generation Choices",
  "<b style='color:#4F0E99'>SSF Only</b>",
  "<b style='color:#E87D13'>Area Based Only</b>")
)
hjust <- c(1,0)

annotationDF <- cbind(labelLocation, labelLocationText, labelText, arrowAdj,
                      facetSplit, hjust)


print("\u25B2")

modelLabels <- tribble(
  ~x, ~y, ~text, ~hjust, ~vjust, ~facetSplit, ~labCol,
  3.05,   0.2, "<b style='color:#AD6DED'>\u25CF = Two-step</b>", 0.5, 0, "Step Generation Choices", "#AD6DED",
  2.45, -0.42, "<b style='color:#7D26D4'>\u25B2 = Poisson</b>", 0.5, 1, "Step Generation Choices", "#7D26D4",
  1,   -0.5, "<b style='color:#4F0E99'>\u25C6 = Step Selection</b>", 0.5, 0, "<b style='color:#4F0E99'>SSF Only</b>", "#4F0E99",
  4,   0.25, "<b style='color:#E87D13'>\u25BC = Area Based</b>", 0.5, 1, "<b style='color:#E87D13'>Area Based Only</b>", "#E87D13"
)
modelLabels <- modelLabels %>% 
  mutate(facetSplit = factor(facetSplit, levels = c(
    "Sampling Choices",
    "Step Generation Choices",
    "<b style='color:#4F0E99'>SSF Only</b>",
    "<b style='color:#E87D13'>Area Based Only</b>")))

arrowsDF <- tribble(
  ~x, ~y, ~xend, ~yend, ~facetSplit, ~labCol,
  3.05,   0.2, 1.4, 0.24, "Step Generation Choices", "#AD6DED",
  2.45,   -0.42, 4, -0.58, "Step Generation Choices", "#7D26D4",
  1,   -0.5, 1, -0.32, "<b style='color:#4F0E99'>SSF Only</b>", "#4F0E99",
  4,   0.25, 5, 0.02, "<b style='color:#E87D13'>Area Based Only</b>", "#E87D13"
)
arrowsDF <- arrowsDF %>% 
  mutate(facetSplit = factor(facetSplit, levels = c(
    "Sampling Choices",
    "Step Generation Choices",
    "<b style='color:#4F0E99'>SSF Only</b>",
    "<b style='color:#E87D13'>Area Based Only</b>")))

allEffectsPlot <- betasOutputsPlotData %>% 
  ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.5, alpha = 0.9, colour = "#403F41",
             linetype = 1) +
  geom_errorbar(aes(x = .variable, ymin = .lower, ymax = .upper, 
                    colour = model),
                position = position_dodge(0.75), width = 0) +
  geom_vline(xintercept = seq(0.5,20.5,1), linewidth = 0.25, alpha = 0.5, colour = "#403F41",
             linetype = 2) +
  geom_richtext(data = modelLabels,
            aes(x = x, y = y, label = text,
                label.color = labCol,
                hjust = hjust, vjust = vjust),
            label.size = 0.75, label.r = unit(0.1, "lines"),
            fill = "#FFFFFF", alpha = 1) +
  geom_curve(data = arrowsDF,
            aes(x = x, xend = xend, y = y, yend = yend,
                colour = labCol),
            curvature = 0.35,
            arrow = arrow(angle = 30, type = "closed", length = unit(2, "mm")))+
  geom_segment(data = annotationDF,
               aes(x = -0.2, xend = -0.2,
                   y = 0.02, yend = gradIndent + arrowAdj),
               colour = "#9F9FA0", 
               arrow = arrow(angle = 30, type = "closed", length = unit(2, "mm")),
               linewidth = 1.25) +
  geom_text(data = annotationDF,
            aes(x = -0.1, y = gradIndent,
                label = labelText, hjust = hjust),
            colour = "#9F9FA0", vjust = 0.5, lineheight = 0.95,
            size = 3, fontface = 4) +
  geom_point(aes(x = .variable, y = .value, 
                 colour = model, shape = model, fill = model),
             position = position_dodge(0.75)) +
  geom_point(data = annotationDF,
             aes(x = -0.75, y = 0)) +
  scale_fill_manual(values = modelPalette,
                    breaks = c(
                      "<b style='color:#AD6DED'>Two-step</b>",
                      "<b style='color:#7D26D4'>Poisson</b>",
                      "<b style='color:#4F0E99'>Step Selection</b>",
                      "<b style='color:#E87D13'>Area Based</b>")) +
  scale_colour_manual(values = modelPalette,
                      breaks = c(
                        "<b style='color:#AD6DED'>Two-step</b>",
                        "<b style='color:#7D26D4'>Poisson</b>",
                        "<b style='color:#4F0E99'>Step Selection</b>",
                        "<b style='color:#E87D13'>Area Based</b>")) +
  scale_shape_manual(values = c(21,24,23,25),
                     breaks = c(
                       "<b style='color:#AD6DED'>Two-step</b>",
                       "<b style='color:#7D26D4'>Poisson</b>",
                       "<b style='color:#4F0E99'>Step Selection</b>",
                       "<b style='color:#E87D13'>Area Based</b>")) +
  facet_grid(vars(facetSplit), scales = "free", space = "free", switch = "y") +
  coord_flip() +
  labs(x = "", y = "Beta point estimate and 95% HDCI",
       colour = "Model", shape = "Model", fill = "Model") +
  theme_bw() +
  theme(
    line = element_line(colour = palette["coreGrey"]),
    text = element_text(colour = palette["coreGrey"]),
    strip.background = element_blank(),
    strip.text.y.left = element_markdown(face = 4, hjust = 0, vjust = 1, angle = 0,
                                         margin = margin(-5,10,0,-152)),
    # strip.text.y.left = element_text(angle = 0, margin = margin(-8,10,0,0)),
    # axis.text.y.left = element_text(margin = margin(0,-119,0,80)), # 2nd value needed to alligns with facet, 4th gives space left
    axis.title.x = element_text(margin = margin(5,0,0,0)),
    axis.ticks.y.left = element_blank(),
    axis.line.x = element_line(),
    strip.clip = "off",
    panel.border = element_blank(),
    panel.spacing = unit(18, "pt"),
    legend.position = "none",
    # legend.position = c(0.1, 0.2),
    # legend.title = element_text(face = 4),
    # legend.text = element_markdown(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

allEffectsPlot

ggsave(allEffectsPlot,
       filename = here("notebook", "figures", "_allEffectsPlot.png"),
       dpi = 300, width = 260, height = 200,
       units = "mm")



###


resultsData %>% 
  mutate(key = paste0(sampleID, trackFreq, trackDura, modelFormula, availablePerStep, stepDist,
                      turnDist)) %>% 
  filter(key == "samp1115mf.is15expvonmises")

prefDiffDF <- resultsData %>% 
  mutate(key = paste0(sampleID, trackFreq, trackDura, modelFormula, availablePerStep, stepDist,
                      turnDist)) %>% 
  group_by(key) %>% 
  summarise(prefDiff = diff(mean))

resultsData %>% 
  filter(term == "layerc2") %>% 
  mutate(key = paste0(sampleID, trackFreq, trackDura, modelFormula, availablePerStep, stepDist,
                      turnDist)) %>% 
  left_join(prefDiffDF)

targets::tar_load("poisResults")
targets::tar_load("poisBrms")
targets::tar_load("ssfBrms")
targets::tar_load("areaBrms")

summary(areaBrms)


tar_load("sampDuraFreqData_15_1")

optionsList_pois <- list(
  MethodPois_mf = c("mf.is"),
  MethodPois_sd = c("gamma"),
  MethodPois_td = c("vonmises"),
  MethodPois_as = c(2, 10)
)

allIndividualData = sampDuraFreqData_15_1
sampleGroups = optionsList_samples[1]
optionsList = optionsList_pois

start <- Sys.time()
tsRes <- method_twoStep(sampDuraFreqData_15_1, sampleGroups = optionsList_samples[7],
                        optionsList = optionsList_pois)
end <- Sys.time()
difftime(end, start)

inlaRes <- method_pois_inla(sampDuraFreqData_15_1, sampleGroups = optionsList_samples[1],
                            optionsList = optionsList_pois)

inlaRes




##



paste(round(modelExtracts$betasOutputs[modelExtracts$betasOutputs$model == "ssfBrms" &
                                         modelExtracts$betasOutputs$.variable == "b_sampleSizeScaled",
                                       c(".lower", ".upper")],
            digits = 2), collapse = " - ")






targets::tar_load("ssfOUT_15_1")
targets::tar_load("ssfResults")

allIndividualData <- sampDuraFreqData_15_1
sampleGroups <- optionsList_samples
optionsList <- optionsList_pois

length(ssfOUT_15_1$simData_i001)
ssfOUT_15_1$simData_i001[[3]]$model$model
ssfOUT_15_1$simData_i001[[2]]$model$model
ssfOUT_15_1$simData_i001[[1]]$model$model


ssfOUT$model$residuals <- NULL
ssfOUT$model$y <- NULL
ssfOUT$sl_ <- NULL
ssfOUT$ta_ <- NULL
ssfOUT$more <- NULL
ssfOUT$model$model <- NULL
ssfOUT$model$terms <- NULL
ssfOUT$model$xlevels <- NULL
ssfOUT$model$contrasts <- NULL
ssfOUT$model$call <- NULL
ssfOUT$model$linear.predictors <- NULL
ssfOUT$model$formula

MuMIn::model.avg(ssfOUT$model,
                 ssfOUT$model)

object.size(ssfOUT)
length(serialize(ssfOUT, NULL))
saveRDS(ssfOUT, file = "ssfOUT.Rds")
qs::qsave(ssfOUT$model, file = "ssfOUT")

lapply(ssfOUT$model, function(x){
  object.size(x)
  length(serialize(x, NULL))
})

names(ssfOUT_15_1)

allIndividualData <- sampDuraFreqData_60_6

sum(is.na(allIndividuals$simData_i001$locations$datetime))
sum(is.na(as.POSIXct(allIndividuals$simData_i001$locations$datetime, 
                     format = "%Y-%m-%d %H:%M:%S")))
sum(is.na(sampDuraFreqData_60_6$simData_i001$locations$datetime))

movementData <- allIndividualData$simData_i001$locations
movementData[which(movementData$minute == 0 |
                     movementData$minute == 30),]

as.POSIXct(allIndividuals$simData_i001$locations)

library(amt)
movementDataNest$data[[1]]$t
movementDataNest$id

movementData %>% 
  tidyr::nest(moveData = -id)



as.POSIXct(movementData$datetime, format = "%y-%m-%d %H:%M:%S",
           tz = "UTC")


a <- try(
  sum(f)
)

tar_load("poisOUT_60_6")

poisOUT_60_6

# inlaResults <- inlaOUT$summary.fixed[2,]


