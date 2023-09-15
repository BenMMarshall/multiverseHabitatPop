

palette <- multiverseHabitat::get_palette()

library(dplyr)
library(here)
library(ggplot2)
# library(ggridges)
library(tidyr)
library(patchwork)

targets::tar_load("ssfSampled")

palette <- multiverseHabitat::get_palette()

# SSF ---------------------------------------------------------------------

outputResults <- ssfSampled %>% 
  mutate("estimate" = modelAvg)

# convert tf to points/hour to help interpretation
outputResults$trackFreq <- 1/as.numeric(outputResults$trackFreq)
outputResults$trackFreq <- round(outputResults$trackFreq, digits = 2)

levelOrdering <- unique(c(
  unique(outputResults$averagingMethod),
  unique(outputResults$sampleID),
  sort(unique(outputResults$sampleSize)),
  sort(unique(outputResults$trackFreq)),
  sort(unique(outputResults$trackDura)),
  sort(unique(outputResults$modelFormula)),
  sort(unique(outputResults$availableSteps)),
  sort(unique(outputResults$stepDist)),
  sort(unique(outputResults$turnDist))))

outputPlotData <- outputResults %>% 
  select(-method, -sampleID) %>% 
  mutate(across(1:8, as.character)) %>% 
  tidyr::pivot_longer(cols = 1:8, names_to = "variable") %>% 
  dplyr::mutate(
    variable = case_when(
      variable == "trackDura" ~ "Tracking Duration (days)",
      variable == "trackFreq" ~ "Tracking Frequency (points/hour)",
      variable == "sampleSize" ~ "Sample Size (n)",
      variable == "modelFormula" ~ "Model Formula (SSF or iSSF)",
      variable == "availableSteps" ~ "Available Points per Step",
      variable == "stepDist" ~ "Distribution of Step Lengths",
      variable == "turnDist" ~ "Distribution of Turn Angles",
      variable == "averagingMethod" ~ "Model Averaging Method"
    ),
    # sampleID = as.factor(sampleID),
    value = factor(value, levels = levelOrdering)) %>%
  dplyr::group_by(variable, value) %>%
  dplyr::mutate(d_medEst = modelAvg - median(outputResults$modelAvg, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  mutate("estimate" = modelAvg)

# Area based --------------------------------------------------------------

targets::tar_load("areaBasedResults")

outputResults <- areaBasedResults %>% 
  mutate("estimate" = companaLambda)

outputResults$trackFreq <- 1/as.numeric(outputResults$trackFreq)
outputResults$trackFreq <- round(outputResults$trackFreq, digits = 2)

levelOrdering <- unique(c(
  unique(outputResults$sampleID),
  sort(unique(outputResults$sampleSize)),
  sort(unique(outputResults$trackFreq)),
  sort(unique(outputResults$trackDura)),
  unique(outputResults$type),
  unique(outputResults$areaMethod),
  sort(unique(outputResults$contour)),
  sort(unique(outputResults$aPoints)),
  sort(unique(outputResults$spSamp)),
  sort(unique(outputResults$test))))

outputPlotData <- outputResults %>% 
  select(-sampleID) %>% 
  mutate(across(1:9, as.character)) %>% 
  tidyr::pivot_longer(cols = 1:9, names_to = "variable") %>% 
  dplyr::mutate(
    variable = case_when(
      variable == "trackDura" ~ "Tracking Duration (days)",
      variable == "trackFreq" ~ "Tracking Frequency (points/hour)",
      variable == "sampleSize" ~ "Sample Size (n)",
      variable == "type" ~ "Type II or Type III Habitat Selection",
      variable == "aPoints" ~ "Available Points Multiplier",
      variable == "spSamp" ~ "Sampling Pattern",
      variable == "test" ~ "Compana Test Method",
      variable == "contour" ~ "Available Area Contour (%)",
      variable == "areaMethod" ~ "Available Area Method"
    ),
    # sampleID = as.factor(sampleID),
    value = factor(value, levels = levelOrdering)) %>%
  dplyr::group_by(variable, value) %>%
  dplyr::mutate(d_medEst = companaLambda - median(outputResults$companaLambda, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  mutate("estimate" = companaLambda)

# Poisson models ----------------------------------------------------------

targets::tar_load("poisResults")

outputResults <- poisResults %>% 
  mutate("estimate" = mean)

outputResults$trackFreq <- 1/as.numeric(outputResults$trackFreq)
outputResults$trackFreq <- round(outputResults$trackFreq, digits = 2)

levelOrdering <- unique(c(
  unique(outputResults$sampleID),
  sort(unique(outputResults$sampleSize)),
  sort(unique(outputResults$trackFreq)),
  sort(unique(outputResults$trackDura)),
  sort(unique(outputResults$modelForm)),
  sort(unique(outputResults$availablePerStep)),
  sort(unique(outputResults$stepDist)),
  sort(unique(outputResults$turnDist))))

outputPlotData <- outputResults %>% 
  select(-analysis, -sampleID) %>% 
  mutate(across(1:7, as.character)) %>% 
  tidyr::pivot_longer(cols = 1:7, names_to = "variable") %>% 
  dplyr::mutate(
    variable = case_when(
      variable == "trackDura" ~ "Tracking Duration (days)",
      variable == "trackFreq" ~ "Tracking Frequency (points/hour)",
      variable == "sampleSize" ~ "Sample Size (n)",
      variable == "modelForm" ~ "Model Formula (SSF or iSSF)",
      variable == "availablePerStep" ~ "Available Points per Step",
      variable == "stepDist" ~ "Distribution of Step Lengths",
      variable == "turnDist" ~ "Distribution of Turn Angles"
    ),
    # sampleID = as.factor(sampleID),
    value = factor(value, levels = levelOrdering)) %>%
  dplyr::group_by(variable, value) %>%
  dplyr::mutate(d_medEst = mean - median(outputResults$mean, na.rm = TRUE)) %>%
  dplyr::ungroup() %>% 
  mutate("estimate" = mean)


# PLOTS -------------------------------------------------------------------

medData <- outputPlotData %>%
  dplyr::group_by(variable, value) %>%
  dplyr::summarise(modelMedEst = median(estimate, na.rm = TRUE))

(splitSpecCurve <- outputPlotData %>%
    ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.5, alpha = 0.9, colour = "#403F41",
               linetype = 1) +
    geom_point(aes(x = estimate, y = value, colour = d_medEst),
               position = position_jitter(width = 0, height = 0.2), alpha = 0.25,
               pch = 3, size = 0.75) +
    geom_point(data = medData, aes(x = modelMedEst, y = value),
               alpha = 1, size = 1.5, colour = "#FFFFFF") +
    geom_point(data = medData, aes(x = modelMedEst, y = value),
               alpha = 1, size = 1, colour = "#403F41") +
    geom_hline(yintercept = seq(0.5,10.5,1), linewidth = 0.5, alpha = 0.25, colour = "#403F41",
               linetype = 2) +
    facet_grid(variable~., scales = "free_y", space = "free", switch = "y") +
    labs(y = "", x = "Estimate") +
    # scale_colour_gradient2(low = palette["BADGER"],
    #                        mid = palette["coreGrey"],
    #                        high = palette["2"]) +
    theme_bw() +
    theme(
      line = element_line(colour = palette["coreGrey"]),
      text = element_text(colour = palette["coreGrey"]),
      strip.background = element_blank(),
      strip.text = element_text(face = 4, hjust = 1, vjust = 1),
      strip.text.y.left = element_text(angle = 0, margin = margin(-8.5,12,0,0)),
      axis.text.y.left = element_text(margin = margin(0,-165,0,80)), # 2nd value needed to alligns with facet, 4th gives space left
      axis.ticks.y.left = element_blank(),
      axis.line.x = element_line(),
      strip.clip = "off",
      legend.position = "none",
      panel.border = element_blank(),
      panel.spacing = unit(18, "pt"),
      panel.grid = element_blank())
  )

overallMed <- data.frame("medEst" = median(outputResults$estimate, na.rm = TRUE),
                         "indexLoc" = round(nrow(outputResults)/2, digits = 0))

(overallSpecCurve <- outputResults %>%
    dplyr::arrange(estimate) %>%
    dplyr::mutate(index = row_number(),
                  d_medEst = estimate - median(outputResults$estimate, na.rm = TRUE)) %>%
    ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.25, alpha = 0.9, colour = "#403F41",
               linetype = 1) +
    # coord_cartesian(xlim = c(-35, 20)) +
    geom_point(aes(x = estimate, y = index), alpha = 0.25,
               pch = 3, size = 0.75)+
    geom_point(data = data.frame("medEst" = median(outputResults$estimate, na.rm = TRUE),
                                 "indexLoc" = round(nrow(outputResults)/2, digits = 0)),
               aes(x = medEst, y = indexLoc),
               alpha = 1, size = 2.5, colour = "#FFFFFF") +
    geom_point(data = data.frame("medEst" = median(outputResults$estimate, na.rm = TRUE),
                                 "indexLoc" = round(nrow(outputResults)/2, digits = 0)),
               aes(x = medEst, y = indexLoc),
               alpha = 1, size = 2, colour = "#403F41") +
    annotate("text", x = overallMed$medEst + (overallMed$medEst + max(outputResults$estimate))/2,
             y = overallMed$indexLoc, label = "Median",
             fontface = 4, size = 5, colour = palette["coreGrey"],
             hjust = 1, vjust = -0.2) +
    annotate("segment", x = overallMed$medEst + (overallMed$medEst + max(outputResults$estimate))/2,
             xend = overallMed$medEst,
             y = overallMed$indexLoc, yend = overallMed$indexLoc,
             linewidth = 0.75, colour = palette["coreGrey"]) +
    # scale_colour_gradient2(low = palette["BADGER"], mid = palette["coreGrey"], high = palette["2"]) +
    labs(y = "", x = "Estimate") +
    theme_bw() +
    theme(
      line = element_line(colour = palette["coreGrey"]),
      text = element_text(colour = palette["coreGrey"]),
      strip.background = element_blank(),
      strip.text = element_text(face = 4, hjust = 1, vjust = 1),
      strip.text.y.left = element_text(angle = 0, margin = margin(-8,10,0,0)),
      axis.text.y.left = element_blank(),
      axis.ticks.y.left = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_line(),
      strip.clip = "off",
      legend.position = "none",
      panel.border = element_blank(),
      panel.spacing = unit(18, "pt"),
      panel.grid = element_blank())
  )

(specComplete <- overallSpecCurve / splitSpecCurve +
    plot_layout(heights = c(1, 3), guides = "collect"))

ggsave(filename = here("notebook", "figures", "specCurve.png"),
       plot = specComplete,
       width = 360, height = 240, units = "mm", dpi = 300)
# ggsave(filename = here("notebook", "figures", "rsfSpecCurve.pdf"),
#        plot = rsfSpecComplete,
#        width = 360, height = 240, units = "mm", device = cairo_pdf)

return(specComplete)
