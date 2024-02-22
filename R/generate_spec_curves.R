#' Generate specification curves
#'
#' @name generate_spec_curves
#' @description A
#' @param outputResults 
#' @param method c("ssf", "area", "pois")
#' @return a
#'
#' @export
generate_spec_curves <- function(outputResults, method){
  
  palette <- multiverseHabitat::get_palette()
  
  if(method == "ssf"){
    
    outputResults <- outputResults %>% 
      mutate("estimate" = modelAvg)
    
    # convert tf to points/hour to help interpretation
    outputResults$trackFreq <- 1/as.numeric(outputResults$trackFreq)
    outputResults$trackFreq <- round(outputResults$trackFreq, digits = 2)
    
    outputResults$stepDist <- ifelse(outputResults$stepDist == "gamma",
                                     "Gamma",
                                     "Exponential")
    outputResults$turnDist <- ifelse(outputResults$turnDist == "unif",
                                     "Uniform",
                                     "Von Mises")
    outputResults$modelFormula <- ifelse(outputResults$modelForm == "mf.is",
                                         "Integrated",
                                         "Standard")
    
    levelOrdering <- unique(c(
      unique(outputResults$averagingMethod),
      unique(outputResults$sampleID),
      sort(unique(outputResults$sampleSize)),
      sort(unique(outputResults$trackFreq)),
      sort(unique(outputResults$trackDura)),
      sort(unique(outputResults$modelFormula)),
      sort(unique(outputResults$availablePerStep)),
      sort(unique(outputResults$stepDist)),
      sort(unique(outputResults$turnDist))))
    
    outputPlotData <- outputResults %>% 
      dplyr::select(-analysis, -sampleID) %>% 
      dplyr::mutate(across(c(1:3, 5:9), as.character)) %>% 
      tidyr::pivot_longer(cols = c(1:3, 5:9), names_to = "variable") %>% 
      dplyr::mutate(
        variable = case_when(
          variable == "trackDura" ~ "Tracking Duration (days)",
          variable == "trackFreq" ~ "Tracking Frequency (points/hour)",
          variable == "sampleSize" ~ "Sample Size (n)",
          variable == "modelFormula" ~ "Model Formula (SSF or iSSF)",
          variable == "availablePerStep" ~ "Available Points per Step",
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
    
  } else if(method == "area"){
    
    # Area based --------------------------------------------------------------

    outputResults <- outputResults %>% 
      mutate("estimate" = companaHabDiff)
    
    outputResults$trackFreq <- 1/as.numeric(outputResults$trackFreq)
    outputResults$trackFreq <- round(outputResults$trackFreq, digits = 2)
    
    outputResults$samplingPattern <- ifelse(outputResults$samplingPattern == "rd",
                                            "Random",
                                            "Stratified")
    outputResults$test <- ifelse(outputResults$test == "randomisation",
                                            "Randomisation",
                                            "Parametric")
    
    levelOrdering <- unique(c(
      unique(outputResults$sampleID),
      sort(unique(outputResults$sampleSize)),
      sort(unique(outputResults$trackFreq)),
      sort(unique(outputResults$trackDura)),
      unique(outputResults$type),
      unique(outputResults$areaMethod),
      sort(unique(outputResults$contour)),
      sort(unique(outputResults$availablePoints)),
      sort(unique(outputResults$samplingPattern)),
      sort(unique(outputResults$test))))
    
    outputPlotData <- outputResults %>% 
      dplyr::select(-sampleID) %>% 
      dplyr::mutate(across(c(1:3, 5:10), as.character)) %>% 
      tidyr::pivot_longer(cols = c(1:3, 5:10), names_to = "variable") %>% 
      dplyr::mutate(
        variable = case_when(
          variable == "trackDura" ~ "Tracking Duration (days)",
          variable == "trackFreq" ~ "Tracking Frequency (points/hour)",
          variable == "sampleSize" ~ "Sample Size (n)",
          variable == "type" ~ "Type II or Type III Habitat Selection",
          variable == "availablePoints" ~ "Available Points Multiplier",
          variable == "samplingPattern" ~ "Sampling Pattern",
          variable == "test" ~ "Compana Test Method",
          variable == "contour" ~ "Available Area Contour (%)",
          variable == "areaMethod" ~ "Available Area Method"
        ),
        # sampleID = as.factor(sampleID),
        value = factor(value, levels = levelOrdering)) %>%
      dplyr::group_by(variable, value) %>%
      dplyr::mutate(d_medEst = companaHabDiff - median(outputResults$companaHabDiff, na.rm = TRUE)) %>%
      dplyr::ungroup() %>% 
      mutate("estimate" = companaHabDiff)
    
  } else if(method == "pois"){
    
    # Poisson models ----------------------------------------------------------
    
    # outputResults <- poisResults
    
    prefDiffDF <- outputResults %>% 
      mutate(key = paste0(sampleID, trackFreq, trackDura, modelFormula, availablePerStep, stepDist,
                          turnDist, classLandscape)) %>% 
      group_by(key) %>% 
      summarise(prefDiff = diff(mean))
    
    outputResults <- outputResults %>% 
      filter(term == "layerc2") %>% 
      mutate(key = paste0(sampleID, trackFreq, trackDura, modelFormula, availablePerStep, stepDist,
                          turnDist, classLandscape)) %>% 
      left_join(prefDiffDF) %>% 
      dplyr::mutate(medEst = median(mean),
                    absDeltaEst = mean - medEst) %>% 
      dplyr::mutate(
        sampleSizeScaled = (sampleSize - mean(sampleSize))/sd(sampleSize),
        trackFreqScaled = (trackFreq - mean(trackFreq))/sd(trackFreq),
        trackDuraScaled = (trackDura - mean(trackDura))/sd(trackDura),
        availablePerStepScaled  = (availablePerStep - mean(availablePerStep))/sd(availablePerStep)
      )
    
    outputResults <- outputResults %>% 
      mutate("estimate" = prefDiff)
    
    outputResults$trackFreq <- 1/as.numeric(outputResults$trackFreq)
    outputResults$trackFreq <- round(outputResults$trackFreq, digits = 2)
    
    outputResults$stepDist <- ifelse(outputResults$stepDist == "gamma",
                                     "Gamma",
                                     "Exponential")
    outputResults$turnDist <- ifelse(outputResults$turnDist == "unif",
                                     "Uniform",
                                     "Von Mises")
    outputResults$modelFormula <- ifelse(outputResults$modelForm == "mf.is",
                                         "Integrated",
                                         "Standard")
    
    levelOrdering <- unique(c(
      unique(outputResults$sampleID),
      sort(unique(outputResults$sampleSize)),
      sort(unique(outputResults$trackFreq)),
      sort(unique(outputResults$trackDura)),
      sort(unique(outputResults$modelFormula)),
      sort(unique(outputResults$availablePerStep)),
      sort(unique(outputResults$stepDist)),
      sort(unique(outputResults$turnDist))))
    
    outputPlotData <- outputResults %>% 
      dplyr::select(-analysis, -sampleID) %>% 
      dplyr::mutate(across(c(1:3, 5:8), as.character)) %>% 
      tidyr::pivot_longer(cols = c(1:3, 5:8), names_to = "variable") %>% 
      dplyr::mutate(
        variable = case_when(
          variable == "trackDura" ~ "Tracking Duration (days)",
          variable == "trackFreq" ~ "Tracking Frequency (points/hour)",
          variable == "sampleSize" ~ "Sample Size (n)",
          variable == "modelFormula" ~ "Model Formula (SSF or iSSF)",
          variable == "availablePerStep" ~ "Available Points per Step",
          variable == "stepDist" ~ "Distribution of Step Lengths",
          variable == "turnDist" ~ "Distribution of Turn Angles"
        ),
        # sampleID = as.factor(sampleID),
        value = factor(value, levels = levelOrdering)) %>%
      dplyr::group_by(variable, value) %>%
      dplyr::mutate(d_medEst = prefDiff - median(outputResults$prefDiff, na.rm = TRUE)) %>%
      dplyr::ungroup() %>% 
      mutate("estimate" = prefDiff)
    
  } else if(method == "twoStep"){
    
    # TwoStep models ----------------------------------------------------------
    
    outputResults <- outputResults %>% 
      mutate("estimate" = twoStepBeta)
    
    outputResults$trackFreq <- 1/as.numeric(outputResults$trackFreq)
    outputResults$trackFreq <- round(outputResults$trackFreq, digits = 2)
    
    outputResults$stepDist <- ifelse(outputResults$stepDist == "gamma",
                                     "Gamma",
                                     "Exponential")
    outputResults$turnDist <- ifelse(outputResults$turnDist == "unif",
                                     "Uniform",
                                     "Von Mises")
    outputResults$modelFormula <- ifelse(outputResults$modelForm == "mf.is",
                                         "Integrated",
                                         "Standard")
    
    levelOrdering <- unique(c(
      unique(outputResults$sampleID),
      sort(unique(outputResults$sampleSize)),
      sort(unique(outputResults$trackFreq)),
      sort(unique(outputResults$trackDura)),
      sort(unique(outputResults$modelFormula)),
      sort(unique(outputResults$availablePerStep)),
      sort(unique(outputResults$stepDist)),
      sort(unique(outputResults$turnDist))))
    
    outputPlotData <- outputResults %>% 
      dplyr::select(-analysis, -sampleID) %>% 
      dplyr::mutate(across(c(1:3, 5:8), as.character)) %>% 
      tidyr::pivot_longer(cols = c(1:3, 5:8), names_to = "variable") %>% 
      dplyr::mutate(
        variable = case_when(
          variable == "trackDura" ~ "Tracking Duration (days)",
          variable == "trackFreq" ~ "Tracking Frequency (points/hour)",
          variable == "sampleSize" ~ "Sample Size (n)",
          variable == "modelFormula" ~ "Model Formula (SSF or iSSF)",
          variable == "availablePerStep" ~ "Available Points per Step",
          variable == "stepDist" ~ "Distribution of Step Lengths",
          variable == "turnDist" ~ "Distribution of Turn Angles"
        ),
        # sampleID = as.factor(sampleID),
        value = factor(value, levels = levelOrdering)) %>%
      dplyr::group_by(variable, value) %>%
      dplyr::mutate(d_medEst = twoStepBeta - median(outputResults$twoStepBeta, na.rm = TRUE)) %>%
      dplyr::ungroup() %>% 
      mutate("estimate" = twoStepBeta)
    
  }
  
  # PLOTS -------------------------------------------------------------------
  
  outputPlotData$classLandscape <- ifelse(str_detect(outputPlotData$classLandscape, "Scram"),
                                          "Scrambled Habitat Layer (i.e., No selection)",
                                          "Correct Habitat Layer (i.e., Positive selection)")
  outputResults$classLandscape <- ifelse(str_detect(outputResults$classLandscape, "Scram"),
                                          "Scrambled Habitat Layer (i.e., No selection)",
                                          "Correct Habitat Layer (i.e., Positive selection)")
  
  medData <- outputPlotData %>%
    dplyr::group_by(variable, value, classLandscape) %>%
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
      facet_grid(variable~classLandscape, scales = "free_y", space = "free", switch = "y") +
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
  
  overallMed <- outputResults %>%
    mutate(n = n()) %>%
    group_by(classLandscape) %>%
    summarise(medEst = median(estimate, na.rm = TRUE),
              n = n[1])
  
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
      geom_segment(data = overallMed,
                   aes(x = medEst, xend = medEst, y = Inf,
                       yend = -Inf),
                   alpha = 1, linewidth = 0.45, linetype = 1,
                   colour = palette["BADGER"]) +
      geom_text(data = overallMed, aes(x = medEst, y = 0,
                                       label = paste0(" Median = ",
                                                      round(medEst, digits = 2))),
                hjust = 0, vjust = 0, fontface = 4, colour = palette["BADGER"]) +
      # annotate("text", x = overallMed$medEst + (overallMed$medEst + max(outputResults$estimate))/2,
      #          y = overallMed$indexLoc, label = "Median",
      #          fontface = 4, size = 5, colour = palette["coreGrey"],
      #          hjust = 1, vjust = -0.2) +
      # annotate("segment", x = overallMed$medEst + (overallMed$medEst + max(outputResults$estimate))/2,
      #          xend = overallMed$medEst,
      #          y = overallMed$indexLoc, yend = overallMed$indexLoc,
      #          linewidth = 0.75, colour = palette["coreGrey"]) +
      # scale_colour_gradient2(low = palette["BADGER"], mid = palette["coreGrey"], high = palette["2"]) +
      facet_grid(.~classLandscape, space = "free", switch = "y") +
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
  
  (specComplete <- wrap_plots(overallSpecCurve, splitSpecCurve) +
      plot_layout(heights = c(1, 3), guides = "collect"))
  
  # (specComplete <- overallSpecCurve / splitSpecCurve +
  #     plot_layout(heights = c(1, 3), guides = "collect"))
  
  ggsave(filename = here("notebook", "figures",
                         paste0(method, "_specCurve.png")),
         plot = specComplete,
         width = 360, height = 240, units = "mm", dpi = 300)
  # ggsave(filename = here("notebook", "figures",
  #                        paste0(method, "_specCurve.pdf")),
  #        plot = specComplete,
  #        width = 360, height = 240, units = "mm", device = cairo_pdf)
  
  return(specComplete)
  
}