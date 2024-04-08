#' Generate presentation curves
#'
#' @name generate_presentation_curves
#' @description A
#' @param outputResults 
#' @param method c("ssf", "area", "pois", "twoStep")
#' @return a
#'
#' @export
generate_presentation_curves <- function(outputResults, method){
  
  palette <- multiverseHabitat::get_palette()
  # targets::tar_load("ssfResults")
  # outputResults <- ssfResults
  # method <- "ssf"
  
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
      unique(outputResults$sampleID),
      sort(c(unique(outputResults$sampleSize),
             unique(outputResults$trackFreq),
             unique(outputResults$trackDura),
             unique(outputResults$availablePerStep))),
      sort(unique(outputResults$modelFormula)),
      sort(unique(outputResults$stepDist)),
      sort(unique(outputResults$turnDist)),
      sort(unique(outputResults$averagingMethod))
    ))
    
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
      mutate("estimate" = modelAvg) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(variable = factor(variable, levels = c(
        "Tracking Duration (days)",
        "Tracking Frequency (points/hour)",
        "Sample Size (n)",
        "Model Formula (SSF or iSSF)",
        "Available Points per Step",
        "Distribution of Step Lengths",
        "Distribution of Turn Angles",
        "Model Averaging Method")
      ))
    
  } else if(method == "area"){
    return(NULL)
    
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
      sort(c(unique(outputResults$sampleSize),
             unique(outputResults$trackFreq),
             unique(outputResults$trackDura),
             unique(outputResults$availablePerStep))),
      sort(unique(outputResults$modelFormula)),
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
      mutate("estimate" = prefDiff) %>% 
      dplyr::ungroup() %>%
      dplyr::mutate(variable = factor(variable, levels = c(
        "Tracking Duration (days)",
        "Tracking Frequency (points/hour)",
        "Sample Size (n)",
        "Model Formula (SSF or iSSF)",
        "Available Points per Step",
        "Distribution of Step Lengths",
        "Distribution of Turn Angles")
      ))
    
  } else if(method == "twoStep"){
    return(NULL)
  }
  
  # PLOTS -------------------------------------------------------------------
  
  outputPlotData$classLandscape <- ifelse(str_detect(outputPlotData$classLandscape, "Scram"),
                                          "Scrambled Habitat Layer (i.e., No selection)",
                                          "Correct Habitat Layer (i.e., Positive selection)")
  outputResults$classLandscape <- ifelse(str_detect(outputResults$classLandscape, "Scram"),
                                         "Scrambled Habitat Layer (i.e., No selection)",
                                         "Correct Habitat Layer (i.e., Positive selection)")
  
  if(method %in% c("twoStep", "area")){
    xlimits <- as.vector(quantile(outputPlotData$estimate, probs = c(.001, .999)))
    upperOutliers <- outputResults %>% 
      group_by(classLandscape) %>% 
      filter(estimate < xlimits[1]) %>% 
      count()
    lowerOutliers <- outputResults %>% 
      group_by(classLandscape) %>% 
      filter(estimate > xlimits[2]) %>% 
      count()
  } else {
    xlimits <- c(NA, NA)
  }
  
  plotTypes <- c("trackDura7", "trackDura60", "modelFormI", "modelFormS")
  
  plotList <- vector("list", length = length(plotTypes))
  names(plotList) <- plotTypes
  for(pType in plotTypes){
    # pType <- plotTypes[1]
    overallMed_present <- outputPlotData %>%
      filter(classLandscape == "Correct Habitat Layer (i.e., Positive selection)") %>% 
      filter(variable == "Tracking Duration (days)") %>% 
      mutate(n = n()) %>%
      summarise(medEst = median(estimate, na.rm = TRUE),
                n = n[1])
    
    if(str_detect(pType, "Dura")){
      outputPlotData_present <- outputPlotData %>%
        filter(classLandscape == "Correct Habitat Layer (i.e., Positive selection)") %>% 
        filter(variable == "Tracking Duration (days)")
      
      if(str_detect(pType, "7")){
        outputPlotData_present <- outputPlotData_present %>% 
          mutate(value = factor(case_when(
            value == "7" ~ "7",
            TRUE ~ "Remaining"
          ), levels = c("Remaining", "7")
          ))
      } else if(str_detect(pType, "60")){
        outputPlotData_present <- outputPlotData_present %>% 
          mutate(value = factor(case_when(
            value == "60" ~ "60",
            TRUE ~ "Remaining"
          ), levels = c("Remaining", "60")
          ))
      }
      
    } else if(str_detect(pType, "Form")){
      outputPlotData_present <- outputPlotData %>%
        filter(classLandscape == "Correct Habitat Layer (i.e., Positive selection)") %>% 
        filter(variable == "Model Formula (SSF or iSSF)")
    }
    
    presentationFigure <- outputPlotData_present %>% 
      ggplot() +
      geom_vline(xintercept = 0, linewidth = 0.5, alpha = 0.9, colour = "#403F41",
                 linetype = 1) +
      geom_point(aes(x = estimate, y = value, colour = d_medEst),
                 position = position_jitter(width = 0, height = 0.45), alpha = 0.25,
                 pch = 3, size = 0.75) +
      # geom_point(data = medData, aes(x = modelMedEst, y = value),
      #            alpha = 1, size = 1.5, colour = "#FFFFFF") +
      # geom_point(data = medData, aes(x = modelMedEst, y = value),
      #            alpha = 1, size = 1, colour = "#403F41") +
      geom_segment(data = overallMed_present,
                   aes(x = medEst, xend = medEst, y = Inf,
                       yend = -Inf),
                   alpha = 1, linewidth = 0.45, linetype = 1,
                   colour = palette["BADGER"]) +
      geom_hline(yintercept = seq(0.5,10.5,1), linewidth = 0.5, alpha = 0.25, colour = "#403F41",
                 linetype = 2) +
      facet_grid(variable~classLandscape, scales = "free_y", space = "free", switch = "y") +
      labs(y = "", x = "Estimate") +
      scale_x_continuous(limits = xlimits) +
      scale_colour_gradient2(low = palette["BADGER"], mid = palette["coreGrey"], high = palette["2"]) +
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
    
    ggsave(filename = here("notebook", "figures",
                           paste0("presentationDraft_", method, "_", pType, ".pdf")),
           plot = presentationFigure,
           width = 360, height = 240, units = "mm", device = cairo_pdf)
    
    plotList[[pType]] <-  presentationFigure
  }
  
  return(plotList)
  
}

