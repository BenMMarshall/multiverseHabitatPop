#' Generate effect plots
#'
#' @name generate_effect_plots
#' @description A
#' @param modelsList output tar_targets resulting in ssfBrms or areaBrms
#' @return a
#'
#' @export
generate_effect_plots <- function(modelsList){
  
  palette <- multiverseHabitat::get_palette()
  
  effectsPlotList <- lapply(names(modelsList), function(x){
    
    model <- modelsList[[x]]
    
    print(x)
    
    gradLimits <- model %>%
      spread_draws(`b_.*`, regex = TRUE) %>%
      dplyr::select(-.chain, -.iteration, -.draw) %>%
      tidyr::pivot_longer(cols = everything(), names_to = "variable") %>% 
      filter(!variable %in% c("b_Intercept"),
             !str_detect(variable, ":")) %>%
      summarise(
        minX = min(value),
        maxX = max(value)
      )
    
    modelEffectData <- model %>%
      spread_draws(`b_.*`, regex = TRUE) %>%
      dplyr::select(-.chain, -.iteration, -.draw) %>%
      tidyr::pivot_longer(cols = everything(), names_to = "variable") %>% 
      filter(!variable %in% c("b_Intercept"),
             !str_detect(variable, ":")) %>%
      mutate(
        variable = case_when(
          variable == "b_trackDuraScaled" ~ "\u03B2 Tracking Duration",
          variable == "b_trackFreqScaled" ~ "\u03B2 Tracking Frequency",
          variable == "b_sampleSizeScaled" ~ "\u03B2 Available Area: AKDE",
          variable == "b_modelFormulamf.ss" ~ "\u03B2 Model Formula: Not Integrated",
          variable == "b_stepDistgamma" ~ "\u03B2 Step Distribution: Gamma",
          variable == "b_turnDistvonmises" ~ "\u03B2 Turn Distribution: Von Mises",
          variable == "b_availablePerStepScaled" ~ "\u03B2 Available Area Contour",
          variable == "b_averagingMethodNaiveaverage" ~ "\u03B2 Available Points Per Step",
          variable == "b_areaMethodMCP" ~ "\u03B2 Available Area: MCP",
          variable == "b_samplingPatternst" ~ "\u03B2 Sampling Pattern: Stratified",
          variable == "b_typeIII" ~ "\u03B2 Desigen Type: III",
          variable == "b_contourScaled" ~ "\u03B2 Available Area Contour",
          variable == "b_availablePointsScaled" ~ "\u03B2 Available Points Multipiler",
          variable == "b_testrandomisation" ~ "\u03B2 Compana Test: Randomisation"
        )
      )
    
    # labelLocation <- c(-1.5, 1.15)
    # labelLocationText <- c(-1.45, 1.1)
    labelLocation <- data.frame(gradLimits)
    labelLocationText <- data.frame(gradLimits + c(0.05, -0.05))
    labelText <- c("Lower \npreference estimates",
                   "Higher \npreference estimates")
    gradAdj <- c(0, 0)
    
    (effectsPlot <- modelEffectData %>%
        ggplot(aes(x = value, y = variable)) +
        geom_vline(xintercept = 0, linewidth = 0.5, alpha = 0.9, colour = "#403F41",
                   linetype = 1) +
        tidybayes::stat_slab(aes(fill = after_stat(x)), fill_type = "gradient",
                             alpha = 0.85) +
        stat_pointinterval(position = position_dodge(width = 0.4, preserve = "single"),
                           point_interval = median_hdci, .width = c(.66, .95),
                           stroke = 1.25, colour = palette[c("coreGrey")]) +
        stat_summary(aes(colour = after_stat(x) > 0),
                     position = position_dodge(width = 0.2, preserve = "single"),
                     fun = median, size = 0.25) +
        scale_fill_gradient2(low = palette["BADGER"], mid = palette["coreGrey"], high = palette["2"],
                             # limits = c(gradLimits$minX+gradAdj[1],
                             #            gradLimits$maxX+gradAdj[2]),
                             # limits = c(-0.5, 0.35),
                             midpoint = 0,
                             oob = scales::oob_squish_any) +
        scale_colour_manual(values = unname(palette[c("BADGER", "2")])) +
        annotate("segment", x = -0.1, xend = labelLocation$minX, y = 0.65, yend = 0.65,
                 linewidth = 1.25,
                 arrow = arrow(angle = 30, type = "closed", length = unit(2, "mm")),
                 colour = palette["BADGER"]) +
        annotate("text", x = labelLocationText$minX, y = 0.85,
                 label = labelText[1],
                 colour = palette["BADGER"], hjust = 0, vjust = 0, lineheight = 0.95,
                 size = 3, fontface = 4) +
        annotate("segment", x = 0.1, xend = labelLocation$maxX, y = 0.65, yend = 0.65,
                 linewidth = 1.25,
                 arrow = arrow(angle = 30, type = "closed", length = unit(2, "mm")),
                 colour = palette["2"]) +
        annotate("text", x = labelLocationText$maxX, y = 0.85,
                 label = labelText[2],
                 colour = palette["2"], hjust = 1, vjust = 0, lineheight = 0.95,
                 size = 3, fontface = 4) +
        coord_cartesian(clip = "off", ylim = c(1, NA)) +
        labs(x = "Beta", y = "") +
        theme_bw() +
        theme(
          line = element_line(colour = palette["coreGrey"]),
          text = element_text(colour = palette["coreGrey"]),
          strip.background = element_blank(),
          strip.text = element_text(face = 4, hjust = 1, vjust = 1),
          # strip.text.y.left = element_text(angle = 0, margin = margin(-8,10,0,0)),
          # axis.text.y.left = element_text(margin = margin(0,-119,0,80)), # 2nd value needed to alligns with facet, 4th gives space left
          axis.title.x = element_text(margin = margin(5,0,0,0)),
          axis.ticks.y.left = element_blank(),
          axis.line.x = element_line(),
          strip.clip = "off",
          panel.border = element_blank(),
          panel.spacing = unit(18, "pt"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none")
    )
    
    ggsave(effectsPlot,
           filename = here("notebook", "figures", paste0(x, "_effectsPlot.png")),
           dpi = 300, width = 210, height = 120,
           units = "mm")
    
    print(paste0(x, " --- Done"))
    
  })
  
  return(effectsPlotList)
}