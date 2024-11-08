# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
library(tibble)
library(tidyselect)
library(readr)

# Set target options:
tar_option_set(
  packages = c("tibble",
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
               "ggtext",
               "patchwork",
               "bayesplot",
               "tidyr",
               "tidybayes"
  ), # packages that your targets need to run
  garbage_collection = TRUE,
  format = "qs", # storage format
  storage = "worker",
  retrieval = "worker",
  memory = "transient"
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# Run the R scripts in the R/ folder with your custom functions:
tar_source()
# source("other_functions.R") # Source other scripts as needed. # nolint


# Tibbles of splits -------------------------------------------------------

## every split requires it's own tibble, otherwise it's applied to all of them
values_SimSpecies_B <- tibble(
  # species = c("BADGER", "VULTURE", "KINGCOBRA")
  species = c("BADGER")
)
values_SimSpecies_V <- tibble(
  # species = c("BADGER", "VULTURE", "KINGCOBRA")
  species = c("VULTURE")
)
values_SimSpecies_K <- tibble(
  # species = c("BADGER", "VULTURE", "KINGCOBRA")
  species = c("KINGCOBRA")
)

values_SimIndi <- tibble(
  # individual = paste0("i", sprintf("%03d", 1:5))
  individual = paste0("i", sprintf("%03d", 1:25))
  # individual = paste0("i", 1:50)
  # individual = seq_len(30)
)

values_Regime <- tidyr::expand_grid(
  # td = c(7, 15, 30, 60, 120, 240),
  # tf = c(0.5, 1.0, 2.0, 6.0, 12.0, 24.0, 48.0, 168.0)
  td = c(7, 15, 30, 60),
  tf = c(1.0, 2.0, 6.0, 12.0)
)
# have to filter out certain combos that have too little data to work with
values_Regime <- values_Regime %>%
  dplyr::mutate(datapoints = td*24 * (1/tf)) %>%
  dplyr::filter(datapoints > 30) %>%
  dplyr::select(td, tf)

optionsList_area <- list(
  areaMethod = c("MCP", "AKDE"),
  areaContour = c(95, 99),
  # Method_ap = c(1,4),
  Method_ap = as.integer(round(exp(seq(log(1), log(10), length.out = 4)), digits = 1)),
  Method_sp = c("rd", "st"),
  Method_land = c("classRaster", "classRasterScram")
)

optionsList_areaMethods <- list(
  areaBasedMethod = c("Compana"),
  areaBasedTest = c("randomisation", "parametric")
)

optionsList_sff <- list(
  Method_method = c("ssf"),
  MethodSSF_land = c("classRaster", "classRasterScram"),
  # MethodSSF_as = c(2, 10),
  MethodSSF_as = as.integer(round(exp(seq(log(5), log(50), length.out = 5)), digits = 1)),
  # MethodSSF_mf = c("mf.is"),
  MethodSSF_mf = c("mf.is", "mf.ss"),
  # MethodSSF_sd = c("gamma"),
  # MethodSSF_td = c("vonmises")
  MethodSSF_td = c("vonmises", "unif"),
  MethodSSF_sd = c("gamma", "exp")
)

optionsList_pois <- list(
  MethodPois_land = c("classRaster", "classRasterScram"),
  # MethodPois_as = c(2, 10),
  MethodPois_as = as.integer(round(exp(seq(log(5), log(50), length.out = 5)), digits = 1)),
  # MethodPois_mf = c("mf.is"),
  MethodPois_mf = c("mf.is", "mf.ss"),
  # MethodPois_sd = c("gamma"),
  MethodPois_sd = c("gamma", "exp"),
  # MethodPois_td = c("vonmises")
  MethodPois_td = c("vonmises", "unif")
)


# Sampling set-up ---------------------------------------------------------

# repeats <- 1
repeats <- 2

# values_Sample <-
#   list(sampleSize = rep(c(3,5), each = repeats))
values_Sample <-
  list(sampleSize = rep(c(3,5,10,15,20), each = repeats))

set.seed(2023)

optionsList_samples <- lapply(values_Sample$sampleSize, function(x){
  sample(1:length(values_SimIndi$individual), x, replace = FALSE)
})

# optionsList_samples <- list(c(1,2),
#                             c(1,2,3))

names(optionsList_samples) <- paste0("samp", 1:length(optionsList_samples))

# remove the repeats of the full sample as they are all the same
# optionsList_samples <- optionsList_samples[1:(length(optionsList_samples)-(repeats-1))]

optionsCompleteList <- list(
  # "species" = values_SimSpecies,
  "individuals" = values_SimIndi,
  "regime" = values_Regime,
  "samples" = optionsList_samples,
  "repeats" = repeats,
  "pois" = optionsList_pois,
  "ssf" = optionsList_sff,
  "area" = optionsList_area,
  "areaMethods" = optionsList_areaMethods
)

saveRDS(optionsCompleteList, file = here::here("data", "optionsCompleteList.rds"))

# Targets workflow lists --------------------------------------------------

# BADGER SUB TREE ---------------------------------------------------------

individualSimulationsList_B <- list(
  ## LANDSCAPE SIMULATION
  individualSimulations_B <- tar_map(
    unlist = FALSE,
    values = values_SimSpecies_B,
    tar_target(landscape_B, simulate_landscape(species, 2023)), # FUNCTION simulate_landscape
    
    ## INDIDIVUAL SIMULATION
    tar_map(
      unlist = FALSE,
      values = values_SimIndi,
      tar_target(simData_B, simulate_individual(
        individualNum = individual,
        species = species,
        simSteps = 24*60 *365,
        desOptions = 12,
        options = 15,
        landscapeList = landscape_B,
        seed = 2023),
        priority = 0.93)#, # FUNCTION simulate_individual
    ) # INDI SIM
  ) # LANDSCAPE SIM
  
)

allIndividualsList_B <- list(
  tar_combine(
    allIndividuals_B,
    use_names = FALSE,
    individualSimulations_B,
    # command = dplyr::bind_rows(!!!.x))
    command = list(!!!.x))
)

# coreMultiverse_B <- list(
#   tar_map(
#     unlist = TRUE,
#     values = values_Regime,
#     tar_target(sampDuraFreqData_B,
#                subset_duration(
#                  allIndividualData = subset_frequency(
#                    allIndividualData = allIndividuals_B,
#                    freqPreset = tf),
#                  daysDuration = td),
#                priority = 0.92),
#     tar_target(areaBasedAvailUse_B,
#                area_based_extraction(
#                  allIndividualData = sampDuraFreqData_B,
#                  optionsList = optionsList_area
#                ),
#                priority = 0.9),
#     tar_target(areaBasedOUT_B,
#                area_based_calculations(
#                  availUseData = areaBasedAvailUse_B,
#                  sampleGroups = optionsList_samples,
#                  optionsList = optionsList_area,
#                  optionsListArea = optionsList_areaMethods
#                ),
#                priority = 0.9),
#     # tar_target(ssfOUT,
#     #            wrapper_indi_ssf(
#     #              allIndividualData = sampDuraFreqData,
#     #              optionsList = optionsList_sff
#     #            ),
#     #            priority = 0.9),
#     tar_target(ssfSampled_B,
#                sample_ssf_results(
#                  sampDuraFreqData_B,
#                  sampleGroups = optionsList_samples,
#                  optionsList = optionsList_sff
#                ),
#                priority = 0.9),
#     tar_target(poisOUT_B,
#                method_pois_inla(
#                  allIndividualData = sampDuraFreqData_B,
#                  sampleGroups = optionsList_samples,
#                  optionsList = optionsList_pois),
#                priority = 0.9),
#     tar_target(twoStepOUT_B,
#                method_twoStep(
#                  allIndividualData = sampDuraFreqData_B,
#                  sampleGroups = optionsList_samples,
#                  optionsList = optionsList_pois),
#                priority = 0.9)
#   )
# )

# VULTURE SUB TREE --------------------------------------------------------

individualSimulationsList_V <- list(
  ## LANDSCAPE SIMULATION
  individualSimulations_V <- tar_map(
    unlist = FALSE,
    values = values_SimSpecies_V,
    tar_target(landscape_V, simulate_landscape(species, 2023)), # FUNCTION simulate_landscape
    
    ## INDIDIVUAL SIMULATION
    tar_map(
      unlist = FALSE,
      values = values_SimIndi,
      tar_target(simData_V, simulate_individual(
        individualNum = individual,
        species = species,
        simSteps = 24*60 *365,
        desOptions = 12,
        options = 15,
        landscapeList = landscape_V,
        seed = 2023),
        priority = 0.93)#, # FUNCTION simulate_individual
    ) # INDI SIM
  ) # LANDSCAPE SIM
  
)

allIndividualsList_V <- list(
  tar_combine(
    allIndividuals_V,
    use_names = FALSE,
    individualSimulations_V,
    # command = dplyr::bind_rows(!!!.x))
    command = list(!!!.x))
)

# KING COBRA SUB TREE -----------------------------------------------------

individualSimulationsList_K <- list(
  ## LANDSCAPE SIMULATION
  individualSimulations_K <- tar_map(
    unlist = FALSE,
    values = values_SimSpecies_K,
    tar_target(landscape_K, simulate_landscape(species, 2023)), # FUNCTION simulate_landscape
    
    ## INDIDIVUAL SIMULATION
    tar_map(
      unlist = FALSE,
      values = values_SimIndi,
      tar_target(simData_K, simulate_individual(
        individualNum = individual,
        species = species,
        simSteps = 24*60 *365,
        desOptions = 12,
        options = 15,
        landscapeList = landscape_K,
        seed = 2023),
        priority = 0.93)#, # FUNCTION simulate_individual
    ) # INDI SIM
  ) # LANDSCAPE SIM
  
)

allIndividualsList_K <- list(
  tar_combine(
    allIndividuals_K,
    use_names = FALSE,
    individualSimulations_K,
    # command = dplyr::bind_rows(!!!.x))
    command = list(!!!.x))
)


# CORE MULTIVERSE ANALYSIS ------------------------------------------------

coreMultiverse <- list(
  tar_map(
    unlist = TRUE,
    values = values_Regime,
    ## BADGER
    tar_target(sampDuraFreqData_B,
               subset_duration(
                 allIndividualData = subset_frequency(
                   allIndividualData = allIndividuals_B,
                   freqPreset = tf),
                 daysDuration = td),
               priority = 0.92),
    tar_target(areaBasedAvailUse_B,
               area_based_extraction(
                 allIndividualData = sampDuraFreqData_B,
                 optionsList = optionsList_area
               ),
               priority = 0.9),
    tar_target(areaBasedOUT_B,
               area_based_calculations(
                 availUseData = areaBasedAvailUse_B,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_area,
                 optionsListArea = optionsList_areaMethods
               ),
               priority = 0.9),
    # tar_target(ssfOUT,
    #            wrapper_indi_ssf(
    #              allIndividualData = sampDuraFreqData,
    #              optionsList = optionsList_sff
    #            ),
    #            priority = 0.9),
    tar_target(ssfSampled_B,
               sample_ssf_results(
                 sampDuraFreqData_B,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_sff
               ),
               priority = 0.9),
    tar_target(poisOUT_B,
               method_pois_inla(
                 allIndividualData = sampDuraFreqData_B,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9),
    tar_target(twoStepOUT_B,
               method_twoStep(
                 allIndividualData = sampDuraFreqData_B,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9),
    ## VULTURE
    tar_target(sampDuraFreqData_V,
               subset_duration(
                 allIndividualData = subset_frequency(
                   allIndividualData = allIndividuals_V,
                   freqPreset = tf),
                 daysDuration = td),
               priority = 0.92),
    tar_target(areaBasedAvailUse_V,
               area_based_extraction(
                 allIndividualData = sampDuraFreqData_V,
                 optionsList = optionsList_area
               ),
               priority = 0.9),
    tar_target(areaBasedOUT_V,
               area_based_calculations(
                 availUseData = areaBasedAvailUse_V,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_area,
                 optionsListArea = optionsList_areaMethods
               ),
               priority = 0.9),
    tar_target(ssfSampled_V,
               sample_ssf_results(
                 sampDuraFreqData_V,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_sff
               ),
               priority = 0.9),
    tar_target(poisOUT_V,
               method_pois_inla(
                 allIndividualData = sampDuraFreqData_V,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9),
    tar_target(twoStepOUT_V,
               method_twoStep(
                 allIndividualData = sampDuraFreqData_V,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9),
    ## KING COBRA
    tar_target(sampDuraFreqData_K,
               subset_duration(
                 allIndividualData = subset_frequency(
                   allIndividualData = allIndividuals_K,
                   freqPreset = tf),
                 daysDuration = td),
               priority = 0.92),
    tar_target(areaBasedAvailUse_K,
               area_based_extraction(
                 allIndividualData = sampDuraFreqData_K,
                 optionsList = optionsList_area
               ),
               priority = 0.9),
    tar_target(areaBasedOUT_K,
               area_based_calculations(
                 availUseData = areaBasedAvailUse_K,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_area,
                 optionsListArea = optionsList_areaMethods
               ),
               priority = 0.9),
    tar_target(ssfSampled_K,
               sample_ssf_results(
                 sampDuraFreqData_K,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_sff
               ),
               priority = 0.9),
    tar_target(poisOUT_K,
               method_pois_inla(
                 allIndividualData = sampDuraFreqData_K,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9),
    tar_target(twoStepOUT_K,
               method_twoStep(
                 allIndividualData = sampDuraFreqData_K,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9)
  )
)

# Combining all results ---------------------------------------------------

ssfCompiled <- list(
  tar_combine(
    ssfResults,
    coreMultiverse[[1]][grep("ssfSampled", names(coreMultiverse[[1]]))],
    command = rbind(!!!.x),
    priority = 0.8
  ),
  tar_target(
    ssfEstimateOutputs,
    write.csv(ssfResults,
              here::here("data", "ssfEstimateOutputs_uncom.csv"), row.names = FALSE),
    format = "file"
  ),
  tar_target(
    ssfSpecCurve,
    generate_spec_curves(
      outputResults = ssfResults,
      method = "ssf"
    )
  ),
  # tar_target(
  #   ssfPresentCurve,
  #   generate_presentation_curves(
  #     outputResults = ssfResults,
  #     method = "ssf"
  #   )
  # ),
  tar_target(
    ssfBrms,
    run_brms(
      resultsData = ssfResults,
      iter = 4000,
      warmup = 750,
      thin = 4
    )
  )
)

poisCompiled <- list(
  tar_combine(
    poisResults,
    coreMultiverse[[1]][grep("poisOUT", names(coreMultiverse[[1]]))],
    command = rbind(!!!.x),
    priority = 0.8
  ),
  tar_target(
    poisEstimateOutputs,
    write.csv(poisResults,
              here::here("data", "poisEstimateOutputs_uncom.csv"), row.names = FALSE),
    format = "file"
  ),
  tar_target(
    poisSpecCurve,
    generate_spec_curves(
      outputResults = poisResults,
      method = "pois"
    )
  ),
  # tar_target(
  #   poisPresentCurve,
  #   generate_presentation_curves(
  #     outputResults = poisResults,
  #     method = "pois"
  #   )
  # ),
  tar_target(
    poisBrms,
    run_brms(
      resultsData = poisResults,
      iter = 20000,
      warmup = 8000,
      thin = 20
    )
  )
)

twoStepCompiled <- list(
  tar_combine(
    twoStepResults,
    coreMultiverse[[1]][grep("twoStep", names(coreMultiverse[[1]]))],
    command = rbind(!!!.x),
    priority = 0.8
  ),
  tar_target(
    twoStepEstimateOutputs,
    write.csv(twoStepResults,
              here::here("data", "twoStepEstimateOutputs_uncom.csv"), row.names = FALSE),
    format = "file"
  ),
  tar_target(
    twoStepSpecCurve,
    generate_spec_curves(
      outputResults = twoStepResults,
      method = "twoStep"
    )
  ),
  tar_target(
    twoStepBrms,
    run_brms(
      resultsData = twoStepResults,
      iter = 4000,
      warmup = 750,
      thin = 4
    )
  )
)

areaBasedCompiled <- list(
  tar_combine(
    areaBasedResults,
    coreMultiverse[[1]][grep("areaBasedOUT", names(coreMultiverse[[1]]))],
    command = rbind(!!!.x),
    priority = 0.8
  ),
  tar_target(
    areaBasedEstimateOutputs,
    write.csv(areaBasedResults,
              here::here("data", "areaBasedEstimateOutputs_uncom.csv"), row.names = FALSE),
    format = "file"
  ),
  tar_target(
    areaSpecCurve,
    generate_spec_curves(
      outputResults = areaBasedResults,
      method = "area"
    )
  ),
  tar_target(
    areaBrms,
    run_brms(
      resultsData = areaBasedResults,
      iter = 4000,
      warmup = 750,
      thin = 4
    )
  )
)

brmModelOutputs <- list(
  tar_combine(
    modelsBrms,
    # manually pull out the brms model outputs
    list(
      ssfCompiled[[4]],
      areaBasedCompiled[[4]],
      poisCompiled[[4]],
      twoStepCompiled[[4]]),
    command = list(!!!.x),
    priority = 0.5
  ),
  tar_target(
    diagnosticPlots,
    diagnostics_brms(modelsList = modelsBrms),
    priority = 0.4
  ),
  tar_target(
    modelExtracts,
    extract_model_values(modelsList = modelsBrms),
    priority = 0.4
  ),
  tar_target(
    effectPlots,
    generate_effect_plots(modelsList = modelsBrms),
    priority = 0.4
  ),
  tar_target(
    allEffectPlots,
    generate_allEffect_plots(modelExtracts = modelExtracts),
    priority = 0.4
  ),
  tar_target(
    rmdRender,
    render_rmd(modelExtracts,
               effectPlots,
               allEffectPlots,
               areaBasedEstimateOutputs,
               areaSpecCurve,
               twoStepSpecCurve,
               poisSpecCurve,
               ssfSpecCurve),
    cue = tar_cue(mode = "always"),
    priority = 0.01
  )
)

# odd issue, that targets have problem with writing to .gz compression during
# tar_targets, figures sidestepping this with a post compression is ok. Added
# uncompressed files to .gitignore
if(all(list.files(here::here("data"), pattern = "EstimateOutputs_uncom.csv") %in% 
       c("areaBasedEstimateOutputs_uncom.csv",
         "poisEstimateOutputs_uncom.csv", 
         "ssfEstimateOutputs_uncom.csv",
         "twoStepEstimateOutputs_uncom.csv"))){
  # zip(zipfile = here::here("data", "EstimateOutputs"),
  #     files = list.files(here::here("data"), pattern = "EstimateOutputs_uncom.csv",
  #                        full.names = TRUE))
  for(file in list.files(here::here("data"), pattern = "EstimateOutputs_uncom.csv")){
    # file <- "areaBasedEstimateOutputs_uncom.csv"
    fileCom <- sub("_uncom.csv$", ".csv.gz", file)
    write_csv(read_csv(here::here("data", file)),
              here::here("data", fileCom))
  }
}

# All targets lists -------------------------------------------------------

list(individualSimulationsList_B,
     individualSimulationsList_V,
     individualSimulationsList_K,
     allIndividualsList_B,
     allIndividualsList_V,
     allIndividualsList_K,
     coreMultiverse,
     ssfCompiled,
     poisCompiled,
     twoStepCompiled,
     areaBasedCompiled,
     brmModelOutputs
)