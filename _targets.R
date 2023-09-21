# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
library(tarchetypes) # Load other packages as needed. # nolint
library(tibble)
library(tidyselect)

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
               "patchwork",
               "bayesplot",
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
values_SimSpecies <- tibble(
  # species = c("BADGER", "VULTURE", "KINGCOBRA")
  species = c("BADGER")
)
values_SimIndi <- tibble(
  individual = paste0("i", sprintf("%03d", 1:3))
  # individual = paste0("i", 1:50)
  # individual = seq_len(30)
)

values_Regime <- tidyr::expand_grid(
  # td = c(7, 15, 30, 60, 120, 240),
  # tf = c(0.5, 1.0, 2.0, 6.0, 12.0, 24.0, 48.0, 168.0)
  td = c(15, 60),
  tf = c(1, 6)
)
# have to filter out certain combos that have too little data to work with
values_Regime <- values_Regime %>%
  dplyr::mutate(datapoints = td*24 * (1/tf)) %>%
  dplyr::filter(datapoints > 30) %>%
  dplyr::select(td, tf)


optionsList_area <- list(
  areaMethod = c("MCP", "AKDE"),
  areaContour = c(95, 99),
  Method_ap = as.integer(round(exp(seq(log(1), log(10), length.out = 2)), digits = 1)),
  Method_sp = c("rd", "st")
)

optionsList_areaMethods <- list(
  areaBasedMethod = c("Compana"),
  areaBasedTest = c("randomisation", "parametric")
)

optionsList_sff <- list(
  Method_method = c("ssf"),
  # MethodSSF_mf = c("mf.is"),
  # MethodSSF_sd = c("gamma"),
  # MethodSSF_td = c("vonmises"),
  MethodSSF_as = c(2, 10),
  MethodSSF_mf = c("mf.is", "mf.ss"),
  MethodSSF_sd = c("gamma", "exp"),
  MethodSSF_td = c("vonmises", "unif")
  # MethodSSF_as = as.integer(round(exp(seq(log(5), log(500), length.out = 5)), digits = 1))
)

optionsList_pois <- list(
  # MethodPois_mf = c("mf.is"),
  # MethodPois_sd = c("gamma"),
  # MethodPois_td = c("vonmises"),
  MethodPois_as = c(2, 10),
  MethodPois_mf = c("mf.is", "mf.ss"),
  MethodPois_sd = c("gamma", "exp"),
  MethodPois_td = c("vonmises", "unif")
)


# Sampling set-up ---------------------------------------------------------

# values_Sample <-
#   list(sampleSize = rep(c(5,10,15,25,45), each = 3))
# 
# set.seed(1)
# 
# optionsList_samples <- lapply(values_Sample$sampleSize, function(x){
#   sample(1:50, x, replace = FALSE)
# })

optionsList_samples <- list(c(1,2),
                            c(1,2,3))

names(optionsList_samples) <- paste0("samp", 1:length(optionsList_samples))

# Targets workflow lists --------------------------------------------------

individualSimulationsList <- list(
  ## LANDSCAPE SIMULATION
  individualSimulations <- tar_map(
    unlist = FALSE,
    values = values_SimSpecies,
    tar_target(landscape, simulate_landscape(species, 2023)), # FUNCTION simulate_landscape
    
    ## INDIDIVUAL SIMULATION
    tar_map(
      unlist = FALSE,
      values = values_SimIndi,
      tar_target(simData, simulate_individual(
        individualNum = individual,
        species = species,
        simSteps = 24*60 *365,
        desOptions = 12,
        options = 15,
        landscapeList = landscape,
        seed = 2023),
        priority = 0.93)#, # FUNCTION simulate_individual
    ) # INDI SIM
  ) # LANDSCAPE SIM
  
)

allIndividualsList <- list(
  tar_combine(
    allIndividuals,
    use_names = FALSE,
    individualSimulations,
    # command = dplyr::bind_rows(!!!.x))
    command = list(!!!.x))
)

coreMultiverse <- list(
  tar_map(
    unlist = TRUE,
    values = values_Regime,
    tar_target(sampDuraFreqData,
               subset_duration(
                 allIndividualData = subset_frequency(
                   allIndividualData = allIndividuals,
                   freqPreset = tf),
                 daysDuration = td),
               priority = 0.92),
    tar_target(areaBasedAvailUse,
               area_based_extraction(
                 allIndividualData = sampDuraFreqData,
                 optionsList = optionsList_area
               ),
               priority = 0.9),
    tar_target(areaBasedOUT,
               area_based_calculations(
                 availUseData = areaBasedAvailUse,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_areaMethods
               ),
               priority = 0.9),
    tar_target(ssfOUT,
               wrapper_indi_ssf(
                 allIndividualData = sampDuraFreqData,
                 optionsList = optionsList_sff
               ),
               priority = 0.9),
    tar_target(poisOUT,
               method_pois_inla(
                 allIndividualData = sampDuraFreqData,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9),
    tar_target(twoStepOUT,
               method_twoStep(
                 allIndividualData = sampDuraFreqData,
                 sampleGroups = optionsList_samples,
                 optionsList = optionsList_pois),
               priority = 0.9)
  )
)

ssfCompiled <- list(
  tar_combine(
    ssfResults,
    coreMultiverse[[1]][grep("ssfOUT", names(coreMultiverse[[1]]))],
    command = list(!!!.x),
    priority = 0.8
  ),
  tar_target(
    ssfSampled,
    sample_ssf_results(
      ssfResults,
      sampleGroups = optionsList_samples,
      optionsList = optionsList_sff
    )
  ),
  tar_target(
    ssfSpecCurve,
    generate_spec_curves(
      outputResults = ssfSampled,
      method = "ssf"
    )
  ),
  tar_target(
    ssfBrms,
    run_brms(
      resultsData = ssfSampled
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
    poisSpecCurve,
    generate_spec_curves(
      outputResults = poisResults,
      method = "pois"
    )
  ),
  tar_target(
    poisBrms,
    run_brms(
      resultsData = poisResults
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
    twoStepSpecCurve,
    generate_spec_curves(
      outputResults = twoStepResults,
      method = "twoStep"
    )
  ),
  tar_target(
    twoStepBrms,
    run_brms(
      resultsData = twoStepResults
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
    areaSpecCurve,
    generate_spec_curves(
      outputResults = areaBasedResults,
      method = "area"
    )
  ),
  tar_target(
    areaBrms,
    run_brms(
      resultsData = areaBasedResults
    )
  )
)

brmModelOutputs <- list(
  tar_combine(
    modelsBrms,
    # manually pull out the brms model outputs
    list(ssfCompiled[[4]],
         areaBasedCompiled[[3]],
         poisCompiled[[3]],
         twoStepCompiled[[3]]),
    command = list(!!!.x),
    priority = 0.5
  ),
  tar_target(
    diagnosticPlots,
    diagnostics_brms(modelsList = modelsBrms)
  ),
  tar_target(
    modelExtracts,
    extract_model_values(modelsList = modelsBrms)
  ),
  tar_target(
    effectPlots,
    generate_effect_plots(modelsList = modelsBrms)
  ),
  tar_target(
    rmdRender,
    render_rmd(modelExtracts, effectPlots),
    cue = tar_cue(mode = "always"),
    priority = 0.1
  )
)


# All targets lists -------------------------------------------------------

list(individualSimulationsList,
     allIndividualsList,
     coreMultiverse,
     ssfCompiled,
     poisCompiled,
     twoStepCompiled,
     areaBasedCompiled,
     brmModelOutputs
)

# Examine -----------------------------------------------------------------

# targets::tar_make(starts_with("sampDuraFreq"))
# targets::tar_visnetwork(names = starts_with("combined"))
# mani <- targets::tar_manifest()
# pattern <- paste0("ssfOUT.*(",
#                   paste(sample(1:4, 2, replace = FALSE), collapse = "|"),
#                   ")_BADGER")
# mani$name[grep("ssfOUT.*4_BADGER", mani$name)]
# mani$name[grep(pattern, mani$name)]

