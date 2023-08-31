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
  packages = c("tibble", "qs"), # packages that your targets need to run
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
  individual = paste0("i", 1:50)
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

optionsList_sff <- list(
  Method_method = c("ssf"),
  MethodSSF_mf = c("mf.is"),
  MethodSSF_sd = c("gamma"),
  MethodSSF_td = c("vonmises"),
  MethodSSF_as = 10
  # MethodSSF_mf = c("mf.is", "mf.ss"),
  # MethodSSF_sd = c("gamma", "exp"),
  # MethodSSF_td = c("vonmises", "unif"),
  # MethodSSF_as = as.integer(round(exp(seq(log(5), log(500), length.out = 5)), digits = 1))
)

# optionsList_pois <- list(
#   
# )

optionsList_ssfCombine <- list(
  weighted = c("weighted", "notweigthed")
)

# Targets workflow lists --------------------------------------------------

allIndividualEstimatesList <- list(
  ## LANDSCAPE SIMULATION
  regimeData <- tar_map(
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
        priority = 0.93), # FUNCTION simulate_individual
      
      ## DURATION + FREQUENCY MAP
      tar_map(
        unlist = TRUE,
        values = values_Regime,
        tar_target(sampDuraFreqData,
                   subset_duration(
                     movementData = subset_frequency(movementData = simData$locations,
                                                     freqPreset = tf),
                     daysDuration = td),
                   priority = 0.92),
        
        
        ## SSF
        tar_target(ssfOUT,
                   wrapper_indi_ssf(
                     movementData = sampDuraFreqData,
                     landscape = landscape,
                     optionsList = optionsList_sff
                   ),
                   priority = 0.9)
        # 
        ## POIS BRANCH POINT
        # tar_target(poisData,
        #            pop_sample(
        #              movementData = sampDuraFreqData,
        #              landscape = landscape,
        #              optionsList = sampleSize,
        #            ),
        #            priority = 0.9),
      ) # DURA FREQ
    ) # INDI SIM
  ) # LANDSCAPE SIM
)

values_Sample <-
  list(sampleSize = c(5,10,15,25,50))
# c(5,10,20,40)

set.seed(1)

sampleIDs <- lapply(values_Sample$sampleSize, function(x){
  sample(1:50, x, replace = FALSE)
})

######## WILL NEED AN EXTRA COLUMN IN THE DATA TO SEPERATE OUT REGIMES
sampledData_01 <- tar_combine(
  combinedMovementData_01,
  use_names = FALSE,
  regimeData[eval_select(matches(paste("sampDuraFreqData.*(",
                                       paste0(
                                         paste0("\\_i", sampleIDs[[1]]),
                                         collapse = "|"), ")",
                                       sep = "")), regimeData)],
  command = dplyr::bind_rows(!!!.x)
)

sampledData_02 <- tar_combine(
  combinedMovementData_02,
  use_names = FALSE,
  regimeData[eval_select(matches(paste("sampDuraFreqData.*(",
                                       paste0(
                                         paste0("\\_i", sampleIDs[[2]]),
                                         collapse = "|"), ")",
                                       sep = "")), regimeData)],
  command = dplyr::bind_rows(!!!.x)
)
# 
# sampledData_03 <- tar_combine(
#   combinedMovementData_03,
#   use_names = FALSE,
#   regimeData[eval_select(matches(paste("sampDuraFreqData.*(",
#                                        paste0(
#                                          paste0("\\_i", sampleIDs[[3]]),
#                                          collapse = "|"), ")",
#                                        sep = "")), regimeData)],
#   command = dplyr::bind_rows(!!!.x)
# )
# 
# sampledData_04 <- tar_combine(
#   combinedMovementData_04,
#   use_names = FALSE,
#   regimeData[eval_select(matches(paste("sampDuraFreqData.*(",
#                                        paste0(
#                                          paste0("\\_i", sampleIDs[[4]]),
#                                          collapse = "|"), ")",
#                                        sep = "")), regimeData)],
#   command = dplyr::bind_rows(!!!.x)
# )
# 
sampledData_05 <- tar_combine(
  combinedMovementData_05,
  use_names = FALSE,
  regimeData[eval_select(matches(paste("sampDuraFreqData.*(",
                                       paste0(
                                         paste0("\\_i", sampleIDs[[5]]),
                                         collapse = "|"), ")",
                                       sep = "")), regimeData)],
  command = dplyr::bind_rows(!!!.x)
)

# for(f in c(7,15)){
#   assign(paste0("sampledData_", f), 
#     tar_combine(
#       movementData,
#       allIndividualEstimatesList[[1]][grep(paste0("sampDuraFreqData_", f),
#                                            names(allIndividualEstimatesList[[1]]))],
#       command = rbind(!!!.x))
#     )
# }

# sampledData <- list(mget(ls(pattern = "sampledData_")))

# ssfCompiled <- tar_combine(
#   ssfResults,
#   allIndividualEstimatesList[[1]][grep("ssfOUT", names(allIndividualEstimatesList[[1]]))],
#   # command = list(!!!.x),
#   command = rbind(!!!.x),
#   priority = 0.8
# )


#### ASSINGMENT TO SAMPLES CAN OCCUR HERE WITHIN THE FUNCTION, SIMPLY SUBSET
#### FROM THE COMPILED DATAFRAME
# ssfPopulation <- list(
#   tar_map(
#     values = optionsList_ssfCombine,
#     tar_target(ssfPopulation,
#                combine_indi_estimates(
#                  indiSSF = ssfResults,
#                  # samples = c() # SAMPLES HERE
#                  method = "SSF"),
#                priority = 0.7)
#   )
# )

# pattern <- paste0("sampDuraFreqData_15_0.5_.*(",
#                   paste(sample(1:4, 2, replace = FALSE), collapse = "|"),
#                   ")")
# 
# values_Sample <- 
#   list(sampleSize = c(3,3,4))
# c(5,10,20,40)

# sample(1:4, 2, replace = FALSE)

# movementCombined <- list(
#   tar_map(
#     values = values_Sample,
#     tar_target(
#       name = movementDataSample,
#       command = example_function(x, sampleSize),
#       pattern = paste0("sampDuraFreqData_15_0.5_.*(",
#                        paste(sample(1:4, sampleSize, replace = FALSE), collapse = "|"),
#                        ")")
#     )
#       # tar_combine(
#       #   movementData,
#       #   allIndividualEstimatesList[[1]][grep(paste0("sampDuraFreqData_15_0.5_.*(",
#       #                                   paste(sample(1:4, 2, replace = FALSE),
#       #                                         collapse = "|"),
#       #                                   ")"),
#       #                            names(allIndividualEstimatesList[[1]]))],
#       #   # command = list(!!!.x),
#       #   command = rbind(!!!.x),
#       #   priority = 0.8
#       # )
#   )
# )

# ## Poisson Model
# tar_map(
#   values = optionsList_poisSamples,
#     tar_target(poisOUT,
#            wrapper_pop_pois(
#              movementData = sampDuraFreqData,
#              landscape = landscape,
#              optionsList = optionsList_pois
#            ),
#            priority = 0.9)
# )

# All targets lists -------------------------------------------------------

# sampledDataList <- list(mget(ls(pattern = "sampledData_")))

list(allIndividualEstimatesList,
     # sampledDataList
     sampledData_01,
     sampledData_02,
     # sampledData_03,
     # sampledData_04,
     sampledData_05
     # ssfCompiled,
     # ssfPopulation
)

# Examine -----------------------------------------------------------------

# targets::tar_visnetwork(names = starts_with("combined"))
# mani <- targets::tar_manifest()
# pattern <- paste0("ssfOUT.*(",
#                   paste(sample(1:4, 2, replace = FALSE), collapse = "|"),
#                   ")_BADGER")
# mani$name[grep("ssfOUT.*4_BADGER", mani$name)]
# mani$name[grep(pattern, mani$name)]

