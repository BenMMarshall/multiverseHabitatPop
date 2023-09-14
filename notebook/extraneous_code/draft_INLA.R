
# Code orignally from:
#' @article{doi:10.1111/1365-2656.13087,
#'   author = {Muff, Stefanie and Signer, Johannes and Fieberg, John},
#'   title = {Accounting for individual-specific variation in habitat-selection studies: Efficient estimation of mixed-effects models using Bayesian or frequentist computation},
#'   journal = {Journal of Animal Ecology},
#'   volume = {89},
#'   number = {1},
#'   pages = {80-92},
#'   doi = {10.1111/1365-2656.13087},
#'   url = {https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/1365-2656.13087},
#'   eprint = {https://besjournals.onlinelibrary.wiley.com/doi/pdf/10.1111/1365-2656.13087},
#'   year = {2020}
#' }
#' 

# Edited by Benjamin Michael Marshall and Samantha Nicole Smith 2020-09-02 for publication of: 
# title = {Restricted constrictors: Space use and habitat selection of Burmese pythons in Northeast Thailand}
## authors= {Samantha Nicole Smith, Max Dolton Jones, Benjamin Michael Marshall, Surachit Waengsothorn, George A. Gale and Colin Thomas Strine}
## year = {2020}

library(dplyr)
library(raster)
library(amt)
library(ggplot2)

# install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#install.packages("INLA", lib = .libPaths()[3], 
#                 repos=c(getOption("repos"),
#                 INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)

targets::tar_load("sampDuraFreqData_15_1")

allIndividualData <- sampDuraFreqData_15_1

landscape <- allIndividualData$landscape

allIndividualData[-1]

movementData <- do.call(rbind, lapply(allIndividualData[-1], function(x){
  x$locations
}))

# Read in animal tracking data --------------------------------------------

optionsForm <- optionsList$MethodPois_mf
optionsAPoints <- optionsList$MethodPois_as

listLength <- length(optionsForm)*
  length(optionsAPoints)*
  length(names(sampleGroups))

companaResultsList <- vector("list", length = listLength)
i <- 0

for(sampID in names(sampleGroups)){
  
  IDs <- optionsList_samples[[sampID]]
  IDs <- paste0("simData_i", sprintf("%03d", IDs))
  
  # select only the information we need
  dat <- movementData %>% 
    mutate(datetime = as.POSIXct(datetime, format = "%y-%m-%d %H:%M:%S",
                                 tz = "UTC")) %>% 
    dplyr::select("x" = x, "y" = y,
                  "t" = datetime, "id" = id) %>% 
    group_by(id) %>% 
    arrange(t) %>% 
    filter(id %in% sampleIDs)
  
  # separate the individuals out into a list-column of dataframes, each item an animal
  dat_all <- dat %>% 
    nest(data = -id) 
  
  # map operates a lot like an apply or loop. It repeats the function to each item
  # in the list. In this case we make all the individual dataframes into track
  # objects. Check dat_all to see the list of dataframes now has a second
  # dataframe/tibble for the track.
  dat_all <- dat_all %>% 
    mutate(trk = map(data, function(d){
      make_track(d, .x = x, .y = y, .t = t, crs = 32601)
    }))
  
  # Here the summarize_sampling_rate is repeated on each track object to give you
  # an individual level summary.
  # dat_all %>% 
  #   mutate(sr = lapply(trk, summarize_sampling_rate)) %>% 
  #   dplyr::select(id, sr) %>%
  #   unnest(cols = c(sr))
  
  # had to modify the code and avoid map to make sure the distributions are based
  # on a single individual, issues with the sl_ and ta_ being passed to the
  # fit_distr functions inside map
  allTracksList <- vector("list", length = length(dat_all$id))
  for(indiID in dat_all$id){
    # indiID <- "BADGER_i001"
    which(dat_all$id == indiID)
    
    indiTrack <- dat_all$trk[[which(dat_all$id == indiID)]] %>% 
      steps() %>% 
      filter(sl_>0)
    
    indiTrackCov <- indiTrack  %>% # removing the non-moves, or under GPS error
      random_steps(
        n_control = 10,
        sl_distr = amt::fit_distr(x = indiTrack$sl_, dist_name = "gamma"),
        ta_distr = amt::fit_distr(x = indiTrack$ta_, dist_name = "vonmises")
      ) %>% 
      extract_covariates(landscape$classRaster) %>% 
      mutate(id = indiID)
    
    allTracksList[[indiID]] <- indiTrackCov
  }
  poisModelData <- do.call(rbind, allTracksList)
  
  poisModelData <- poisModelData %>% 
    mutate(
      y = as.numeric(case_),
      id = as.numeric(factor(id)), 
      step_id = paste0(id, step_id_, sep = "-"),
      cos_ta = cos(ta_), 
      log_sl = log(sl_),
      layer = factor(paste0("c", layer)))
  
  
  
  
  
  
  
}


# if you want to add meta data, like sex, do so here. We had a second dataframe
# that included metadata

#dat_all <- left_join(dat_all,
#                     snake.meta %>% 
#                       select(id = deployment.id, sex = animal.sex))


# So dat_all is set up ready for the generation of random steps and the
# extraction of raster covariate values

# Step-Selection function -------------------------------------------------

# To make the random steps we use the map function again to do this on an
# individual basis. We calculate the step lengths so we have data to draw from
# to generate the random locations. Then we drop out the zero step lengths that
# break things, or we can drop the sub-GPS error points because we cannot be
# confident they were moves (ie not a new selection choice from the animal).
# Then the random steps, can be generous here because we aren't using GPS data
# with the high resolution that often entails. We have use more randoms to pick
# up smaller changes and maximise the high res rasters we have. Final step is to
# get the covariate values. After the mutate+map bit is a few lines to compile
# the data into a single v large dataframe for the model.


# dat_ssf <- dat_all %>% 
#   mutate(stps = map(trk, ~ .x %>%
#                       steps() %>% 
#                       filter(sl_>0) %>% # removing the non-moves, or under GPS error
#                       random_steps(
#                         n_control = 10
#                       ) %>% 
#                       extract_covariates(landscape$classRaster))) %>% 
#   dplyr::select(id, stps) %>%
#   unnest(cols = c(stps)) %>% 
#   group_by(id) %>% 
#   arrange(t1_) %>% 
#   ungroup() %>% 
#   mutate(
#     y = as.numeric(case_),
#     id = as.numeric(factor(id)), 
#     step_id = paste0(id, step_id_, sep = "-"),
#     cos_ta = cos(ta_), 
#     log_sl = log(sl_))

# dat_ssf$layer <- paste0("c", dat_ssf$layer)
# dat_ssf$layer <- factor(dat_ssf$layer)

# So dat_ssf is the correct format for the modelling. We have the id, a heap of
# movement columns, and the critical case_ that describes whether they used a
# location or not.

# Running the INLA model --------------------------------------------------

# We can run the INLA model using the priors and set-up from Muff et al.
# Precision for the priors of slope coefficients
prec.beta.trls <- 1e-4

# "In the model formula for INLA, we set the stratum-specific intercept variance
# to $10^6$ (or rather: the precision to $10^{ ^'6}$) by fixing it (`fixed=T`)
# to an initial value. The other precision is given a PC(1,0.05) prior:"

###### For BUCA data
### making formulas for each habitat feature

#Natural
formula.random.n <- y ~  -1 + 
  layer + # fixed covariate effect
  layer:log_sl + # covar iteractions
  layer:cos_ta + 
  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  f(id, layer, values = 1:length(unique(dat_ssf$id)), model="iid",
    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                              prior = "pc.prec", param = c(1, 0.05)))) 

# Last two parts of the model are the gaussian processes to deal with step_id
# (each move), and the id of the animal. So steps_id is the covariate, then in
# the second one is id that is weighted by the raster cov.

# Fit the models -----------------------------------------------------------

# natural model
inla.ssf.n <- inla(formula.random.n,
                   family = "Poisson",
                   data = poisModelData, #verbose=TRUE,
                   control.fixed = list(
                     mean = 0,
                     prec = list(default = prec.beta.trls)))

# Model results -----------------------------------------------------------

# We can see a dataframe of the fixed effects here, so each covar and any interaction terms
inla.ssf.n$summary.fixed

# "Since variances are parameterized and treated as precisions, the summary of
# the respective posterior distributions is given for the precisions:"
inla.ssf.n$summary.hyperpar


inla_mmarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.mmarginal(inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mode of variance"))
  results
}
inla_emarginal <- function(r.out){ 
  results <- sapply(r.out$marginals.hyperpar, 
                    function(y) 
                      inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mean of variance"))
  results
}

# "Posterior mean and mode are obtained as"
inla_emarginal(inla.ssf.n)
inla_mmarginal(inla.ssf.n)

fixed.df.n <- inla.ssf.n$summary.fixed
fixed.df.n$term <- row.names(fixed.df.n)
names(fixed.df.n) <- c("mean", "sd", "q025", "q50", "q975",
                       "mode", "kld", "term")

