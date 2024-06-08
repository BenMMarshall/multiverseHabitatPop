#' Calculate population selection using a Poisson model
#'
#' @name method_pois_inla
#' @description A
#' @param allIndividualData All simulated movement data
#' @param sampleGroups list if IDs for the sample groups
#' @param optionsList options list that will be used to loop through options
#' @return Population estimates using a Poisson model.
#'
#' @export
method_pois_inla <- function(allIndividualData, sampleGroups, optionsList){
  # allIndividualData <- sampDuraFreqData_15_1
  # optionsList <- optionsList_pois
  # sampleGroups <- optionsList_samples
  landscape <- allIndividualData$landscape
  
  optionsForm <- optionsList$MethodPois_mf
  optionsLand <- optionsList$MethodPois_land
  optionsASteps <- optionsList$MethodPois_as
  optionsStepD <- optionsList$MethodPois_sd
  optionsTurnD <- optionsList$MethodPois_td
  
  listLength <- length(optionsForm) *
    length(optionsLand) *
    length(optionsASteps) *
    length(optionsStepD) *
    length(optionsTurnD) *
    length(names(sampleGroups))
  
  inlaOUTList <- vector("list", length = listLength)
  i <- 0
  for(sampID in names(sampleGroups)){
    print(sampID)
    # sampID <- "samp1"
    IDs <- optionsList_samples[[sampID]]
    speciesCodeLetter <- stringr::str_extract(names(allIndividualData)[1], ".$")
    IDs <- paste0("simData_",
                  speciesCodeLetter, # add in species code letter
                  "_i", sprintf("%03d", IDs))
    
    sampledIndividualData <- allIndividualData[names(allIndividualData) %in% IDs]
    
    movementData <- do.call(rbind, lapply(sampledIndividualData, function(x){
      x$locations
    }))
    
    print("--- movementData subsampled")
    
    # select only the information we need
    movementData <- movementData %>% 
      dplyr::mutate(datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S",
                                          tz = "UTC")) %>% 
      dplyr::select("x" = x, "y" = y,
                    "t" = datetime, "id" = id) %>% 
      dplyr::group_by(id) %>% 
      dplyr::arrange(t)
    
    # separate the individuals out into a list-column of dataframes, each item an animal
    movementDataNest <- movementData %>% 
      tidyr::nest(moveData = -id)
    
    print(sum(is.na(movementData$t)))
    
    # map operates a lot like an apply or loop. It repeats the function to each item
    # in the list. In this case we make all the individual dataframes into track
    # objects. Check dat_all to see the list of dataframes now has a second
    # dataframe/tibble for the track.
    movementDataNest <- movementDataNest %>% 
      mutate(trk = purrr::map(.x = moveData, .f = function(d){
        make_track(d, .x = x, .y = y, .t = t, crs = 32601)
      }))
    
    # Here the summarize_sampling_rate is repeated on each track object to give you
    # an individual level summary.
    # movementData %>%
    #   mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
    #   dplyr::select(id, sr) %>%
    #   unnest(cols = c(sr))
    
    for(as in optionsASteps){
      # as <- optionsASteps[1]
      for(sd in optionsStepD){
        # sd <- optionsStepD[1]
        for(td in optionsTurnD){
          # td <- optionsTurnD[1]
          for(land in optionsLand){
            # land <- optionsLand[1]
            # had to modify the code and avoid map to make sure the distributions are
            # based on a single individual, issues with the sl_ and ta_ being passed to
            # the fit_distr functions inside map
            allTracksList <- vector("list", length = length(movementDataNest$id))
            names(allTracksList) <- movementDataNest$id
            for(indiID in movementDataNest$id){
              print(indiID)
              # indiID <- "BADGER_i009"
              # which(movementDataNest$id == indiID)
              
              indiTrack <- movementDataNest$trk[[which(movementDataNest$id == indiID)]] %>% 
                steps() %>% 
                filter(sl_>0)
              
              indiTrackCov <- indiTrack %>% # removing the non-moves, or under GPS error
                random_steps(
                  n_control = as,
                  sl_distr = amt::fit_distr(x = indiTrack$sl_, dist_name = sd),
                  ta_distr = amt::fit_distr(x = indiTrack$ta_, dist_name = td)
                ) %>% 
                extract_covariates(landscape[[land]]) %>% 
                mutate(id = indiID)
              
              # print(unique(indiTrackCov$layer))
              # need a while loop to dodge the very very rare instances of NA from generated steps
              j <- 0
              while(any(is.na(unique(indiTrackCov$layer)))){
                j <- j+1
                indiTrackCov <- indiTrack %>% # removing the non-moves, or under GPS error
                  random_steps(
                    n_control = as,
                    sl_distr = amt::fit_distr(x = indiTrack$sl_, dist_name = sd),
                    ta_distr = amt::fit_distr(x = indiTrack$ta_, dist_name = td)
                  ) %>% 
                  extract_covariates(landscape[[land]]) %>% 
                  mutate(id = indiID)
                
                print(unique(indiTrackCov$layer))
                print(sum(is.na(indiTrackCov$layer)))
                if(j > 50){
                  print("Over 50 attempts for non NAs in extracted values, dropping NA instead when fewer than 10 instances")
                  
                  if(sum(is.na(indiTrackCov$layer)) < 10){
                    indiTrackCov <- indiTrackCov %>% 
                      filter(!is.na(layer))
                    {break}
                  }
                }
              }
              
              allTracksList[[indiID]] <- indiTrackCov
            }
            poisModelData <- do.call(rbind, allTracksList)
            
            print("--- poisModelData generated")
            
            poisModelData <- poisModelData %>% 
              mutate(
                y = as.numeric(case_),
                id = as.numeric(factor(id)), 
                step_id = paste0(id, step_id_, sep = "-"),
                cos_ta = cos(ta_), 
                log_sl = log(sl_),
                layer = factor(paste0("c", layer),
                               levels = c("c0", "c2")))
            
            # We can run the INLA model using the priors and set-up from Muff et al.
            # Precision for the priors of slope coefficients
            prec.beta.trls <- 1e-4
            
            for(form in optionsForm){
              # form <- optionsForm[1]
              if(form == "mf.is"){
                
                inlaFormula <- y ~ -1 + 
                  layer + # fixed covariate effect
                  layer:log_sl + # covar iteractions
                  layer:cos_ta + 
                  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
                  f(id, layer, values = 1:length(unique(poisModelData$id)), model="iid",
                    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                                              prior = "pc.prec", param = c(1, 0.05)))) 
                
              } else if(form == "mf.ss"){
                
                inlaFormula <- y ~ -1 + 
                  layer + # fixed covariate effect
                  f(step_id, model="iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
                  f(id, layer, values = 1:length(unique(poisModelData$id)), model="iid",
                    hyper = list(theta = list(initial = log(1), fixed = FALSE, 
                                              prior = "pc.prec", param = c(1, 0.05)))) 
              } # if else end
              
              # print(sum(poisModelData$y[poisModelData$layer == "c0"]))
              # print(sum(poisModelData$y[poisModelData$layer == "c2"]))
              # print(sum(is.na(poisModelData$y)))
              # print(sum(is.na(poisModelData$log_sl)))
              # print(sum(is.na(poisModelData$cos_ta)))
              # print(sum(is.na(poisModelData$id)))
              
              # natural model
              inlaOUT <- try(
                inla(inlaFormula,
                     family = "Poisson",
                     data = poisModelData, #verbose=TRUE,
                     control.fixed = list(
                       mean = 0,
                       prec = list(default = prec.beta.trls)),
                     control.inla = list(control.vb = list(emergency = 30)))
              )
              
              if(class(inlaOUT)[1] == "try-error"){
                
                inlaResults <- as.data.frame(matrix(NA, 2, 7))
                inlaResults$term <- c("layerc2", "layerc0")
                names(inlaResults) <- c("mean", "sd", "q025", "q50", "q975",
                                        "mode", "kld", "term")
                inlaResults$mmarginal <- NA
                inlaResults$emarginal <- NA
                
              } else {
                
                inlaResults <- inlaOUT$summary.fixed[1:2,]
                inlaResults$term <- row.names(inlaResults)
                names(inlaResults) <- c("mean", "sd", "q025", "q50", "q975",
                                        "mode", "kld", "term")
                inlaResults$mmarginal <- inla_mmarginal(inlaOUT)
                inlaResults$emarginal <- inla_emarginal(inlaOUT)
                
              }
              
              print(inlaResults$mean)
              
              optionsInfo <-
                data.frame(
                  sampleID = sampID,
                  species = speciesCodeLetter,
                  sampleSize = length(IDs),
                  trackFreq = allIndividualData[[2]]$trackFreq,
                  trackDura = allIndividualData[[2]]$trackDura,
                  analysis = "Poisson",
                  classLandscape = land,
                  modelFormula = form,
                  availablePerStep = as,
                  stepDist = sd,
                  turnDist = td
                )
              
              print(optionsInfo)
              print("INLA Results")
              
              combinedResults <- cbind(optionsInfo,
                                       inlaResults)
              i <- i + 1
              inlaOUTList[[i]] <- combinedResults
              
            } #form
          } #land
        }
      }
    }
  }
  return(do.call(rbind, inlaOUTList))
}
