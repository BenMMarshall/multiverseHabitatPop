#' Simulate individual
#'
#' @name simulate_individual
#' @description A
#' @param individualNum number for the repeats of a given sim style
#' @param species c("BADGER", "VULTURE", "KINGCOBRA"), add other sim styles here
#' @param simSteps = 24\*60 \*365
#' @param desOptions = 10
#' @param options = 12
#' @param landscapeList
#' @param seed
#' @return a
#'
#' @export
simulate_individual <- function(
    individualNum,
    species,
    simSteps = 24*60 *365,
    desOptions = 10,
    options = 12,
    landscapeList,
    seed
){
  # suppressWarning about raster version
  suppressWarnings({
    pcks <- c("abmAnimalMovement", "raster")
    missing <- !unlist(lapply(
      pcks,
      require,
      character.only = TRUE))
  })

  if(any(missing)){
    stop(paste0(pcks[missing], " not installed"))
  }

  # # example of where branch sim approach could be added, maybe or better approach would be new tree
  # if(species == "RSF_SIM"){
  #   # atm RSF simulation approach here, but still needs a SPECIES connection to get the same landscape
  # }

  ## SPECIES INPUTS

  b0 <- c(0.97, 0.01, 0.001) # rest
  b1 <- c(0.0002, 0.95, 0.0008) # explore/move
  b2 <- c(0.001, 0.00001, 0.99) # forage

  Default_behaveMatrix <- rbind(b0, b1, b2)

  colnames(Default_behaveMatrix) <- c("b0", "b1", "b2")

  Default_behaveMatrix


  ## BADGER SERVES AS DEFAULT
  # sampledShelters <- raster::sampleRandom(raster::raster(landscapeList$shelter), 2,
  #                                         ext = raster::extent(0.45, 0.65, 0.45, 0.65),
  #                                         rowcol = TRUE)
  extractedMovementValues <- c(0, 0)
  while(!all(extractedMovementValues == 2)){
    sampledShelters <- raster::sampleRandom(raster::raster(landscapeList$shelter), 2,
                                            ext = raster::extent(0.45, 0.65, 0.45, 0.65),
                                            rowcol = TRUE)

    spMovementPoints <- sp::SpatialPoints(sampledShelters[,c("col", "row")],
                                          sp::CRS(SRS_string = "EPSG:32601"))
    extractedMovementValues <- raster::extract(landscapeList$classRaster, spMovementPoints)
  }

  IN_shelterLocs <- data.frame(
    "x" = sampledShelters[,2],
    "y" = sampledShelters[,1])
  IN_shelterSize <- 8
  IN_k_step <- c(0.3*60, 1.25*60, 0.25*60)
  IN_s_step <- c(0.8, 0.25, 0.5)
  IN_mu_angle <- c(0, 0, 0)
  IN_k_angle <- c(0.6, 0.99, 0.6)
  IN_destinationRange <- c(3, 120)
  IN_destinationDirection <- c(0, 0.01)
  IN_destinationTransformation <- 2
  IN_destinationModifier <- 2
  IN_rescale <- 5
  sampledAvoid <- raster::sampleRandom(raster::raster(landscapeList$forage), 3,
                                       ext = raster::extent(0.45, 0.65, 0.45, 0.65),
                                       rowcol = TRUE)
  IN_avoidLocs <- data.frame(
    "x" = sampledAvoid[,2],
    "y" = sampledAvoid[,1])
  IN_avoidTransformation <- 2
  IN_avoidModifier <- 4
  IN_rest_Cycle <- c(0.12, 0, 24, 24)
  c0 <- c(0.075, 0, 24* (365/2), 24* 365) # seasonal
  IN_additional_Cycles <- rbind(c0)
  IN_behaveMatrix <- Default_behaveMatrix

  ## VULTURE CHANGES
  if(species == "VULTURE"){

    # sampledShelters <- raster::sampleRandom(raster::raster(landscapeList$shelter), 3,
    #                                         ext = raster::extent(0.45, 0.65, 0.45, 0.65),
    #                                         rowcol = TRUE)
    extractedMovementValues <- c(0, 0)
    while(!all(extractedMovementValues == 2)){
      sampledShelters <- raster::sampleRandom(raster::raster(landscapeList$shelter), 2,
                                              ext = raster::extent(0.45, 0.65, 0.45, 0.65),
                                              rowcol = TRUE)

      spMovementPoints <- sp::SpatialPoints(sampledShelters[,c("col", "row")],
                                            sp::CRS(SRS_string = "EPSG:32601"))
      extractedMovementValues <- raster::extract(landscapeList$classRaster, spMovementPoints)
    }

    IN_shelterLocs <- data.frame(
      "x" = sampledShelters[,2],
      "y" = sampledShelters[,1])
    IN_shelterSize <- 5
    IN_k_step <- c(2, 2.2*60, 1.5*60)
    IN_s_step <- c(40, 1.2, 1)
    IN_mu_angle <- c(0, 0, 0)
    IN_k_angle <- c(0.6, 0.99, 0.6)
    IN_destinationRange <- c(50, 120)
    IN_destinationDirection <- c(0, 0.01)
    IN_destinationTransformation <- 2
    IN_destinationModifier <- 2
    IN_rescale <- 20
    IN_rest_Cycle <- c(0.1, 0, 24, 24)
    c0 <- c(0.025, 0, 24* (365/2), 24* 365) # seasonal
    IN_additional_Cycles <- rbind(c0)
    IN_behaveMatrix <- Default_behaveMatrix
    IN_behaveMatrix[2,3] <- 0.0002
    IN_behaveMatrix[3,2] <- 0.000015
  }

  ## KINGCOBRA CHANGES
  if(species == "KINGCOBRA"){
    IN_k_step <- c(30, 40, 20)
    IN_s_step <- c(0.75, 1.2, 1.75)
    IN_mu_angle <- c(0, 0, 0)
    IN_k_angle <- c(0.6, 0.99, 0.6)
    IN_destinationRange <- c(50, 10)
    IN_destinationDirection <- c(0, 0.01)
    IN_destinationTransformation <- 2
    IN_destinationModifier <- 1.5
    IN_rescale <- 10
    IN_rest_Cycle <- c(0.14, 0, 24, 24)
    c0 <- c(0.12, 0, 24, 24*4) # digestion
    c1 <- c(0.05, 0, 24 * (365/2), 24* 365 ) # seasonal
    IN_additional_Cycles <- rbind(c0, c1)
    IN_behaveMatrix <- Default_behaveMatrix
    IN_behaveMatrix[1,1] <- 0.95
    IN_behaveMatrix[1,2] <- 0.005
    IN_behaveMatrix[3,1] <- 0.00025
    IN_behaveMatrix[3,2] <- 0.000001
    IN_behaveMatrix[3,3] <- 0.999

    # sampledShelters <- raster::sampleRandom(raster::raster(landscapeList$shelter), 12,
    #                                         ext = raster::extent(0.35, 0.65, 0.42, 0.65),
    #                                         rowcol = TRUE)
    extractedMovementValues <- c(0, 0)
    while(!all(extractedMovementValues == 2)){
      sampledShelters <- raster::sampleRandom(raster::raster(landscapeList$shelter), 2,
                                              ext = raster::extent(0.45, 0.65, 0.45, 0.65),
                                              rowcol = TRUE)

      spMovementPoints <- sp::SpatialPoints(sampledShelters[,c("col", "row")],
                                            sp::CRS(SRS_string = "EPSG:32601"))
      extractedMovementValues <- raster::extract(landscapeList$classRaster, spMovementPoints)
    }

    IN_shelterLocs <- data.frame(
      "x" = sampledShelters[,2],
      "y" = sampledShelters[,1])
    IN_shelterSize <- 10
    sampledAvoid <- raster::sampleRandom(raster::raster(landscapeList$forage), 3,
                                         ext = raster::extent(0.45, 0.65, 0.45, 0.65),
                                         rowcol = TRUE)
    IN_avoidLocs <- data.frame(
      "x" = sampledAvoid[,2],
      "y" = sampledAvoid[,1])
    # and specifying a weak avoidance of the points
    IN_avoidTransformation <- 2
    IN_avoidModifier <- 1
  }

  ## MAIN FUNCTION
  startLocation <- sample(800:1200, 2, replace = TRUE)

  simResults <- abm_simulate(
    start = startLocation,
    timesteps = simSteps,
    des_options = desOptions,
    options = options,
    k_step = IN_k_step,
    s_step = IN_s_step,
    mu_angle = IN_mu_angle,
    k_angle = IN_k_angle,
    rescale_step2cell = IN_rescale,
    shelterLocations = IN_shelterLocs,
    shelterSize = IN_shelterSize,
    avoidPoints = IN_avoidLocs,
    destinationRange = IN_destinationRange,
    destinationDirection = IN_destinationDirection,
    destinationTransformation = IN_destinationTransformation,
    destinationModifier = IN_destinationModifier,
    avoidTransformation = IN_avoidTransformation,
    avoidModifier = IN_avoidModifier,
    behave_Tmat = IN_behaveMatrix,
    rest_Cycle = IN_rest_Cycle,
    additional_Cycles = IN_additional_Cycles,
    shelteringMatrix = landscapeList$shelter,
    foragingMatrix = landscapeList$forage,
    movementMatrix = landscapeList$movement
  )

  simResults$locations$datetime <- as.POSIXct(simResults$locations$timestep * 60,
                                              origin = "2022-01-02 00:00:00")

  simResults$locations$hour <- as.numeric(substr(simResults$locations$datetime, 12, 13))
  simResults$locations$minute <- as.numeric(substr(simResults$locations$datetime, 15, 16))
  simResults$locations$yday <- as.numeric(format(simResults$locations$datetime,"%j"))

  simResults$locations$id <- paste0(species, "_", individualNum)

  return(simResults)

}
