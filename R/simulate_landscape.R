#' Simulate landscapes
#'
#' @name simulate_landscape
#' @description A
#' @param species
#' @param seed
#' @return a
#'
#' @export
simulate_landscape <- function(
    species,
    seed
){

  pcks <- c("raster", "RandomFields", "NLMR")
  missing <- !unlist(lapply(
    pcks,
    require,
    character.only = TRUE))

  if(any(missing)){
    stop(paste0(pcks[missing], " not installed"))
  }

  seed <- 2022
  row = 2000; col = 2000
  RandomFields::RFoptions(install="no")
  baseLandscape <- suppressMessages(NLMR::nlm_gaussianfield(ncol = col,
                                                  nrow = row,
                                                  resolution = 1,
                                                  autocorr_range = 40,
                                                  mag_var = 5,
                                                  nug = 0.2,
                                                  mean = 0.5,
                                                  user_seed = seed,
                                                  rescale = TRUE))

  # duplicate the base landscape to be modified for each of the simulating layers
  forageQual <- baseLandscape
  moveQual <- baseLandscape
  shelterQual <- baseLandscape

  forageQual[forageQual[] < 0.4] <- 0
  # set min 0 max 1, normalise the values between 1 and 0
  forageQual[] <- (forageQual[] - min(forageQual[], na.rm = TRUE)) /
    (max(forageQual[], na.rm = TRUE) - min(forageQual[], na.rm = TRUE))

  # areas with high resources are accessible (> 0.6, increased by 0.5), but the
  # fastest least resistance routes are actually edge habitat areas (0.6 to 0.3,
  # increased by 1). Core areas of low resrouce are also difficult to move
  # through.
  moveQual[moveQual[] > 0.6] <- moveQual[moveQual[] > 0.6] + 0.5
  moveQual[moveQual[] < 0.6 & moveQual[] > 0.3] <-
    moveQual[moveQual[] < 0.6 & moveQual[] > 0.3] + 1
  moveQual[] <- (moveQual[] - min(moveQual[], na.rm = TRUE)) /
    (max(moveQual[], na.rm = TRUE) - min(moveQual[], na.rm = TRUE))


  # shelter sites are best found near the edge of high resource areas, but
  # deeper than the best movement routes
  shelterQual[shelterQual[] < 0.7 & shelterQual[] > 0.5] <-
    shelterQual[shelterQual[] < 0.7 & shelterQual[] > 0.5] + 1
  shelterQual[shelterQual[] < 0.5] <- 0
  shelterQual[] <- (shelterQual[] - min(shelterQual[], na.rm = TRUE)) /
    (max(shelterQual[], na.rm = TRUE) - min(shelterQual[], na.rm = TRUE))

  landscapeLayersList <- list(
    "shelter" = matrix(data = raster::getValues(shelterQual),
                       nrow = row,
                       ncol = col),
    "forage" = matrix(data = raster::getValues(forageQual),
                      nrow = row,
                      ncol = col),
    "movement" = matrix(data = raster::getValues(moveQual),
                        nrow = row,
                        ncol = col))
  # landscapeLayersList <- list(
  #   "shelter" = raster::as.matrix(shelterQual),
  #   "forage" = raster::as.matrix(forageQual),
  #   "movement" = raster::as.matrix(moveQual))


  if(species == "VULTURE"){
    # maximising movement ease by changing the default movement matrix
    VULTURE_movementMatrix <- landscapeLayersList$movement
    VULTURE_movementMatrix[] <- 1
    landscapeLayersList$movement <- VULTURE_movementMatrix
    # reducing foraging quality in the West of the landscape
    VULTURE_forageMatrix <- landscapeLayersList$forage
    VULTURE_forageMatrix[1:950,1:2000] <- VULTURE_forageMatrix[1:950,1:2000] - 0.6
    VULTURE_forageMatrix[VULTURE_forageMatrix[] < 0] <- 0
    landscapeLayersList$forage <- VULTURE_forageMatrix

  }

  if(species == "KINGCOBRA"){

    # extracting the default matrices ready for changes
    KINGCOBRA_shelteringMatrix <- landscapeLayersList$shelter
    KINGCOBRA_forageMatrix <- landscapeLayersList$forage
    KINGCOBRA_movementMatrix <- landscapeLayersList$movement
    # defining the start and end of strong intersections hampering movement
    roadMin_x <- 1360
    roadMax_x <- roadMin_x + 40
    roadMin_y <- 660
    roadMax_y <- roadMin_y + 40

    # applying the change, dramatically reducing the weighting for all matrices
    KINGCOBRA_shelteringMatrix[roadMin_x:roadMax_x,1:2000] <-
      KINGCOBRA_shelteringMatrix[roadMin_x:roadMax_x,1:2000] - 90
    KINGCOBRA_shelteringMatrix[1:2000,roadMin_y:roadMax_y] <-
      KINGCOBRA_shelteringMatrix[1:2000,roadMin_y:roadMax_y] - 90
    KINGCOBRA_shelteringMatrix[!KINGCOBRA_shelteringMatrix >= -99.9] <- -99

    KINGCOBRA_forageMatrix[roadMin_x:roadMax_x,1:2000] <-
      KINGCOBRA_forageMatrix[roadMin_x:roadMax_x,1:2000] - 90
    KINGCOBRA_forageMatrix[1:2000,roadMin_y:roadMax_y] <-
      KINGCOBRA_forageMatrix[1:2000,roadMin_y:roadMax_y] - 90
    KINGCOBRA_forageMatrix[!KINGCOBRA_forageMatrix >= -99.9] <- -99

    KINGCOBRA_movementMatrix[roadMin_x:roadMax_x,1:2000] <-
      KINGCOBRA_movementMatrix[roadMin_x:roadMax_x,1:2000] - 90
    KINGCOBRA_movementMatrix[1:2000,roadMin_y:roadMax_y] <-
      KINGCOBRA_movementMatrix[1:2000,roadMin_y:roadMax_y] - 90
    KINGCOBRA_movementMatrix[!KINGCOBRA_movementMatrix >= -99.9] <- -99

    # swap back into the landlayers list for export
    landscapeLayersList$shelter <- KINGCOBRA_shelteringMatrix
    landscapeLayersList$forage <- KINGCOBRA_forageMatrix
    landscapeLayersList$movement <- KINGCOBRA_movementMatrix

  }

  ## output should be classified too, to save compute later
  classLandscape <- baseLandscape
  # > 0.6
  # <= 0.6 & > 0.3
  # <= 0.3
  classLandscape[classLandscape[] > 0.5] <- 2
  classLandscape[classLandscape[] <= 0.5 & classLandscape[] > 0.3] <- 1
  classLandscape[classLandscape[] <= 0.3] <- 0


  classLandscapeList <- list(
    "classified" = matrix(data = raster::getValues(classLandscape),
                          nrow = row,
                          ncol = col))

  landscapeLayersList$classified <- classLandscapeList$classified

  classRaster <- raster::raster(nrows = nrow(landscapeLayersList$classified),
                                ncols = ncol(landscapeLayersList$classified),
                                xmn = 0, xmx = nrow(landscapeLayersList$classified),
                                ymn = 0, ymx = ncol(landscapeLayersList$classified),
                                crs = CRS(SRS_string = "EPSG:32601"),
                                # need to transpose cos matrix and raster deal with rows and col differently
                                vals = t(landscapeLayersList$classified))
  # and flip to full match the raster with the matrix used in the sims
  classRaster <- raster::flip(classRaster)


  classRasterList <- list(
    "classRaster" = raster::flip(classRaster))

  landscapeLayersList$classRaster <- classRasterList$classRaster

  rBase <- landscapeLayersList$classRaster
  rAll <- raster::projectRaster(from = rBase,
                                to = projectExtent(rBase, crs = sp::CRS(SRS_string = "EPSG:4326")))
  rAll[] <- as.factor(paste0("c", round(rAll[], digits = 0)))

  classRasterLatLonList <- list(
    "classRasterLatLon" = rAll)

  landscapeLayersList$classRasterLatLon <- classRasterLatLonList$classRasterLatLon

  return(landscapeLayersList)

}
