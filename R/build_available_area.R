#' Generate the prerequisite for building a polygon for use as an available area
#'
#' @name build_available_area
#' @description A
#' @param movementData must have a x and y column for locations, and a datetime
#'   column for timestamps ("%Y-%m-%d %H:%M:%S")
#' @param method "MCP", "KDE_LSCV", "KDE_href", "AKDE", "dBBMM"
#' @param SRS_string = "EPSG:32601"
#' @param dBBMMsettings time duration of window and margin in hours, function
#'   translate duration to datapoints
#' @return An item that is the fit or UD that the a polygon can be
#'   extracted from using build_available_polygon
#'
#' @export
build_available_area <- function(movementData,
                                 method = c("MCP", "KDElscv", "KDEhref", "AKDE", "dBBMM"),
                                 SRS_string = "EPSG:32601",
                                 dBBMMsettings = NULL){


  if(method == "MCP"){

    suppressWarnings({
      area_OUT <- sp::SpatialPoints(movementData[,c("x", "y")], sp::CRS(SRS_string = "EPSG:32601"))
    })

  } else if(method == "KDElscv"){

    suppressWarnings({
      spPoints <- sp::SpatialPoints(movementData[,c("x", "y")], sp::CRS(SRS_string = "EPSG:32601"))
    })

    # area_OUT <- vector("list", 2)

    suppressWarnings({
      kdeLSCV_UD <- adehabitatHR::kernelUD(spPoints,
                                           h = "LSCV",
                                           grid = 240,
                                           same4all = FALSE,
                                           hlim = c(0.001, 2000), # might need to play with the limits to help convergence
                                           kern = "bivnorm",
                                           extent = 4,
                                           boundary = NULL)
    })
    area_OUT <- kdeLSCV_UD

  } else if(method == "KDEhref"){

    suppressWarnings({
      spPoints <- sp::SpatialPoints(movementData[,c("x", "y")], sp::CRS(SRS_string = "EPSG:32601"))
    })
    # area_OUT <- vector("list", 2)
    suppressWarnings({
      kdehref_UD <- adehabitatHR::kernelUD(spPoints,
                                           h = "href",
                                           grid = 240, # needs to be large enough to be smooth-ish
                                           same4all = FALSE,
                                           hlim = c(0.1, 1.5),
                                           kern = "bivnorm",
                                           extent = 4,
                                           boundary = NULL)
    })
    area_OUT <- kdehref_UD

  } else if(method == "AKDE"){

    suppressWarnings({
      spPoints <- sp::SpatialPoints(movementData[,c("x", "y")], sp::CRS(SRS_string = "EPSG:32601"))
    })

    suppressWarnings({
      spLL <- sp::spTransform(spPoints, sp::CRS(SRS_string = "EPSG:4326"))
    })

    movementData$lon <- spLL@coords[,1]
    movementData$lat <- spLL@coords[,2]

    # area_OUT <- vector("list", 2)

    teleObj <- ctmm::as.telemetry(movementData,
                                  timeformat = "%Y-%m-%d %H:%M:%S",
                                  timezone="UTC",
                                  projection = sp::CRS(SRS_string = "EPSG:32601"))

    print(attributes(teleObj)$info$identity)
    print("teleObj")
    # can do the slower one for the real deal
    # varioDataVar <- variogram(teleObj, fast = FALSE, CI = "Gauss")
    varioDataVar <- ctmm::variogram(teleObj, fast = TRUE)
    print("vario")
    guess <- ctmm::ctmm.guess(teleObj, interactive = FALSE)
    print("guess")
    # need to specify more cores???
    fits <- try(
      ctmm::ctmm.select(teleObj, guess, verbose = FALSE,
                        cores = 4, method = "pHREML")
    )
    print("fit")

    if(class(fits) == "try-error"){
      area_OUT <- fits
    } else {
      # akdeRes <- ctmm::akde(teleObj, fits[[1]],
      #                       weights = TRUE)
      # needed to catch weird instances with limited data
      akdeRes <- try(
        ctmm::akde(teleObj, fits,
                   weights = TRUE)
      )
      print("area")
      area_OUT <- akdeRes
    }

  } else if(method == "dBBMM"){

    if(is.null(dBBMMsettings)){
      stop("dBBMMsettings required for dBBMM running, 2 length vector of time duration of ws and mrg (in hours)")
    }

    # As our most infrequent tracking is 168 hours (1 week), we will set the
    # window to the number of data points collected over 168 hours, and a margin
    # of 48 hours. In cases where that is not possible set window and margin to minimum
    windowSize <- nrow(movementData[movementData$timestep <= dBBMMsettings[1]*60,])
    if(windowSize %% 2 == 0){
      windowSize <- windowSize - 1
    }
    if(windowSize < 7){
      windowSize <- 7
    }

    marginSize <- nrow(movementData[movementData$timestep <= dBBMMsettings[2]*60,])
    if(marginSize %% 2 == 0){
      marginSize <- marginSize - 1
    }
    if(marginSize < 3){
      marginSize <- 3
    }

    # area_OUT <- vector("list", 2)

    suppressWarnings({
      moveObj <- move::move(x = movementData$x, y = movementData$y,
                            time = movementData$datetime,
                            proj = sp::CRS(SRS_string = "EPSG:32601"))
    })

    set_grid.ext <- 4
    set_dimsize <- 640
    suppressWarnings({
      dbbmm <- move::brownian.bridge.dyn(object = moveObj,
                                         location.error = 0.1,
                                         window.size = windowSize,
                                         margin = marginSize,
                                         ext = set_grid.ext,
                                         dimSize = set_dimsize,
                                         verbose = FALSE)
    })
    area_OUT <- dbbmm

  }
  return(area_OUT)

}

# library(here)
#
# movementData <- read.csv(here("notebook", "prereg", "prelimData.csv"))
# movementData <- movementData[seq(0, 100000, by = 200),]
# movementData$datetime <- as.POSIXct(movementData$timestep * 60,
#            origin = "2022-01-01")



# # MCP ---------------------------------------------------------------------
#
# library(adehabitatHR)
#
# spPoints <- SpatialPoints(movementData[,c("x", "y")], CRS(SRS_string = "EPSG:32601"))
#
# mcp_OUT <- mcp(spPoints, percent = 100, unin = "m",
#                unout = "m2")
#
# plot(mcp_OUT)
# points(spPoints)
#
# # KDE LSCV ----------------------------------------------------------------
#
# spPoints <- SpatialPoints(movementData[,c("x", "y")], CRS(SRS_string = "EPSG:32601"))
#
# kdeLSCV_UD <- kernelUD(spPoints,
#                        h = "LSCV",
#                        grid = 120,
#                        same4all = FALSE,
#                        hlim = c(0.001, 2000), # might need to play with the limits to help convergence
#                        kern = "bivnorm",
#                        extent = 2,
#                        boundary = NULL)
#
# kdeLSCV_OUT <- getverticeshr(kdeLSCV_UD, 95)
#
# plot(kdeLSCV_OUT, add = TRUE)
#
# # KDE href ----------------------------------------------------------------
#
# kdehref_UD <- kernelUD(spPoints,
#                        h = "href",
#                        grid = 120, # needs to be large enough to be smooth-ish
#                        same4all = FALSE,
#                        hlim = c(0.1, 1.5),
#                        kern = "bivnorm",
#                        extent = 2,
#                        boundary = NULL)
#
# kdehref_OUT <- getverticeshr(kdehref_UD, 95)
#
# plot(kdehref_OUT, add = TRUE)
#
# # AKDE --------------------------------------------------------------------
#
# library(sp)
#
# spPoints <- SpatialPoints(movementData[,c("x", "y")], CRS(SRS_string = "EPSG:32601"))
# spLL <- spTransform(spPoints, CRS(SRS_string = "EPSG:4326"))
# movementData$lon <- spLL@coords[,1]
# movementData$lat <- spLL@coords[,2]
#
# library(ctmm)
#
# teleObj <- as.telemetry(movementData,
#                         timeformat = "%Y-%m-%d %H:%M:%S",
#                         timezone="UTC",
#                         projection = CRS(SRS_string = "EPSG:32601"))
#
# # can do the slower one for the real deal
# # varioDataVar <- variogram(teleObj, fast = FALSE, CI = "Gauss")
# varioDataVar <- variogram(teleObj, fast = TRUE)
# guess <- ctmm.guess(teleObj, interactive = FALSE)
# # need to specify more cores???
# fits <- ctmm.select(teleObj, guess, verbose = TRUE, cores = 2, method = "pHREML")
#
# akdeRes <- akde(teleObj, fits[[1]],
#                 weights = TRUE)
#
# ## as.sf() might be a better way of doing this???
# akde_OUT <- SpatialPolygonsDataFrame.UD(akdeRes, level.UD = 0.95)
#
# plot(akde_OUT, add = TRUE)
#
# # dBBMM -------------------------------------------------------------------
#
# library(move)
#
# moveObj <- move(x = movementData$x, y = movementData$y,
#                 time = movementData$datetime,
#                 proj = CRS(SRS_string = "EPSG:32601"))
#
# ws <- 25
# mrg <- 5
# set_grid.ext <- 2
# set_dimsize <- 400
# dbbmm <- brownian.bridge.dyn(object = moveObj,
#                              location.error = 5,
#                              margin = mrg,
#                              window.size = ws,
#                              ext = set_grid.ext,
#                              dimSize = set_dimsize,
#                              verbose = FALSE)
# library(rgeos)
#
# dbbmmSP <- as(dbbmm, "SpatialPixelsDataFrame")
# dbbmmSP_UD <- new("estUD", dbbmmSP)
# dbbmmSP_UD@vol = FALSE
# dbbmmSP_UD@h$meth = "dBBMM"
# dbbmm_UD <- getvolumeUD(dbbmmSP_UD, standardize = TRUE)
# dbbmm_OUT <- getverticeshr(dbbmm_UD, percent = 95)
#
# plot(dbbmm_OUT, add = TRUE)



