
targets::tar_load("sampDuraFreqData_15_1")

landscape <- sampDuraFreqData_15_1$landscape$classRaster

# c("MCP", "KDElscv", "KDEhref", "AKDE", "dBBMM"),

# Settings for variations / choices ---------------------------------------

# spSamp <- "random"
# availablePointsPer <- 2

# areaMethod = c("MCP", "AKDE")
areaMethod = "MCP"
areaContour = c(95, 99)
Method_ap = as.integer(round(exp(seq(log(1), log(10), length.out = 2)), digits = 1))
Method_sp = c("rd", "st")
Method_we = exp(seq(log(1), log(10000000), length.out = 2))


movementData <- sampDuraFreqData_15_1$simData_i001$locations

allIndividualData <- sampDuraFreqData_15_1


# Loop to create individual polygons --------------------------------------

i <- 0
usedAvailList <- vector("list",
                        length = 
                          length(areaMethod) *
                          length(areaContour) *
                          length(Method_sp) *
                          length(Method_ap))
for(method in areaMethod){
  
  print(method)
  
  for(contour in areaContour){
    print(contour)
    
    polygonList <- vector("list", length = length(names(allIndividualData)[-1]))
    names(polygonList) <- names(allIndividualData)[-1]
    for(indiID in names(allIndividualData)){
      
      if(indiID == "landscape"){
        {next}
      }
      
      print(indiID)
      
      movementData <- allIndividualData[[indiID]]$locations
      
      baa <- build_available_area(movementData = movementData,
                                  method = method,
                                  SRS_string = "EPSG:32601",
                                  dBBMMsettings = NULL)
      
      bap <- build_available_polygon(areaResource = baa,
                                     method = method,
                                     contour = contour,
                                     SRS_string = "EPSG:32601")
      
      print("polygon complete")
      
      polygonList[[indiID]] <- bap
    } # end indiID polygon generation
    polygonList
    # Loop to create TypeII polygon -------------------------------------------
    
    popPolygon <- polygonList[[1]]
    for(p in names(polygonList[-1])){
      popPolygon <- rgeos::gUnion(popPolygon, polygonList[[p]])
    }
    popPolygon
    print("pop polygon complete")

    for(spSamp in Method_sp){
      print(spSamp)
      for(aPoints in Method_ap){
        print(aPoints)
        
        
        # Individual use / available ----------------------------------------------
        # usedAvailList <- vector("list", length = length(names(allIndividualData)[-1]))
        # names(usedAvailList) <- names(allIndividualData)[-1]
        for(indiID in names(allIndividualData)){
          
          if(indiID == "landscape"){
            {next}
          }
          
          print(indiID)
          
          movementData <- allIndividualData[[indiID]]$locations
          
          indiPolygon <- polygonList[[indiID]]
          
          # generate points based on the availableArea and the number of points
          suppressWarnings({
            availPoints <- sp::spsample(indiPolygon,
                                        n = nrow(movementData) * aPoints,
                                        type = ifelse(spSamp == "rd", "random", "stratified"))
          })
          
          # extract the habitat types each point is located within
          availValues <- raster::extract(landscape, availPoints)
          
          availValues_DF <- data.frame(rbind(table(availValues)))
          names(availValues_DF) <- sub("X", "c", names(availValues_DF))
          
          suppressWarnings({
            usedValues <- raster::extract(landscape,
                                          sp::SpatialPoints(movementData[,c("x", "y")],
                                                            sp::CRS(SRS_string = "EPSG:32601")))
          })
          usedValues_DF <- data.frame(rbind(table(usedValues)))
          names(usedValues_DF) <- sub("X", "c", names(usedValues_DF))
          
          aClass <- names(availValues_DF)
          uClass <- names(usedValues_DF)
          
          if(length(aClass) > length(uClass)){
            toAdd <- as.data.frame(matrix(0, nrow = 1, ncol = length(aClass[!aClass %in% uClass])))
            names(toAdd) <- aClass[!aClass %in% uClass]
            usedValues_DF <- cbind(usedValues_DF, toAdd)
          } else if(length(uClass) > length(aClass)){
            toAdd <- as.data.frame(matrix(0, nrow = 1, ncol = length(uClass[!uClass %in% aClass])))
            names(toAdd) <- uClass[!uClass %in% aClass]
            availValues_DF <- cbind(availValues_DF, toAdd)
          }
          
          ## TYPE II ##
          suppressWarnings({
            availPopPoints <- sp::spsample(popPolygon,
                                           n = nrow(movementData) * aPoints,
                                           type = ifelse(spSamp == "rd", "random", "stratified"))
          })
          
          availPopValues <- raster::extract(landscape, availPopPoints)
          
          availPopValues_DF <- data.frame(rbind(table(availPopValues)))
          names(availPopValues_DF) <- sub("X", "c", names(availPopValues_DF))
          
          aPopClass <- names(availPopValues_DF)
          if(length(aPopClass) > length(uClass)){
            toAdd <- as.data.frame(matrix(0, nrow = 1, ncol = length(aPopClass[!aPopClass %in% uClass])))
            names(toAdd) <- aPopClass[!aPopClass %in% uClass]
            usedValues_DF <- cbind(usedValues_DF, toAdd)
          } else if(length(uClass) > length(aPopClass)){
            toAdd <- as.data.frame(matrix(0, nrow = 1, ncol = length(uClass[!uClass %in% aPopClass])))
            names(toAdd) <- uClass[!uClass %in% aPopClass]
            availPopValues_DF <- cbind(availPopValues_DF, toAdd)
          }
          
          usedValues_DF <- usedValues_DF[,sort(names(usedValues_DF))]
          availValues_DF <- availValues_DF[,sort(names(availValues_DF))]
          availPopValues_DF <- availPopValues_DF[,sort(names(availPopValues_DF))]
          
          names(usedValues_DF) <- paste0("used_", names(usedValues_DF))
          names(availValues_DF) <- paste0("avail_", names(availValues_DF))
          names(availPopValues_DF) <- paste0("avail_", names(availPopValues_DF))
          
          usedAvailable <- cbind(usedValues_DF, availValues_DF)
          usedAvailable$id <- indiID
          usedAvailable$type <- "III"
          usedAvailable$method <- method
          usedAvailable$contour <- contour
          usedAvailable$aPoints <- aPoints
          usedAvailable$spSamp <- spSamp
          
          usedAvailable
          
          usedAvailablePop <- cbind(usedValues_DF, availPopValues_DF)
          usedAvailablePop$id <- indiID
          usedAvailablePop$type <- "II"
          usedAvailablePop$method <- method
          usedAvailablePop$contour <- contour
          usedAvailablePop$aPoints <- aPoints
          usedAvailablePop$spSamp <- spSamp
          
          usedAvailableAll <- rbind(usedAvailable, usedAvailablePop)
          
          i <- i+1
          
          usedAvailList[[i]] <- usedAvailableAll
          
        }
      } # aPoints
    } # spSamp
  } # area contour
} # area method
do.call(rbind, usedAvailList)



# -------------------------------------------------------------------------

# generate points based on the availableArea and the number of points
### POSSIBLE NEW NODE, RANDOM VERSUS SYSTEMATIC???
suppressWarnings({
  availPoints <- sp::spsample(bap,
                              n = nrow(movementData) * availablePointsPer,
                              type = ifelse(spSamp == "rd", "random", "stratified"))
})

# extract the habitat types each point is located within
availValues <- raster::extract(landscape, availPoints)

availValues_DF <- data.frame(rbind(table(availValues)))
names(availValues_DF) <- sub("X", "c", names(availValues_DF))

suppressWarnings({
  usedValues <- raster::extract(landscape,
                                sp::SpatialPoints(movementData[,c("x", "y")],
                                                  sp::CRS(SRS_string = "EPSG:32601")))
})
usedValues_DF <- data.frame(rbind(table(usedValues)))
names(usedValues_DF) <- sub("X", "c", names(usedValues_DF))

aClass <- names(availValues_DF)
uClass <- names(usedValues_DF)

if(length(aClass) > length(uClass)){
  
  toAdd <- as.data.frame(matrix(0, nrow = 1, ncol = length(aClass[!aClass %in% uClass])))
  names(toAdd) <- aClass[!aClass %in% uClass]
  usedValues_DF <- cbind(usedValues_DF, toAdd)
  
} else if(length(uClass) > length(aClass)){
  
  toAdd <- as.data.frame(matrix(0, nrow = 1, ncol = length(uClass[!uClass %in% aClass])))
  names(toAdd) <- uClass[!uClass %in% aClass]
  availValues_DF <- cbind(availValues_DF, toAdd)
  
}

usedValues_DF <- usedValues_DF[,sort(names(usedValues_DF))]
availValues_DF <- availValues_DF[,sort(names(availValues_DF))]

names(usedValues_DF) <- paste0("used_", names(usedValues_DF))
names(availValues_DF) <- paste0("avail_", names(availValues_DF))

usedAvailable <- cbind(usedValues_DF, availValues_DF)
usedAvailable$id <- indiID
usedAvailable$type <- "III"

# usedAvailable$method <- "III"
# usedAvailable$contour <- "III"
# usedAvailable$ap <- "III"
# usedAvailable$spSamp <- spSamp

# so because the difference is use the available, we can do III design with the
# II set up just with different availabilities. Is in a try() function because
# of the instances where habitats are used, but not available bceause of the
# random sampling of available points