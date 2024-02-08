#' Extract all required values for area methods
#'
#' @name area_based_extraction
#' @description A
#' @param allIndividualData The full list of all simulated individuals, subsampled into regimes
#' @param optionsList 
#' @return a
#'
#' @export
area_based_extraction <- function(allIndividualData, optionsList){
  # allIndividualData <- sampDuraFreqData_60_1
  areaMethod <- optionsList$areaMethod
  areaContour <- optionsList$areaContour
  Method_ap <- optionsList$Method_ap
  Method_sp <- optionsList$Method_sp
  
  landscape <- list()
  landscape[[1]] <- allIndividualData$landscape$classRaster
  landscape[[2]] <- allIndividualData$landscape$classRasterScram
  names(landscape) <- c("classRaster", "classRasterScram")
  
  # Loop to create individual polygons --------------------------------------
  
  i <- 0
  usedAvailList <- vector("list",
                          length = 
                            length(areaMethod) *
                            length(landscape) *
                            length(areaContour) *
                            length(Method_sp) *
                            length(Method_ap))
  
  for(method in areaMethod){
    
    resourceList <- vector("list", length = length(names(allIndividualData)[-1]))
    names(resourceList) <- names(allIndividualData)[-1]
    for(indiID in names(allIndividualData)){
      if(indiID == "landscape"){
        {next}
      }
      
      movementData <- allIndividualData[[indiID]]$locations
      
      baa <- build_available_area(movementData = movementData,
                                  method = method,
                                  SRS_string = "EPSG:32601",
                                  dBBMMsettings = NULL)
      
      resourceList[[indiID]] <- baa
      
      print(paste0(indiID, " --- Area Resource --- ", method))
    }
    
    for(contour in areaContour){
      print(contour)
      
      polygonList <- vector("list", length = length(names(allIndividualData)[-1]))
      names(polygonList) <- names(allIndividualData)[-1]
      for(indiID in names(allIndividualData)){
        if(indiID == "landscape"){
          {next}
        }
        
        print(indiID)
        
        bap <- build_available_polygon(areaResource = resourceList[[indiID]],
                                       method = method,
                                       contour = contour,
                                       SRS_string = "EPSG:32601")
        
        
        polygonList[[indiID]] <- bap
        
        print(paste0(indiID, " --- Polygon created --- ", method, " - ", contour))
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
            
            for(land in names(landscape)){
              # land <- names(landscape)
              # landscape[[land]]
              movementData <- allIndividualData[[indiID]]$locations
              
              indiPolygon <- polygonList[[indiID]]
              
              # generate points based on the availableArea and the number of points
              suppressWarnings({
                availPoints <- sp::spsample(indiPolygon,
                                            n = nrow(movementData) * aPoints,
                                            type = ifelse(spSamp == "rd", "random", "stratified"))
              })
              
              # extract the habitat types each point is located within
              availValues <- raster::extract(landscape[[land]], availPoints)
              
              print("--- Available / Used start")
              availValues_DF <- data.frame(rbind(table(availValues)))
              names(availValues_DF) <- sub("X", "c", names(availValues_DF))
              
              suppressWarnings({
                usedValues <- raster::extract(landscape[[land]],
                                              sp::SpatialPoints(movementData[,c("x", "y")],
                                                                sp::CRS(SRS_string = "EPSG:32601")))
              })
              usedValues_DF <- data.frame(rbind(table(usedValues)))
              names(usedValues_DF) <- sub("X", "c", names(usedValues_DF))
              
              print("-- Break 1")
              
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
              
              print("-- Break 2")
              ## TYPE II ##
              suppressWarnings({
                availPopPoints <- sp::spsample(popPolygon,
                                               n = nrow(movementData) * aPoints,
                                               type = ifelse(spSamp == "rd", "random", "stratified"))
              })
              
              availPopValues <- raster::extract(landscape[[land]], availPopPoints)
              
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
              
              print("-- Break 3")
              
              usedValues_DF <- usedValues_DF[,sort(names(usedValues_DF))]
              availValues_DF <- availValues_DF[,sort(names(availValues_DF))]
              availPopValues_DF <- availPopValues_DF[,sort(names(availPopValues_DF))]
              
              names(usedValues_DF) <- paste0("used_", names(usedValues_DF))
              names(availValues_DF) <- paste0("avail_", names(availValues_DF))
              names(availPopValues_DF) <- paste0("avail_", names(availPopValues_DF))
              
              usedAvailable <- cbind(usedValues_DF, availValues_DF)
              usedAvailable$id <- movementData$id[1]
              usedAvailable$classLandscape <- land
              usedAvailable$type <- "III"
              usedAvailable$method <- method
              usedAvailable$contour <- contour
              usedAvailable$availablePoints <- aPoints
              usedAvailable$samplingPattern <- spSamp
              usedAvailable$trackFreq <- allIndividualData[[indiID]]$trackFreq
              usedAvailable$trackDura <- allIndividualData[[indiID]]$trackDura
              
              usedAvailable
              
              usedAvailablePop <- cbind(usedValues_DF, availPopValues_DF)
              usedAvailablePop$id <- movementData$id[1]
              usedAvailablePop$classLandscape <- land
              usedAvailablePop$type <- "II"
              usedAvailablePop$method <- method
              usedAvailablePop$contour <- contour
              usedAvailablePop$availablePoints <- aPoints
              usedAvailablePop$samplingPattern <- spSamp
              usedAvailablePop$trackFreq <- allIndividualData[[indiID]]$trackFreq
              usedAvailablePop$trackDura <- allIndividualData[[indiID]]$trackDura
              
              print("-- Break 4")
              
              print(usedAvailable)
              print(usedAvailablePop)
              
              usedAvailableAll <- rbind(usedAvailable, usedAvailablePop)
              
              i <- i+1
              
              usedAvailList[[i]] <- usedAvailableAll
            } # landscape loop
            
          } # indiID
        } # aPoints
      } # spSamp
    } # area contour
  } # area method
  usedAvailFull <- do.call(rbind, usedAvailList)
  
  return(usedAvailFull)
}