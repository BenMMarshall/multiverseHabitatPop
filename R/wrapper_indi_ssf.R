#' Run all ssf options
#'
#' @name wrapper_indi_ssf
#' @description A
#' @param movementData must have a x and y column for locations, and a datetime
#'   column for timestamps ("%Y-%m-%d %H:%M:%S")
#' @param landscape
#' @param optionsList Must have the following items: Method_method,
#'   MethodSSF_mf, MethodSSF_sd, MethodSSF_td, MethodSSF_as
#' @return a
#'
#' @export
wrapper_indi_ssf <- function(
    allIndividualData,
    # landscape,
    optionsList
){
  
  # allIndividualData <- sampledIndividualData
  
  Method_method <- optionsList$Method_method
  MethodSSF_mf <- optionsList$MethodSSF_mf
  MethodSSF_sd <- optionsList$MethodSSF_sd
  MethodSSF_td <- optionsList$MethodSSF_td
  MethodSSF_as <- optionsList$MethodSSF_as
  
  landscape <- allIndividualData$landscape
  
  indiSSFResults <- vector("list", length = length(names(allIndividualData))-1)
  names(indiSSFResults) <- names(allIndividualData)[!names(allIndividualData) == "landscape"]
  for(indiID in names(allIndividualData)){
    if(indiID == "landscape"){
      {next}
    }
    
    # print(indiID)
    # indiID <- "simData_i001"
    movementData <- allIndividualData[[indiID]]$locations
    
    # ssf places
    listSize <- length(MethodSSF_mf) *
      length(MethodSSF_sd) *
      length(MethodSSF_td) *
      length(MethodSSF_as)
    
    listOUT <- vector("list",
                      length = listSize)
    i <- 0
    for(mf in MethodSSF_mf){
      # mf <- MethodSSF_mf[1]
      for(stepD in MethodSSF_sd){
        # stepD <- MethodSSF_sd[1]
        for(turnD in MethodSSF_td){
          # turnD <- MethodSSF_td[1]
          for(as in MethodSSF_as){
            # as <- MethodSSF_as[1]
            
            ssfOUT <- method_indi_ssf(
              movementData = movementData,
              landscape = landscape,
              methodForm = mf,
              stepDist = stepD,
              turnDist = turnD,
              availableSteps = as
            )
            
            ssfDF <- as.data.frame(summary(ssfOUT)$coef)
            method <- rep("ssf", nrow(ssfDF))
            ssfDF <- cbind(ssfDF, method)
            ssfEst <- multiverseHabitat::extract_estimate(ssfDF)
            
            optionsData <- data.frame(
              id = movementData$id[1],
              Estimate = ssfEst$Estimate,
              Lower = ssfEst$Estimate - ssfEst$SE,
              Upper = ssfEst$Estimate + ssfEst$SE,
              analysis = "ssf",
              modelFormula = mf,
              stepDist = stepD,
              turnDist = turnD,
              availablePerStep = as,
              trackFreq = allIndividualData[[indiID]]$trackFreq,
              trackDura = allIndividualData[[indiID]]$trackDura
            )
            
            ssfOUT$options <- optionsData
            
            i <- i+1
            listOUT[[i]] <- ssfOUT
            # print(i)
          } # as
        } # stepD
      } # turnD
    } # mf
    indiSSFResults[[indiID]] <- do.call(list, listOUT)
  }
  
  # allIndiSSFResults <- do.call(list, indiSSFResults)
  
  # return(allIndiSSFResults)
  return(indiSSFResults)
}
