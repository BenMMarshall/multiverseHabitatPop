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
  
  # allIndividualData <- sampDuraFreqData_15_1
    
  Method_method <- optionsList_sff$Method_method
  MethodSSF_mf <- optionsList_sff$MethodSSF_mf
  MethodSSF_sd <- optionsList_sff$MethodSSF_sd
  MethodSSF_td <- optionsList_sff$MethodSSF_td
  MethodSSF_as <- optionsList_sff$MethodSSF_as
    
  landscape <- allIndividualData$landscape
  
  indiSSFResults <- vector("list", length = length(names(allIndividualData))-1)
  names(indiSSFResults) <- names(allIndividualData)[-1]
  for(indiID in names(allIndividualData)){
    if(indiID == "landscape"){
      {next}
    }
    
    print(indiID)
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
      
      for(sd in MethodSSF_sd){
        
        for(td in MethodSSF_td){
          
          for(as in MethodSSF_as){
            
            ssfOUT <- method_indi_ssf(
              movementData = movementData,
              landscape = landscape,
              methodForm = mf,
              stepDist = sd,
              turnDist = td,
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
              stepDist = sd,
              turnDist = td,
              availablePerStep = as,
              trackFreq = allIndividualData[[indiID]]$trackFreq,
              trackDura = allIndividualData[[indiID]]$trackDura
            )
            
            ssfOUT$options <- optionsData
            
            i <- i+1
            listOUT[[i]] <- ssfOUT
            # print(i)
          } # as
        } # sd
      } # td
    } # mf
    indiSSFResults[[indiID]] <- do.call(list, listOUT)
  }
  
  # allIndiSSFResults <- do.call(list, indiSSFResults)
  
  # return(allIndiSSFResults)
  return(indiSSFResults)
}
