
targets::tar_load("modelsBrms")

# brmsOUT <- poisBrms

names(modelsBrms)

# compile all R2 from all models ------------------------------------------
# library(performance)
# as.data.frame(wrsfBrms$modOUT_dEst_r2)
# data.frame(
#   R2 = c(unlist(
#     wrsfBrms$modOUT_dEst_r2[1]),
#     unlist(wrsfBrms$modOUT_dEst_r2[2])
#     ),
# SE = c(attr(wrsfBrms$modOUT_dEst_r2,"SE")$R2_Bayes,
#        attr(wrsfBrms$modOUT_dEst_r2,"SE")$R2_Bayes_marginal),
# rbind(attr(wrsfBrms$modOUT_dEst_r2,"CI")$R2_Bayes,
#       attr(wrsfBrms$modOUT_dEst_r2,"CI")$R2_Bayes_marginal),
# Component = c("conditional", "marginal")
# )

modelsBrms


# Get all r2 --------------------------------------------------------------

r2Outputs <- lapply(names(modelsBrms), function(x){
  # x <- names(modelsBrms)[1]
  model <- modelsBrms[[x]]
  r2OUT <- performance::r2_bayes(model)
  
  r2OUT
  
  r2df_1 <- data.frame(
    R2 = c(unlist(
      r2OUT[1]),
      unlist(r2OUT[2])
    ),
    SE = c(attr(r2OUT,"SE")$R2_Bayes,
           attr(r2OUT,"SE")$R2_Bayes_marginal),
    rbind(attr(r2OUT,"CI")$R2_Bayes,
          attr(r2OUT,"CI")$R2_Bayes_marginal),
    Component = c("conditional", "marginal")
  )
  r2df_1$model <- x
  
  return(r2df_1)
  
})
r2Outputs <- do.call(rbind, r2Outputs)

write.csv(r2Outputs, file = here::here("data", "brmsR2Results.csv"),
          row.names = FALSE)

# get all betas -----------------------------------------------------------

modelsBrms

betasOutputs <- lapply(names(modelsBrms), function(x){
  model <- modelsBrms[[x]]
  if(class(model) == "brmsfit"){
    out <- ggdist::median_hdci(tidybayes::gather_draws(model,
                                                       `b_.*`, regex = TRUE),
                               .width = c(0.95))
    out$model <- x
    return(out)
  } else (
    return(NULL)
  )
})
betasOutputs <- do.call(rbind, betasOutputs)

write.csv(betasOutputs, file = here::here("data", "brmsEstResults.csv"),
          row.names = FALSE)

extractedList <- list(
  "r2Outputs" = r2Outputs,
  "betasOutputs" = betasOutputs)

return(extractedList)