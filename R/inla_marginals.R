#' Extract marginals from INLA
#'
#' @name inla_mmarginal
#' @description abc.
#' @param inlaOutput The output from an INLA model.
#' @return abc
#'
#' @export
inla_mmarginal <- function(inlaOutput){ 
  results <- sapply(inlaOutput$marginals.hyperpar, 
                    function(y) 
                      inla.mmarginal(inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mode of variance"))
  results
}

#' @export
inla_emarginal <- function(inlaOutput){ 
  results <- sapply(inlaOutput$marginals.hyperpar, 
                    function(y) 
                      inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, y)))
  
  names(results) <- sapply(as.vector(as.character(names(results))), function(y) gsub("Precision", x=y, "Mean of variance"))
  results
}