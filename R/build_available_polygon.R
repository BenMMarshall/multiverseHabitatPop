#' Generate a polygon for use as an available area
#'
#' @name build_available_polygon
#' @description A
#' @param areaResource an output from build_available_areas that is used to
#'   extract a polygon.
#' @param method "MCP", "KDE_LSCV", "KDE_href", "AKDE", "dBBMM"
#' @param contour Numeric, > 0 and < 100. Can be a vector.
#' @param SRS_string = "EPSG:32601"
#' @return A polygon.
#'
#' @export
build_available_polygon <- function(areaResource,
                                    method = c("MCP", "KDElscv", "KDEhref", "AKDE", "dBBMM"),
                                    contour,
                                    SRS_string = "EPSG:32601"){


  if(method == "MCP"){

    poly_OUT <- adehabitatHR::mcp(areaResource, percent = contour, unin = "m",
                                  unout = "m2")

  } else if(method == "KDElscv"){

    poly_OUT <- adehabitatHR::getverticeshr(areaResource, contour)

  } else if(method == "KDEhref"){

    poly_OUT <- adehabitatHR::getverticeshr(areaResource, contour)

  } else if(method == "AKDE"){

    print(contour)
    print(class(areaResource) == "UD")
    akdePoly <- ctmm::SpatialPolygonsDataFrame.UD(object = areaResource, level.UD = contour/100)
    # just get the point estimate
    poly_OUT <- akdePoly[akdePoly$name ==
                           akdePoly$name[stringr::str_detect(akdePoly$name, "est")],]

  } else if(method == "dBBMM"){

    suppressWarnings({
      dbbmmSP <- as(areaResource, "SpatialPixelsDataFrame")
      dbbmmSP_UD <- new(getClass("estUD", where = "adehabitatHR"), dbbmmSP)
      dbbmmSP_UD@vol = FALSE
      dbbmmSP_UD@h$meth = "dBBMM"
      dbbmm_UD <- adehabitatHR::getvolumeUD(dbbmmSP_UD, standardize = TRUE)

      poly_OUT <- adehabitatHR::getverticeshr(dbbmm_UD, percent = contour)
    })

  }

  return(poly_OUT)

}

