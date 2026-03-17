#' get_ENM_par
#'
#' get ellipsoid parameters for a Ecollogical fundamental niche model.
#' @param occPts A data frame with ocurrence points. Has species name, longitude and latitude
#' @param M_shp A shape file as M hypothesis
#' @param env A raster stack
#' @param itnmax The maximum number of iterations.
#' @param fastmode logical if is true the L-BFGS-B method is by default
#' @param return_best logical. If true return the best optimization. If false return all the possible optimizations
#' @param method Method selection. By default, optimization of the negloglike function with lower bound
#' @return a list of parameters
#' @export
#'
#' @examples
#' M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
#' Mshp <- terra::unwrap(readr::read_rds(M_path))
#' stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
#' stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
#' get_ENM_par(rawSpOccPnts, stack_1_12, Mshp, method = "mahalanobis")
get_ENM_par <- function(occPts, env,
                        M_shp = NULL,
                        method = c("unbound", "bound", "mahalanobis" ),
                        fastmode = TRUE,
                        itnmax = 100,
                        return_best = TRUE){
  method <- match.arg(method)
  occPts <- get_env_var(occPts, env)

  if(is.null(M_shp)){
    warning("No shp provided. Mahalanobis method by default")
    method = "mahalanobis"
  }

  if(method == "mahalanobis"){
    pars <- get_ellip_par(occPts)
  } else if(method == "unbound"){
    sampleM <- sam_polyM(M_shp, env)
    pars <- get_negloglike_optimr_par(occPts, sampleM, lower = FALSE, fastmode, itnmax )
  } else if (method == "bound"){
    sampleM <- sam_polyM(M_shp, env)
    pars <- get_negloglike_optimr_par(occPts, sampleM, lower = TRUE, fastmode, itnmax )
  }
  return(ellipse(pars))
}
