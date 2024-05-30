#' get_ENM_par
#'
#' get ellipsoid parameters for a Ecollogical fundamental niche model.
#' @param occPts A data frame with ocurrence points. Has species name, longitude and latitude
#' @param M_shp A shape file as M hypothesis
#' @param env A raster stack
#' @param method Method selection. By default, optimization of the negloglike function with lower bound
#'
#' @return a list of parameters
#' @export
#'
#' @examples
#' M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
#' Mshp <- terra::unwrap(readr::read_rds(M_path))
#' stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
#' stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
#' get_ENM_par(rawSpOccPnts, stack_1_12, Mshp, method = "bound")
get_ENM_par <- function(occPts, env, M_shp = NULL,
                        method = c("bound", "mahalanobis", "vmmin", "Nelder-Mead")){
  method <- match.arg(method)
  occPts <- get_env_var(occPts, env)

  if(is.null(M_shp)){
    warning("No shp provided. Mahalanobis method by default")
    method = "mahalanobis"
  }

  if(method == "mahalanobis"){
    pars <- get_ellip_par(occPts)
  }

  if(method == "Nelder-Mead"){
    sampleM <- sam_polyM(M_shp, env)
    pars <- get_negloglike_optim_par(occPts, sampleM )
  }
  if(method == "vmmin"){
    sampleM <- sam_polyM(M_shp, env)
    pars <- get_negloglike_optimr_par(occPts, sampleM )
  }
  if (method == "bound"){
    sampleM <- sam_polyM(M_shp, env)
    pars <- get_negloglike_optimr_par(occPts, sampleM, lower = TRUE )
  }
  return(ellipse(pars))
}
