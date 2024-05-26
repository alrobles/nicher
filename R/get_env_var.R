
#' get_env_var
#' Function to extract environmental values from presence points
#'
#' @param env A raster with environmental layers. Must have a projection string
#' @param df A data frame with presence points
#' @export
#'
#' @examples
#' stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
#' stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
#' get_env_var(rawSpOccPnts, stack_1_12)
get_env_var <- function(df, env){
  if(ncol(df) == 3){
    pts <- df[ , c(2, 3)]
  } else if(ncol(df) == 2) {
    pts = df
  } else {
    stop("Provide a valid data frame")
  }
  occ <- terra::extract(env, pts)
  occ[ ,-1]
}
