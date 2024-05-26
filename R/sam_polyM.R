#' Sample environmental points from M hypothesis
#'
#' @param M.shp A shapefile with an M hypothesis of the species
#' @param N a numeric with the number of points of the sample, default is 10 000
#' @param env an stack of raster with environmental variables to extract information
#' @return a data frame with a sample of environmental values inside the M polygon
#' @export
#'
#' @examples
#' M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
#' Mshp <- terra::unwrap(readr::read_rds(M_path))
#' stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
#' stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
#' sam_polyM(M.shp = Mshp, N = 5, env = stack_1_12)

sam_polyM <- function(M.shp, env, N = 100){

  # crop and mask the environmental layers with the M polygon
  crop.M <- terra::crop(x = env, y = M.shp, mask = TRUE)
  # get ride of cells with NA values
  ind <- which(!is.na(crop.M[[1]][]))
  # get a random sample of indices
  sam <- sample(ind, N, replace = TRUE)
  # choose the points corresponding to the selected indices
  Mpnts <- crop.M[][sam, ]
  return(Mpnts)
}
