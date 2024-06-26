% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sam_polyM.R
\name{sam_polyM}
\alias{sam_polyM}
\title{Sample environmental points from M hypothesis}
\usage{
sam_polyM(M.shp, env, N = 100)
}
\arguments{
\item{M.shp}{A shapefile with an M hypothesis of the species}

\item{env}{an stack of raster with environmental variables to extract information}

\item{N}{a numeric with the number of points of the sample, default is 10 000}
}
\value{
a data frame with a sample of environmental values inside the M polygon
}
\description{
Sample environmental points from M hypothesis
}
\examples{
M_path <- system.file("extdata", "Mshp_test.rds", package="nicher")
Mshp <- terra::unwrap(readr::read_rds(M_path))
stack_path <- system.file("extdata", "stack_1_12_crop.rds", package="nicher")
stack_1_12 <- terra::unwrap(readr::read_rds(stack_path))
sam_polyM(M.shp = Mshp, N = 5, env = stack_1_12)
}
