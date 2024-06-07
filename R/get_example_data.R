#' Function to return mammal shapefiles of Rodentia order
#' @param file An example file for test the r package
#' @export
#' @examples
#' get_example_data(file = "test")
get_example_data <- function(file){
  exampleMatch <- c("Mshp", "stack_1_12", "stack_1_12_19", "stack_1_12_crop", "test", "Mshp_test", "occpts_pellaea_pternifolia", "msample_pellaea_pternifolia", "sciurus_aberti_pts")
  if( ! file %in% exampleMatch )  {
    stop("File not found. \n
          Only Mshp, stack_1_12, stack_1_12_19, occpts_pellaea_pternifolia, msample_pellaea_pternifolia, sciurus_aberti_pts and test are valid. /n
         Check spell")
  }

  if(file == "test"){
    url_test <- "https://examples.fra1.digitaloceanspaces.com/test.txt"
    if (!url_exists(url_test)) {
      return("Can't access data.")
    } else {
      return("Data is online")
    }
  }

  #mddSpList <- get(paste0("mddSpList_", version ))
  if(file == "stack_1_12"){
    url <- "https://examples.fra1.digitaloceanspaces.com/stack_1_12.tif"
    r <- terra::rast(url)
  }

  # check if the url exists
  if(file == "stack_1_12_19"){
    url <- "https://examples.fra1.digitaloceanspaces.com/stack_1_12_19.tif"
    r <- terra::rast(url)
  }

  if(file == "Mshp"){
    url <- "https://examples.fra1.digitaloceanspaces.com/Mshp.rds"
    r <- readr::read_rds(url)
    r <- terra::unwrap(r)
  }

  if(file == "stack_1_12_crop"){
    url <- "https://examples.fra1.digitaloceanspaces.com/stack_1_12_crop.rds"
    r <- readr::read_rds(url)
    r <- terra::unwrap(r)
  }

  if(file == "Mshp_test"){
    url <- "https://examples.fra1.digitaloceanspaces.com/Mshp_test.rds"
    r <- readr::read_rds(url)
    r <- terra::unwrap(r)
  }
  # Examples Pellaea pternifolia
  if(file == "msample_pellaea_pternifolia"){
    url <- "https://examples.fra1.digitaloceanspaces.com/msample_pellaea_pternifolia.rds"
    r <- readr::read_rds(url)
    r <- terra::unwrap(r)
  }

  if(file == "occpts_pellaea_pternifolia"){
    url <- "https://examples.fra1.digitaloceanspaces.com/occpts_pellaea_pternifolia.rds"
    r <- readr::read_rds(url)
    r <- terra::unwrap(r)
  }




  if (!url_exists(url)) {
      return("Can't access data.")
    }
    #Sys.sleep(1)

  return(r)
}
