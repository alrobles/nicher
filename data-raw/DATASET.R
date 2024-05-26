## code to prepare `DATASET` dataset goes here
readr::write_csv(stack_1_12, "data-raw/stack_1_12.csv", )
write.csv(stack_1_12, file=gzfile("data-raw/stack_1_12.csv.gz"))
usethis::use_data(DATASET, overwrite = TRUE)
