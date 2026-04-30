## code to prepare `DATASET` dataset goes here
library(tidyverse)
usethis::use_data(DATASET, overwrite = TRUE)

env_m <- readr::read_csv("data-raw/vicugna_10K_Mpoints_lonlat_bio1_bio12_7nov2025.csv") %>%
  select(bio1, bio12)

coords_m <- readr::read_csv("data-raw/vicugna_10K_Mpoints_lonlat_bio1_bio12_7nov2025.csv") %>%
  select(long, lat) %>%
  rename(lon = long)


env_occ <- readr::read_csv("data-raw/vicugna_occurrences_lonlat_bio1_bio12_7nov2025.csv") %>%
  select(bio1, bio12)

coords_occ <- readr::read_csv("data-raw/vicugna_occurrences_lonlat_bio1_bio12_7nov2025.csv") %>%
  select(long, lat) %>%
  rename(lon = long)


example_vicugna <- list(env_m = env_m,
     env_occ = env_occ,
     coords_m = coords_m,
     coords_occ = coords_occ)

usethis::use_data(example_vicugna)


env_m %>%
  ggplot() + geom_point(aes(long, lat, col = bio12), alpha = 0.01) +
  geom_point(data = env_occ, aes(long, lat))+
  theme_classic()

env_occ %>%
  ggplot() +  geom_histogram(aes(bio1))

env_m %>%
  ggplot() +  geom_histogram(aes(bio1))
