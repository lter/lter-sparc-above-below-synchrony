## ------------------------------------------------------- ##
# NEON Data - Filter to Overlapping Space/Time
## ------------------------------------------------------- ##
# Purpose:
## Filter downloaded NEON data to only data in same time/space between plant & microbe data

# Load libraries
librarian::shelf(tidyverse)

# Do setup tasks
source(file = file.path("00_setup.r"))

# Clear environment
rm(list = ls()); gc()

## -------------------------------------- ##
# Prepare to Check Space/Time Overlap ----
## -------------------------------------- ##

# What files do we have locally already?
local_plant <- dir(path = file.path("data", "raw_neon"), 
                   pattern = "plant-presence-percent-cover_")
local_micro <- dir(path = file.path("data", "raw_neon"), 
                   pattern = "soil-microbe-comm-taxonomy_")

# Identify types of tables we have for each
(type_plant <- unique(gsub(pattern = "plant-presence-percent-cover_|_[[:digit:]]{4}.csv", 
                           replacement = "", x = local_plant)))
(type_micro <- unique(gsub(pattern = "soil-microbe-comm-taxonomy_|_[[:digit:]]{4}.csv", 
                           replacement = "", x = local_micro)))

# Read in one variable table for each
vars_plant_v01 <- read.csv(file = file.path("data", "raw_neon", "plant-presence-percent-cover_variables_10058_2021.csv"))
vars_micro_v01 <- read.csv(file = file.path("data", "raw_neon", "soil-microbe-comm-taxonomy_variables_10081_2021.csv"))

# Subset to only types of data that we have locally
vars_plant_v02 <- dplyr::filter(vars_plant_v01, table %in% type_plant)
vars_micro_v02 <- dplyr::filter(vars_micro_v01, table %in% type_micro)

# Keep only columns found in all included tables
vars_plant_v03 <- vars_plant_v02 |> 
  ## Count tables per variable
  dplyr::group_by(fieldName) |> 
  dplyr::mutate(fieldCount = length(unique(table))) |> 
  dplyr::ungroup() |> 
  ## Filter to only ones where all tables have the column
  dplyr::filter(fieldCount == length(type_plant) - 1) |> 
  ## Streamline the data
  dplyr::select(fieldName) |> 
  dplyr::distinct()
## Same for microbes
vars_micro_v03 <- vars_micro_v02 |> 
  ## Count tables per variable
  dplyr::group_by(fieldName) |> 
  dplyr::mutate(fieldCount = length(unique(table))) |> 
  dplyr::ungroup() |> 
  ## Filter to only ones where all tables have the column
  dplyr::filter(fieldCount == length(type_plant) - 1) |> 
  ## Streamline the data
  dplyr::select(fieldName) |> 
  dplyr::distinct()

# Identify just the variables that are in all of the data tables
## And remove some that we don't care about
vars_both <- intersect(x = vars_plant_v03$fieldName, y = vars_micro_v03$fieldName)

## -------------------------------------- ##
# Identify Overlapped Years ----
## -------------------------------------- ##

# List files of both types






# End ----
