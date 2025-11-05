## ------------------------------------------------------- ##
# NEON Data - Download Above/Below Files
## ------------------------------------------------------- ##
# Purpose:
## Identify & download desired NEON data
## Limited attempt at filtering made so data are quite bulky and this script is time-consuming

# Load libraries
librarian::shelf(tidyverse, neonUtilities)

# Do setup tasks
source(file = file.path("00_setup.r"))

# Clear environment
rm(list = ls()); gc()

# Grab all custom functions
purrr::walk(.x = dir(path = file.path("tools"), pattern = "fxn_"),
            .f = ~ source(file.path("tools", .x)))

## -------------------------------------- ##
# Download Plant Presence/Percent Cover Data ----
## -------------------------------------- ##

# Download plant precence & percent cover data for all years / sites
download_neon(dpi = "DP1.10058.001", 
              ## Start year identified from NEON catalog
              ## End year is final year of non-provisional data
              start_yr = 2013, end_yr = 2023,
              prefix = "plant-presence-percent-cover_", 
              dest = file.path("data", "raw_neon"),
              wanted_files = c("div_1m2Data", "div_10m2Data100m2Data", "variables_10058"))

## -------------------------------------- ##
# Download Soil Microbe Taxonomy Data ----
## -------------------------------------- ##

# Downlaod soil microbe community taxonomy data for all years / sites
download_neon(dpi = "DP1.10081.002", 
              ## Start year identified from NEON catalog
              start_yr = 2016, end_yr = 2025,
              prefix = "soil-microbe-comm-taxonomy_", 
              dest = file.path("data", "raw_neon"),
              wanted_files = c("div_1m2Data", "div_10m2Data100m2Data", "variables_10058"))

# End ----
