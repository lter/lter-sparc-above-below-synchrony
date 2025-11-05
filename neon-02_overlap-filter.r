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

# Identify the key columns used to identify space/time in **ALL** relevant NEON files
spacetime_vars <- c("namedLocation", "domainID", "siteID", "plotID", "endDate", "collectDate")

## -------------------------------------- ##
# Make Index of Spatiotemporal Info ----
## -------------------------------------- ##

# Identify locally-present data files
## Remove an empty file I know exists
local_all <- setdiff(x = dir(path = file.path("data", "raw_neon")),
                     y = "soil-microbe-comm-taxonomy_mct_soilSampleMetadata_16S_2016.csv")

# Remove 'variables' tables
local_data <- local_all[stringr::str_detect(string = local_all, pattern = "variables") != T]

# Identify types of tables (other than variables data)
(local_types <- unique(gsub(pattern = "plant-presence-percent-cover_|soil-microbe-comm-taxonomy_|_[[:digit:]]{4}.csv", 
                            replacement = "", x = local_data)))

# Loop across types of data
for(focal_type in local_types){
  # focal_type <- "div_10m2Data100m2Data"

  # Progress message
  message("Creating spatiotemporal inventory for ", focal_type)

  # Identify just the files of that type
  focal_files <- local_data[stringr::str_detect(string = local_data, pattern = focal_type)]

  # Generate a list for storing outputs
  out_list <- list()

  # Loop across those files
  for(focal_csv in focal_files){
    # focal_csv <- "soil-microbe-comm-taxonomy_mct_soilSampleMetadata_ITS_2016.csv"

    # Read in the data
    focal_v01 <- read.csv(file = file.path("data", "raw_neon", focal_csv))

    # Process this slightly
    focal_v02 <- focal_v01 |> 
      # Streamline this to just the spatiotemporal columns
      dplyr::select(dplyr::contains(spacetime_vars)) |> 
      # Drop non-unique rows
      dplyr::distinct() |> 
      # Add filename as a column
      dplyr::mutate(source = focal_csv, .before = dplyr::everything())

    # Add to list
    out_list[[focal_csv]] <- focal_v02

  } # Close within-type loop

  # Unlist the list to a long table
  out_v01 <- purrr::list_rbind(x = out_list)

  # Wrangle date columns based on which is in the data
  if("collectDate" %in% names(out_v01)){

    # Do needed wrangling
    out_v02 <- out_v01 |> 
      # Strip out date information
      dplyr::mutate(
        collectYear = stringr::str_sub(string = collectDate, start = 1, end = 4),
        collectMonth = stringr::str_sub(string = collectDate, start = 6, end = 7),
        collectDay = stringr::str_sub(string = collectDate, start = 9, end = 10)
      ) |> 
      # Assemble tidy date column
      dplyr::mutate(tidyDate = as.Date(paste(collectYear, collectMonth, collectDay, sep = "/"))) |> 
      # Remove superseded date column
      dplyr::select(-collectDate)

    # Handle data with 'endDate'
  } else if("endDate" %in% names(out_v01)){

    # Do needed wrangling
    out_v02 <- out_v01 |> 
      # Strip out date information
      dplyr::mutate(
        collectYear = stringr::str_sub(string = endDate, start = 1, end = 4),
        collectMonth = stringr::str_sub(string = endDate, start = 6, end = 7),
        collectDay = stringr::str_sub(string = endDate, start = 9, end = 10)
      ) |> 
      # Assemble tidy date column
      dplyr::mutate(tidyDate = as.Date(paste(collectYear, collectMonth, collectDay, sep = "/"))) |> 
      # Remove superseded date column
      dplyr::select(-endDate)

  } else { 
    message("No date wrangling performed; clarify date column")
    out_v02 <- out_v01
  }

  # Export locally
  write.csv(x = out_v02, na = '', row.names = F,
            file = file.path("data", "indices_neon", 
                             paste0("neon-spacetime-index_", focal_type, ".csv")))

} # Close among-type loop


str(focal_v02)





## -------------------------------------- ##
# Identify Overlapped Years ----
## -------------------------------------- ##

# List files of both types






# End ----
