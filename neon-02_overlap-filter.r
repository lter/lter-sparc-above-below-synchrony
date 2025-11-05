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
# Make Spatiotemporal Indices ----
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

# Make a list for storing type outputs
type_list_v01 <- list()

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
      dplyr::select(-dplyr::contains("subplot")) |> 
      # Drop non-unique rows
      dplyr::distinct() |> 
      # Add file name/type as columns
      dplyr::mutate(type = focal_type, 
                    source = focal_csv, 
                    .before = dplyr::everything())

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
      dplyr::mutate(tidyDate = as.Date(paste(collectYear, collectMonth, collectDay, sep = "/")),
                    tidyYearMonth = paste0(collectYear, "-", collectMonth)) |> 
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
      dplyr::mutate(tidyDate = as.Date(paste(collectYear, collectMonth, collectDay, sep = "/")),
                    tidyYearMonth = paste0(collectYear, "-", collectMonth)) |> 
      # Remove superseded date column
      dplyr::select(-endDate)

  } else { 
    message("No date wrangling performed; clarify date column")
    out_v02 <- out_v01
  }

  # And add to list
  type_list_v01[[focal_type]] <- out_v02

} # Close among-type loop

## -------------------------------------- ##
# Assess Year Overlap ----
## -------------------------------------- ##

# Unlist the space/time indices list
overlap_v01 <- purrr::list_rbind(x = type_list_v01)

# Check structure
dplyr::glimpse(overlap_v01)

# Identify years found in all data files
good_years <- overlap_v01 |> 
  # Pare down to only data we need right now
  dplyr::select(type, collectYear) |> 
  dplyr::distinct() |> 
  # Count types of data per year
  dplyr::group_by(collectYear) |> 
  dplyr::summarize(type_ct = length(unique(type)),
                   type_names = paste(unique(type), collapse = "; "),
                   .groups = "keep") |> 
  dplyr::ungroup() |> 
  # Filter to only years with all types of data
  dplyr::filter(type_ct == length(unique(overlap_v01$type))) |> 
  # Grab just the years
  dplyr::pull(collectYear)

# What years had good overlap?
sort(unique(good_years))

# Make a new list
type_list_v02 <- list()

# Loop across bits of the previous type list
for(k in seq_along(type_list_v01)){

  # Grab that list element
  year_v01 <- type_list_v01[[k]]

  # Pare down to only good years
  year_v02 <- dplyr::filter(year_v01, collectYear %in% good_years)

  # Add to the new list
  type_list_v02[[k]] <- year_v02
}

## -------------------------------------- ##
# Assess Space Overlap ----
## -------------------------------------- ##

# Unlist the space/time indices list
overlap_v02 <- purrr::list_rbind(x = type_list_v02)

# Check structure
dplyr::glimpse(overlap_v02)

# Identify years found in all data files
good_spaces <- overlap_v02 |> 
  # Pare down to only data we need right now
  dplyr::select(type, dplyr::ends_with("ID")) |> 
  dplyr::distinct() |> 
  # Pivot into ultra long format
  tidyr::pivot_longer(cols = dplyr::ends_with("ID"),
                      names_to = "space_cat", values_to = "space_level") |> 
  # Count types of data per space category and level
  dplyr::group_by(space_cat, space_level) |> 
  dplyr::summarize(type_ct = length(unique(type)),
                   type_names = paste(unique(type), collapse = "; "),
                   .groups = "keep") |> 
  dplyr::ungroup() |> 
  # Filter to only years with all types of data
  dplyr::filter(type_ct == length(unique(overlap_v02$type)))
  
# Check structure
dplyr::glimpse(good_spaces)

# Make a new list
type_list_v03 <- list()

# Loop across bits of the previous type list
for(k in seq_along(type_list_v02)){

  # Grab that list element
  space_v01 <- type_list_v02[[k]]

  # Pare down to only good spaces
  space_v02 <- space_v01 |> 
    dplyr::filter(domainID %in% dplyr::filter(good_spaces, space_cat == "domainID")$space_level) |> 
    dplyr::filter(siteID %in% dplyr::filter(good_spaces, space_cat == "siteID")$space_level) |> 
    dplyr::filter(plotID %in% dplyr::filter(good_spaces, space_cat == "plotID")$space_level)

  # Add to the new list
  type_list_v03[[k]] <- space_v02
}





  # Export locally
 #  write.csv(x = out_v02, na = '', row.names = F,
  #           file = file.path("data", "indices_neon", 
   #                           paste0("neon-spacetime-index_", focal_type, ".csv")))
  


## -------------------------------------- ##
## -------------------------------------- ##

# List files of both types






# End ----
