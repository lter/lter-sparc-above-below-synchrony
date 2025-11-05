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

# Grab current year & month
this_year <- stringr::str_sub(Sys.Date(), start = 1, end = 4)
this_month <- stringr::str_sub(Sys.Date(), start = 1, end = 7)

## -------------------------------------- ##
# Identify & Download Plant Data ----
## -------------------------------------- ##

# Define some useful objects for making the loop more flexible
focal_dpi <- "DP1.10058.001"

# Loop across potentially-relevant years
for(focal_yr in 2013:this_year){
  # focal_yr <- 2013

  # Identify any data for the focal year that we already have
  local_data <- dir(path = file.path("data", "raw_neon"), 
                    pattern = paste0("_", focal_yr, ".csv"))

  # Skip this year if we already have the expected number of files
  ## Note that the expected number of local files is hard-coded
  if(length(local_data) == 3){ 
    
    # Completion message
    message("Data for ", focal_yr, " already downloaded")

  # Otherwise, download the data!
  } else {
    
    # Progress message
    message("Acquiring data for ", focal_yr)

    # Grab start/end months for the focal year
    focal_start <- paste0(focal_yr, "-01")
    if(focal_yr != this_year){
      focal_end <- paste0(focal_yr, "-12")
    } else { focal_end <- this_month }

    # Identify relevant plant data
    focal_list <- neonUtilities::loadByProduct(
      dpID = focal_dpi, 
      startdate = focal_start, 
      enddate = focal_end, 
      check.size = F)  

    # Export relevant stuff
    for(wanted_file in c("div_1m2Data", "div_10m2Data100m2Data", "variables_10058")){
      # wanted_file <- "div_1m2Data"

      # Actually export files (folders created in `00_setup.r`)
      write.csv(x = focal_list[[wanted_file]], na = '', row.names = F,
                file = file.path("data", "raw_neon", 
                                paste0("plant_", wanted_file, "_", focal_yr, ".csv")))

    } # Close export loop
  } # Close download conditional
} # Close year loop



# Identify relevant plant data
## Note that this operation is quite slow (several mins)
# plants <- neonUtilities::loadByProduct(
 #  dpID = "DP1.10058.001", 
  # startdate = "2018-01", 
  # enddate = stringr::str_sub(Sys.Date(), start = 1, end = 7), 
  # check.size = F)




## -------------------------------------- ##

## -------------------------------------- ##

# Identify desired NEON DPI codes
dpi_micro <- "DP1.10081.002"


# End ----
