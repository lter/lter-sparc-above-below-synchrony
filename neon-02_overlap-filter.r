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
"uid"             "namedLocation"   "domainID"        "siteID"          "plotID"

## -------------------------------------- ##
# Identify Overlapped Years ----
## -------------------------------------- ##

# List files of both types






# End ----
