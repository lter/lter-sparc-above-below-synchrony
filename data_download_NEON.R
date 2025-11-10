# -----------------------------------------------------------------------------#
# Data download from NEON using neonUtilities package
# Original Author: L. McKinley Nevins 
# November 3, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     neonUtilities v 3.0.1
#                     ecocomDP v 1.3.2
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(neonUtilities); packageVersion("neonUtilities")
library(ecocomDP); packageVersion("ecocomDP")

#################################################################################
#                               Main Workflow                                   #
#  Use the neonUtilities package to download data from the NEON sites for       #
#  above and belowground community surveys.                                     #
#                                                                               #
#################################################################################

# Following tutorial:
# https://www.neonscience.org/resources/learning-hub/tutorials/neon-biodiversity-ecocomdp-cyverse

########

# load in NEON token 

wd <- "~/Dropbox/WSU/LTER_SPARC/NEON/"
setwd(wd)

source(paste0(wd, "/neon_token_source.R"))


###############
# Pre-evaluation of the NEON data showed that above ground vegetation data is available 
# for every site repeated for many years. 
   # Plant presence and percent cover: DP1.10058.001 (SUCCESS)

# For belowground, there is less coverage, but it does seem like sampling coincided 
# with above ground sampling. There are a few different options that have the coverage 
# of a few years. Could be combined to get better coverage across sites, potentially.
# Not great temporal coverage though, in general all sampling was 2016-2018/19
  # Soil community taxonomy: DP1.10081.002 (SUCCESS) 
      # Soil community composition is the previous version of this that was depricated 

# ***** For all of these, the date of the file save to the project folder can indicate the date accessed for the citation.

###########################################

# The loadByProduct function will extract the zip files, stack the data, and load 
# it into your R environment. See this cheatsheet 
# (https://www.neonscience.org/sites/default/files/cheat-sheet-neonUtilities.pdf)  
# for more information on neonUtilities.


#########################################################################################################

# DOWNLOAD

# Plant data 

plants <- neonUtilities::loadByProduct(
  dpID = "DP1.10058.001", # the NEON plant presence and percent cover data product
  site = c("BONA", "CLBJ", "CPER", "GUAN", "HARV", "KONZ", "NIWO", "ONAQ",
           "ORNL", "OSBS", "PUUM", "SCBI", "SJER", "SRER", "TALL", "TOOL",
           "UNDE", "WOOD", "WREF", "YELL"), # doing just core sites that have compatible data with belowground 
  startdate = "2016-01", # start year-month - no compatible belowground data earlier than this 
  enddate = "2024-12", # end year-month - data release is through December 2023
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size


# Inspect 
names(plants)

# extract items from list and put in R env. 
plants %>% list2env(.GlobalEnv)


# readme has the same information as what you will find on the landing page on the data portal

## about some of the variables: 

# div_1m2Data - Plant species identifications and cover, and cover of abiotic variables within 1 square meter subplots
# div_10m2Data100m2Data - Plant species identifications within 10 square meter and 100 square meter subplots
# variables_10058 - Defines all of the columns in the data files 

## The tutorial has a lot of useful data cleaning and checking steps, but for now I just want to save the raw files. 

# save some data files for reference 

# following team naming conventions for the data file 

# 1m subplot data 
write.csv(div_1m2Data,"~/Dropbox/WSU/LTER_SPARC/NEON/Plants_1m_NEON.csv", row.names = FALSE)

# 10 and 100 m subplot data 
write.csv(div_10m2Data100m2Data,"~/Dropbox/WSU/LTER_SPARC/NEON/Plants_10m_NEON.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10058,"~/Dropbox/WSU/LTER_SPARC/NEON/Plants_NEON_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(plants,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_plant_data.rds")


#########################################################################################################

# Management and events reporting 

manage <- neonUtilities::loadByProduct(
  dpID = "DP1.10111.001", # the NEON management and event reporting product 
  site = c("BONA", "CLBJ", "CPER", "GUAN", "HARV", "KONZ", "NIWO", "ONAQ",
           "ORNL", "OSBS", "PUUM", "SCBI", "SJER", "SRER", "TALL", "TOOL",
           "UNDE", "WOOD", "WREF", "YELL"), # doing just core sites that have compatible data with belowground 
  startdate = "2016-01", # start year-month - no compatible belowground data earlier than this 
  enddate = "2024-12", # end year-month - data release is through December 2024
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size


# Inspect 
names(manage)

# extract items from list and put in R env. 
manage %>% list2env(.GlobalEnv)


# actual event data 
write.csv(sim_eventData,"~/Dropbox/WSU/LTER_SPARC/NEON/site_events_NEON.csv", row.names = FALSE)

# categorical codes
write.csv(categoricalCodes_10111,"~/Dropbox/WSU/LTER_SPARC/NEON/site_event_codes_NEON.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10111,"~/Dropbox/WSU/LTER_SPARC/NEON/site_events_NEON_meta.csv", row.names = FALSE)


#########################################################################################################

#####################
## Exploratory, Geoff is doing the raw sequences 
#####################


# DOWNLOAD

# Download the soil microbial community taxonomy data

microbes_com <- neonUtilities::loadByProduct(
  dpID = "DP1.10081.002", # the NEON soil microbial community taxonomy data product 
  # site = c(""), # ignoring site because I want all sites that have data available
  startdate = "2018-01", # start year-month - Anything before this is lower quality 
  enddate = "2024-01", # end year-month - through the end of the time that the soil data was collected 
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size

# Inspect 
names(microbes_com)

# extract items from list and put in R env. 
microbes_com %>% list2env(.GlobalEnv)

#about some of the files: 
# mct_soilSampleMetadata_16S - Taxon table metadata for soil microbes from sequence variant data analysis of the 16S marker gene
# mct_soilSampleMetadata_ITS - Taxon table metadata for soil microbes from sequence variant data analysis of the ITS region
# variables_10081 - Defines all of the columns in the data files 

# these two metadata files have links to csv files listed for each of the actual datasets. If you plug one of the csv links into a 
# web browser it allows you to download the data. 

# this is still organized as site-months and the files are huge (16S has 8,765 rows; ITS has 11,232 rows) so will need to figure out
# how to extract and process the data from the urls. 


# save some data files for reference 

# NOT following the naming conventions for the data files because these are in an intermediate stage at the moment

# community composition data - 16S
write.csv(mct_soilSampleMetadata_16S,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_16S_COMP.csv", row.names = FALSE)

# community composition data - ITS
write.csv(mct_soilSampleMetadata_ITS,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_ITS_COMP.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10081,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_B_COMP_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(microbes_com,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_microbe_comp_data_2024.rds")


#########################################################################################################

