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
#                               Main workflow                                   #
#  Use the neonUtilities package to download data from the NEON sites for       #
#  above and belowground community surveys.                                     #
#                                                                               #
#################################################################################

# Following tutorial:
# https://www.neonscience.org/resources/learning-hub/tutorials/neon-biodiversity-ecocomdp-cyverse

########

# load in NEON token 

wd <- "~/Dropbox/WSU/SSECR/ssecr-above-belowground-coupling/" # this will depend on your local machine
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
  # Soil microbe group abundances: DP1.10109.001 (SUCCESS)
  # Soil microbe marker gene sequences: DP1.10108.001 (SUCCESS)
  # Soil microbe metagenome sequences: DP1.10107.001 (SUCCESS)

# Biomass datasets: 
  # Root biomass and chemistry, periodic: DP1.10067.001 (SUCCESS)
        # Not doing the megapit dataset because this was collected just once at site establishment 
  # Soil microbe biomass: DP1.10104.001  (SUCCESS)


# ***** For all of these, the date of the file save to the project folder can indicate the date accessed for the citation.

###########################################

# The loadByProduct function will extract the zip files, stack the data, and load 
# it into your R environment. See this cheatsheet 
# (https://www.neonscience.org/sites/default/files/cheat-sheet-neonUtilities.pdf)  
# for more information on neonUtilities.


#########################################################################################################

# DOWNLOAD

# Test download of the plant data 

plants <- neonUtilities::loadByProduct(
  dpID = "DP1.10058.001", # the NEON plant presence and percent cover data product
  site = c("BONA", "CLBJ", "CPER", "GUAN", "HARV", "KONZ", "NIWO", "ONAQ",
           "ORNL", "OSBS", "PUUM", "SCBI", "SJER", "SRER", "TALL", "TOOL",
           "UNDE", "WOOD", "WREF", "YELL"), # doing just core sites that have compatible data with belowground 
  startdate = "2016-01", # start year-month - no compatible belowground data earlier than this 
  enddate = "2023-12", # end year-month - data release is through December 2023
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
write.csv(div_1m2Data,"~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_001.csv", row.names = FALSE)

# 10 and 100 m subplot data 
write.csv(div_10m2Data100m2Data,"~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_002.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10058,"~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_001.2_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(plants,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_plant_data_2022.rds")

test <- readRDS("~/Dropbox/WSU/LTER_SPARC/NEON/NEON_plant_data_2022.rds")


#########################################################################################################

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

######################
# !! These are available if desired, but haven't been downloaded yet !! 
######################


# DOWNLOAD

# Download the soil microbial group abundance data 

microbes_group <- neonUtilities::loadByProduct(
  dpID = "DP1.10109.001", # the NEON soil microbial group abundance data product 
  # site = c(""), # ignoring site because I want all sites that have data available
  startdate = "2012-12", # start year-month - just before the first microbial sampling was started in 2014
  enddate = "2019-12", # end year-month - through the end of the time that the soil data was collected 
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size

# Inspect 
names(microbes_group)

# extract items from list and put in R env. 
microbes_group %>% list2env(.GlobalEnv)

# data file of interest:
# mga_soilGroupAbundances - Laboratory results of gene copy number data in soil samples
# variables_10109 - Defines all of the columns in the data files 

# save some data files for reference 

# following team naming conventions for the data file 

# group abundance data 
write.csv(mga_soilGroupAbundances,"~/Dropbox/WSU/LTER_SPARC/NEON/B_COM_NEON_001.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10109,"~/Dropbox/WSU/LTER_SPARC/NEON/B_COM_NEON_001.2_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(microbes_group,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_microbe_group_data_2024.rds")


#########################################################################################################

# DOWNLOAD

# Download the soil microbial marker gene sequence data 

marker_group <- neonUtilities::loadByProduct(
  dpID = "DP1.10108.001", # the NEON soil microbial marker gene sequence data 
  # site = c(""), # ignoring site because I want all sites that have data available
  package = "expanded", 
  startdate = "2012-12", # start year-month - just before the first microbial sampling was started in 2014
  enddate = "2024-12", # end year-month - marker gene sampling still seems to be ongoing, but provisional data is only through 2019
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size

# Inspect 
names(marker_group)

# extract items from list and put in R env. 
marker_group %>% list2env(.GlobalEnv)

# data file of interest:
# mmmg_soilRawDataFiles - Raw sequence data files for microbial marker gene sequencing from soil samples
# variables_10108 - Defines all of the columns in the data files 

# Defining the package as 'expanded' solved the problems of the data files being absent 


# marker gene raw data files
write.csv(mmg_soilRawDataFiles,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_marker_gene_raw_data.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10108,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_B_marker_gene_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(marker_group,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_marker_gene_data_2024.rds")


#########################################################################################################

# DOWNLOAD

# Download the soil microbial metagenome sequence data 

microbes_metagen <- neonUtilities::loadByProduct(
  dpID = "DP1.10107.001", # the NEON soil microbial metagenome sequence data product 
  # site = c(""), # ignoring site because I want all sites that have data available
  package = "expanded", 
  startdate = "2012-12", # start year-month - just before the first microbial sampling was started in 2014
  enddate = "2024-12", # end year-month - through the end of the time that the soil data was collected 
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size

# Inspect 
names(microbes_metagen)

# extract items from list and put in R env. 
microbes_metagen %>% list2env(.GlobalEnv)

# data file of interest:
# mms_rawDataFiles - Soil metagenomics raw sequence data
# variables_10107 - Defines all of the columns in the data files 

# save some data files for reference 

# following team naming conventions for the data file 

# metagenome data 
write.csv(mms_rawDataFiles,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_metagenome_raw_data.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10107,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_B_metagenome_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(microbes_metagen,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_metagenome_data_2024.rds")

#########################################################################################################

# DOWNLOAD

# Download the Root biomass and chemistry, periodic data

root_biomass <- neonUtilities::loadByProduct(
  dpID = "DP1.10067.001", # the NEON Root biomass and chemistry, periodic data product 
  # site = c(""), # ignoring site because I want all sites that have data available
  package = "expanded", 
  startdate = "2015-12", # start year-month - started in 2016, sites sampled every 5 years 
  enddate = "2024-12", # end year-month - through the end of the time that the soil data was collected 
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size

# Inspect 
names(root_biomass)

# extract items from list and put in R env. 
root_biomass %>% list2env(.GlobalEnv)

# data file of interest:
# bbc_rootChemistry - C and N and stable isotopes for all root samples 
# bbc_rootmass - Root dry mass, also includes size classes and if mycorrhizae were visible, etc. 
# bbc_percore - Info on root volume per sampled core of soil 
# variables_10067 - Defines all of the columns in the data files 

# save some data files for reference 
# root chemistry data 
write.csv(bbc_rootChemistry,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_rootchem_raw_data.csv", row.names = FALSE)

# root mass data 
write.csv(bbc_rootmass,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_rootmass_raw_data.csv", row.names = FALSE)

# root volume per core data 
write.csv(bbc_percore,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_rootvolume_raw_data.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10067,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_B_rootbiomass_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(root_biomass,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_rootbiomass_data_2024.rds")

# Saved the first three, these could be combined into one file, probably, so this important info is more condensed 

#########################################################################################################

# DOWNLOAD

# Download the Soil microbe biomass data

microbe_biomass <- neonUtilities::loadByProduct(
  dpID = "DP1.10104.001", # the NEON Soil microbe biomass data product 
  # site = c(""), # ignoring site because I want all sites that have data available
  package = "expanded", 
  startdate = "2016-12", # start year-month - started March 2017, sampling happens 1-3 times per year 
  enddate = "2024-12", # end year-month - ongoing sampling, ~240 day lag time 
  token = Sys.getenv("NEON_TOKEN"), # use NEON_TOKEN environmental variable
  check.size = F) # proceed with download regardless of file size

# Inspect 
names(microbe_biomass)

# extract items from list and put in R env. 
microbe_biomass %>% list2env(.GlobalEnv)

# data file of interest:
# sme_microbialBiomass - lab results of microbial biomass data  
# sme_scaledMicrobialBiomass - Laboratory results of microbial biomass data, scaled to an internal standard
# sme_batchResults - Microbial biomass batch-level data, summarizes a total lipid volume combined 
# variables_10104 - Defines all of the columns in the data files 

# save some data files for reference 
# microbial biomass data 
write.csv(sme_microbialBiomass,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_microbe_biomass_raw_data.csv", row.names = FALSE)

# microbial biomass scaled data  
write.csv(sme_scaledMicrobialBiomass,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_microbe_biomass_scaled_data.csv", row.names = FALSE)

# microbial biomass batch-level data  
write.csv(sme_batchResults,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_microbe_biomass_batch_data.csv", row.names = FALSE)

# list of variables 
write.csv(variables_10104,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_B_microbebiomass_meta.csv", row.names = FALSE)

# save the whole stack of files as an R object so I could load it back in later if necessary 
saveRDS(microbe_biomass,"~/Dropbox/WSU/LTER_SPARC/NEON/NEON_microbebiomass_data_2024.rds")

# Saved the first three, these could be combined into one file, probably, so this important info is more condensed 


#########################################################################################################

