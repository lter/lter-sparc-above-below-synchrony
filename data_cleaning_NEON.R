# -----------------------------------------------------------------------------#
# Cleaning and organizing downloaded NEON data
# Original Author: L. McKinley Nevins 
# November 3, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     neonUtilities v 3.0.1
#                     ecocomDP v 1.3.1
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(neonUtilities); packageVersion("neonUtilities")
library(ecocomDP); packageVersion("ecocomDP")

#################################################################################
#                               Main workflow                                   #
#  Explore the data structure of the various data products downloaded from      #
#  NEON. Figure out how to unpack all csv's and make nice stacked data files    #
#  that can then be matched up according to site-months to pair above and below.#
#                                                                               #
#################################################################################

########

# load in NEON token 

wd <- "~/Dropbox/WSU/SSECR/ssecr-above-belowground-coupling/" # this will depend on your local machine
setwd(wd)

source(paste0(wd, "/neon_token_source.R"))

#################################################################################
# PLANT COMMUNITY DATA

# load in the two data files, one at 1m subplots, and one at 10m subplots 

# 1m subplot data - percent cover data 
plant_1m <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_001.csv")

# 10 and 100 m subplot data - presence/absence data 
plant_10m <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_002.csv")

#check variables list 
plant_meta <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_001.2_meta.csv")

# potential variables of interest to retain columns: 
  # - domainID, siteID, decimalLatitude, decimalLongitude, elevation, nlcdClass, plotID, subplotID, endDate, divDataType
  # taxonID, scientificName, taxonRank, family, nativeStatusCode, otherVariables, percentCover, heightPlantSpecies

  # divDataType specifies whether the row is actual plant abundance data, or another type that can be filtered out 
    # plant data will have a 'PlantSpecies' classification 

  # otherVariables lists other things measured in the percent cover measurement, like soil, lichen, etc. 

###########subset 1m subplot data
plant_1m_sub <- select(plant_1m, domainID, siteID, decimalLatitude, decimalLongitude, elevation, nlcdClass, plotID, subplotID, endDate, divDataType,
                       taxonID, scientificName, taxonRank, family, nativeStatusCode, otherVariables, percentCover, heightPlantSpecies)

# narrowed down to 18 more relevant variables 

###########tangent
# pull out a little dataset of the categorical site specifiers to use in other places 
site_info <- select(plant_1m_sub, domainID, siteID, nlcdClass)

# change some data formats 
site_info$domainID <- as.factor(site_info$domainID)
site_info$siteID <- as.factor(site_info$siteID) #20 core terrestrial field sites 
site_info$nlcdClass <- as.factor(site_info$nlcdClass)

#keep only unique, so should just end up with one row for each NEON site 
str(site_info$nlcdClass)

site_info <- site_info[!duplicated(site_info[1:3]),]

# there are many sites that have multiple ecosystem types 

# save file 

write.csv(site_info, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_site_info.csv")
##############end of tangent

# convert some variables to different formats 

# domainID
plant_1m_sub$domainID <- as.factor(plant_1m_sub$domainID) #20 levels 

# siteID
plant_1m_sub$siteID <- as.factor(plant_1m_sub$siteID) #20 levels 

# nlcdClass (like ecosystem type)
plant_1m_sub$nlcdClass <- as.factor(plant_1m_sub$nlcdClass) #11 levels, can use to filter out crop fields, etc. 

# divDataType 
plant_1m_sub$divDataType <- as.factor(plant_1m_sub$divDataType) #2 levels, either otherVariables or PlantSpecies

# otherVariables 
plant_1m_sub$otherVariables <- as.factor(plant_1m_sub$otherVariables) #19 levels, coverage types for percent cover 


# save cleaned 1m subplot data 
write.csv(plant_1m_sub, "~/Dropbox/WSU/LTER_SPARC/NEON/clean_A_COM_NEON_001.csv")

#############subset 10m subplot data
plant_10m_sub <- select(plant_10m, domainID, siteID, decimalLatitude, decimalLongitude, elevation, nlcdClass, plotID, subplotID, endDate,
                       taxonID, scientificName, taxonRank, family, nativeStatusCode)

# narrowed down to 14 more relevant variables 

# convert some variables to different formats 

# domainID
plant_10m_sub$domainID <- as.factor(plant_10m_sub$domainID) #20 levels 

# siteID
plant_10m_sub$siteID <- as.factor(plant_10m_sub$siteID) #47 levels 

# nlcdClass (like ecosystem type)
plant_10m_sub$nlcdClass <- as.factor(plant_10m_sub$nlcdClass) #11 levels, can use to filter out crop fields, etc. 


# save cleaned 1m subplot data 
write.csv(plant_10m_sub, "~/Dropbox/WSU/LTER_SPARC/NEON/clean_A_COM_NEON_002.csv")

#################################################################################
# SOIL COMMUNITY DATA 

# load in two data files 
# these are stacked csv's of taxon tables for bacterial and fungal communities 

# community composition data - 16S
bacteria <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/NEON_16S_COMP.csv")

# community composition data - ITS
fungi <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/NEON_ITS_COMP.csv")

# list of variables 
below_meta <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/NEON_B_COMP_meta.csv")

####bacteria first 

# subsetting variables of interest - domainID, siteID, plotID, collectDate, analysisDate, 
# downloadFileUrl, downloadFileName, sequenceCountQF

# 8 variables of interest 

  # sequenceCountQF specifies if the sequence reads are of sufficient quality. Some are noted as 
  # below threshold and should be excluded 

bacteria_sub <- select(bacteria, domainID, siteID, plotID, collectDate, analysisDate, 
                       downloadFileUrl, downloadFileName, sequenceCountQF)

bacteria_sub$sequenceCountQF <- as.factor(bacteria_sub$sequenceCountQF)

str(bacteria_sub)

# want to read through the 'sequenceCountQF' column and only pull csv's that have 
# sequence counts that passed the threshold as 'OK' 
bacteria_sub <- bacteria_sub[bacteria_sub$sequenceCountQF == 'OK', ]

#drops down to 8,652 files instead 

# each row here is a separate csv file. These will need to be looped through to be unpacked 

# LOOP
# will eventually up the number range to all rows in the bacteria_sub dataframe to cover each of 
# the csv files. This will be 8,652 high quality files, so could take a very long time 
bacteria_list<- list()
for (focal_file in bacteria_sub$downloadFileUrl[1:8652]){
  
  #print(focal_file)

  files <- read.csv(focal_file)
  
  #this lists the separate dataframes out of the csv files 
  bacteria_list[[focal_file]] <- data.frame(files)
  
  #this combines the dataframes into one bigger one. They have the same columns 
  #so this isn't problematic 
  bacteria_df <- purrr::list_rbind(x = bacteria_list)
  
}

bacteria_list

bacteria_df

# doing rest of cleaning operations on a different version
bacteria_df2 <- bacteria_df %>% 
  tidyr::separate_wider_delim(cols = dnaSampleID, 
                              delim = "_",
                              names = c("siteID", "Other"),
                              cols_remove = FALSE)

#this is a dumb way to do it but it works 
bacteria_df2 <- bacteria_df2 %>% 
  tidyr::separate_wider_delim(cols = dnaSampleID, 
                              delim = "-",
                              names = c("plotID", "NA", "NA", "NA", "date", "NA", "NA"),
                              cols_remove = FALSE)


# Check that out
str(bacteria_df2)

#select focal columns 
bacteria_df2 <- select(bacteria_df2, siteID, plotID, date, dnaSampleID, sequenceName, taxonSequence, 
                      completeTaxonomy, domain, kingdom, phylum, class, order, family, genus, specificEpithet, 
                      scientificName, individualCount)
                      
#now need to split the date column into different components 

# formatting date column to have components separated 
bacteria_df2 <- bacteria_df2 %>% mutate(date = as.Date(date, format = "%Y%m%d")) 

# moving hyphenated components into new, separate columns 
bacteria_df2 <- bacteria_df2 %>% 
  tidyr::separate_wider_delim(cols = date, 
                              delim = "-",
                              names = c("year", "month", "day"),
                              cols_remove = TRUE)


# Now would like to add columns that specify the domain. The ecosystem type is a bit more complicated 
# because multiple are present and sampled in each site (at least for the vegetation). Not doing this 
# for now, but can do later once I figure out where the soils are specifically sampled from. 

# add domainID as a column from the site_info file into the bacteria_df file, using siteID to merge 

#remove nlcdClass column from site_info for now 

site_info <- select(site_info, domainID, siteID)
# remove duplicates again 
site_info <- site_info[!duplicated(site_info[1:2]),]

# save bit of site info to use 
write.csv(site_info, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_site_info_minimum.csv")

bacteria_df2 <- merge(bacteria_df2, site_info, by = "siteID", all.x = TRUE)

# this is looking good! I think I can save this cleaned version for now 

# save cleaned bacterial community data 
write.csv(bacteria_df2, "~/Dropbox/WSU/LTER_SPARC/NEON/clean_NEON_16S_COMP.csv")

########fungi 

# subsetting variables of interest - domainID, siteID, plotID, collectDate, analysisDate, geneticSampleID, 
# downloadFileUrl, downloadFileName, sequenceCountQF

# 9 variables of interest 

# sequenceCountQF specifies if the sequence reads are of sufficient quality. Some are noted as 
# below threshold and should be excluded 

fungi_sub <- select(fungi, domainID, siteID, plotID, collectDate, analysisDate, geneticSampleID, 
                     downloadFileUrl, downloadFileName, sequenceCountQF)

#starting with 11,232 csv files 

fungi_sub$sequenceCountQF <- as.factor(fungi_sub$sequenceCountQF)

str(fungi_sub)

# want to read through the 'sequenceCountQF' column and only pull csv's that have 
# sequence counts that passed the threshold as 'OK' (4,000 reads)
fungi_sub <- fungi_sub[fungi_sub$sequenceCountQF == 'OK', ]

# drops down to 3,991 files instead - major bummer! 
# seems like a huge waste of resources to have so much that isn't usable? 

# each row here is a separate csv file. These will need to be looped through to be unpacked 

# LOOP
fungi_list<- list()
for (focal_file in fungi_sub$downloadFileUrl[1:3991]){
  
  #print(focal_file)
  
  files <- read.csv(focal_file)
  
  #this lists the separate dataframes out of the csv files 
  fungi_list[[focal_file]] <- data.frame(files)
  
  #this combines the dataframes into one bigger one. They have the same columns 
  #so this isn't problematic 
  fungi_df <- purrr::list_rbind(x = fungi_list)
  
}

fungi_list

fungi_df

# doing rest of cleaning operations on a different version
fungi_df2 <- fungi_df %>% 
  tidyr::separate_wider_delim(cols = dnaSampleID, 
                              delim = "_",
                              names = c("siteID", "Other"),
                              cols_remove = FALSE)

#this is a dumb way to do it but it works 
fungi_df2 <- fungi_df2 %>% 
  tidyr::separate_wider_delim(cols = dnaSampleID, 
                              delim = "-",
                              names = c("plotID", "NA", "NA", "NA", "date", "NA", "NA"),
                              cols_remove = FALSE)


# Check that out
str(fungi_df2)

#select focal columns 
fungi_df2 <- select(fungi_df2, siteID, plotID, date, dnaSampleID, sequenceName, taxonSequence, 
                       completeTaxonomy, domain, kingdom, phylum, class, order, family, genus, specificEpithet, 
                       scientificName, individualCount)

#now need to split the date column into different components 

# formatting date column to have components separated 
fungi_df2 <- fungi_df2 %>% mutate(date = as.Date(date, format = "%Y%m%d")) 

# moving hyphenated components into new, separate columns 
fungi_df2 <- fungi_df2 %>% 
  tidyr::separate_wider_delim(cols = date, 
                              delim = "-",
                              names = c("year", "month", "day"),
                              cols_remove = TRUE)

# add domainID as a column from the site_info file into the fungi_df2 file, using siteID to merge 

# save bit of site info to use 
site_info <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/NEON_site_info_minimum.csv", row.names = 1)

fungi_df2 <- merge(fungi_df2, site_info, by = "siteID", all.x = TRUE)

# save cleaned fungal community data 
write.csv(fungi_df2, "~/Dropbox/WSU/LTER_SPARC/NEON/clean_NEON_ITS_COMP.csv")








