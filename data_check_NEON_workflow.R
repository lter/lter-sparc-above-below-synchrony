# -----------------------------------------------------------------------------#
# Check and clean downloaded NEON data
# Original Author: L. McKinley Nevins 
# February 5, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     neonUtilities v 2.4.2
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
#  Follow the recommended workflow from NEON to check the data. This is         #
#  organizing taxa into a taxon table, gathering site and observation data, and #
#  making final species x site matrices needed for community analyses in vegan. #
#                                                                               #
#################################################################################

# Following the workflow described here for processing NEON biodiversity data: 
# https://www.neonscience.org/resources/learning-hub/tutorials/neon-biodiversity-ecocomdp-cyverse

## Not doing any other prior cleaning, just following the recommended steps: 

## notes are a mix of my own and existing from the lesson materials 
## some sections have been removed if not relevant 

########
# load in NEON token 
wd <- "~/Dropbox/WSU/SSECR/ssecr-above-belowground-coupling/" # this will depend on your local machine
setwd(wd)

source(paste0(wd, "/neon_token_source.R"))

########

# uses outputs generated from the data_cleaning_NEON.R script

#################################################################################
# PLANT COMMUNITY DATA

# load in the two data files, one at 1m subplots, and one at 10m subplots 

# 1m subplot data - percent cover data 
plant_1m <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_001.csv")

# 10 and 100 m subplot data - presence/absence data 
plant_10m <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_002.csv")

#check variables list 
plant_meta <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/A_COM_NEON_001.2_meta.csv")

# only doing this with the plant 1m percent cover data, not the 10m presence absence data
# could come back and process that later 


###### PLANT 1 M DATA ########

# expand dates, into new columns
plant_1m <- plant_1m %>% 
  tidyr::separate_wider_delim(cols = endDate, 
                              delim = "-",
                              names = c("Year", "Month", "Day"),
                              cols_remove = FALSE)


# extract location data into a separate table
table_location_plant_1m <- plant_1m %>%
  # keep only the columns listed below
  select(siteID, 
         domainID,
         namedLocation, 
         decimalLatitude, 
         decimalLongitude, 
         elevation) %>%
  # keep rows with unique combinations of values, i.e., no duplicate records
  distinct()


# create a taxon table, which describes each taxonID that appears in the data set
table_taxon_plant_1m <- plant_1m %>%
  # keep only the columns listed below
  select(taxonID, taxonRank, scientificName,
         family, identificationQualifier,
         identificationReferences) %>%
  # remove rows with duplicate information
  distinct()


# taxon table information for all taxa in 
# our database can be downloaded here:
# takes 1-2 minutes
 full_taxon_table_from_api <- 
   neonOS::getTaxonList(
     "PLANT", 
     token = Sys.getenv("NEON_TOKEN"))
# downloaded this but then didn't actually do anything with it 
 # can save down below for future use 

# Make the observation table.
# check for repeated taxa within an eventID that need to be added together
taxonomyProcessed_summed_plant_1m <- plant_1m %>% 
  select(eventID,
         taxonID,
         percentCover,
         heightPlantSpecies) %>%
  group_by(eventID, taxonID) %>%
  summarize(
    across(c(percentCover, heightPlantSpecies), ~sum(.x, na.rm = TRUE)))


# join summed taxon counts back with sample and field data
table_observation_plant_1m <- taxonomyProcessed_summed_plant_1m %>%
  # Join relevant sample info back in by eventID
  left_join(plant_1m %>% 
              select(eventID,
                     domainID,
                     siteID,
                     namedLocation,
                     Year, Month, Day,
                     endDate,
                     taxonID,
                     family, 
                     scientificName,
                     taxonRank) %>%
              distinct())


# check for duplicate records, should return a table with 0 rows
table_observation_plant_1m %>% 
  group_by(eventID, taxonID) %>% 
  summarize(n_obs = length(eventID)) %>%
  filter(n_obs > 1)

  # this has a lot of "duplicates" because there were multiple plots sampled for each 
  # "event", so a plant can appear multiple times if it was recorded in different plots. 
  # the plot is provided int he namedLocation, so I think this is okay. They could
  # be summed by species further in the future if we desired. 


## Create a site x species table 

# select only site by species density info and remove duplicate records
table_sample_by_taxon_density_long_plant_1m <- table_observation_plant_1m %>%
  select(eventID, taxonID, percentCover) %>%
  distinct() %>%
  # filter out NA's in any of the columns 
  # pivoting wider down below can't be done if there are missing cases in the 
  # eventID column 
  filter(!is.na(percentCover)) %>%
  filter(!is.na(eventID)) %>%
  filter(!is.na(taxonID))

# this results in 50,238 rows, where each row is a complete observation of 
# a taxon with percent cover for a given sampling event at a site 


# pivot to wide format, sum multiple counts per eventID
table_sample_by_taxon_density_wide_plant_1m <- table_sample_by_taxon_density_long_plant_1m %>%
  tidyr::pivot_wider(id_cols = eventID, 
                     names_from = taxonID,
                     values_from = percentCover,
                     values_fill = list(percentCover = 0),
                     values_fn = list(percentCover = sum)) %>%
  column_to_rownames(var = "eventID") 


# check col and row sums -- mins should all be > 0
colSums(table_sample_by_taxon_density_wide_plant_1m) %>% min()
## [1] 0

rowSums(table_sample_by_taxon_density_wide_plant_1m) %>% min()
## [1] 0.5

# This results in a species x site matrix for all of the NEON sites. 
# the rows are the event ID, which comprise the site and year 

######################################################################
## Key components to save out of that process would be:
  # The wide species x site matrix for the plant communities 
  # "table_sample_by_taxon_density_wide_plant_1m"

write.csv(table_sample_by_taxon_density_wide_plant_1m, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_A_COM_species_x_site.csv")

  # The table of observation data for each eventID and taxon 
  # "table_observation_plant_1m"

write.csv(table_observation_plant_1m, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_A_COM_observation_details.csv")

  # The table of all of the location information for each subplot 
  # sampled at each of the NEON sites 
  # "table_location_plant_1m"

write.csv(table_location_plant_1m, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_A_COM_location_details.csv")

  # The full plant taxon table from API
  # "full_taxon_table_from_api"

write.csv(full_taxon_table_from_api, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_plant_taxonomy_reference.csv")
  
#############################################################################################
## FUNGAL COMMUNITY DATA ##

fungi <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/clean_NEON_ITS_COMP.csv")

# extract location data into a separate table
table_location_fungi <- fungi %>%
  # keep only the columns listed below
  select(siteID, 
         plotID,
         domainID) %>%
  # keep rows with unique combinations of values, i.e., no duplicate records
  distinct()


# create a taxon table, which describes each completeTaxonomy that appears in the data set
  # making the choice to use this because there doesn't seem to be an equivalent to the 
  # taxonID column for the plant community. There are no species codes for fungi, and the 
  # other columns don't fully encapsulate the taxonomic classification
table_taxon_fungi <- fungi %>%
  # keep only the columns listed below
  select(taxonSequence, completeTaxonomy, scientificName,
         phylum, class, order, family, genus,
         scientificName) %>%
  # remove rows with duplicate information
  distinct()

# Make the observation table.
# check for repeated taxa within an dnaSampleID that need to be added together
taxonomyProcessed_summed_fungi <- fungi %>% 
  select(dnaSampleID,
         completeTaxonomy,
         individualCount) %>%
  group_by(dnaSampleID, completeTaxonomy) %>%
  summarize(
    across(c(individualCount), ~sum(.x, na.rm = TRUE)))


# join summed taxon counts back with sample and field data
table_observation_fungi <- taxonomyProcessed_summed_fungi %>%
  # Join relevant sample info back in by dnaSampleID
  left_join(fungi %>% 
              select(dnaSampleID,
                     domainID,
                     siteID,
                     plotID,
                     year, month, day) %>%
              distinct())


# check for duplicate records, should return a table with 0 rows
table_observation_fungi %>% 
  group_by(dnaSampleID, completeTaxonomy) %>% 
  summarize(n_obs = length(dnaSampleID)) %>%
  filter(n_obs > 1)

# no duplicates. This is organized where each unique taxonomic assignment in each 
# unique site-month of sampling has it's own row, so there shouldn't be any duplicates

# 275,827 unique observations 

## Create a site x species table 

# select only site by species density info and remove duplicate records
table_sample_by_taxon_density_long_fungi <- table_observation_fungi %>%
  select(dnaSampleID, completeTaxonomy, individualCount) %>%
  distinct() %>%
  # filter out NA's in any of the columns 
  # pivoting wider down below can't be done if there are missing cases in the 
  # eventID column 
  filter(!is.na(dnaSampleID)) %>%
  filter(!is.na(completeTaxonomy)) %>%
  filter(!is.na(individualCount))

# this results in 275,827 rows, where each row is a complete observation of 
# a taxon with abundance count for a given sampling event at a site 


# pivot to wide format, sum multiple counts per dnaSampleID
table_sample_by_taxon_density_wide_fungi <- table_sample_by_taxon_density_long_fungi %>%
  tidyr::pivot_wider(id_cols = dnaSampleID, 
                     names_from = completeTaxonomy,
                     values_from = individualCount,
                     values_fill = list(individualCount = 0),
                     values_fn = list(individualCount = sum)) %>%
  column_to_rownames(var = "dnaSampleID") 


# check col and row sums -- mins should all be > 0
colSums(table_sample_by_taxon_density_wide_fungi) %>% min()
## [1] 2

rowSums(table_sample_by_taxon_density_wide_fungi) %>% min()
## [1] 4,000

# This results in a species x site matrix for all of the NEON sites. 
# the rows are the dnaSampleID, which contains site and year information

######################################################################
## Key components to save out of that process would be:
# The wide species x site matrix for the fungal communities 
# "table_sample_by_taxon_density_wide_fungi"

write.csv(table_sample_by_taxon_density_wide_fungi, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_ITS_COMP_species_x_site.csv")

# The table of observation data for each dnaSampleID and completeTaxonomy
# "table_observation_fungi"

write.csv(table_observation_fungi, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_ITS_COMP_observation_details.csv")

# The table of all of the location information for each subplot 
# sampled at each of the NEON sites 
# "table_location_fungi"

write.csv(table_location_fungi, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_ITS_COMP_location_details.csv")


#############################################################################################
## BACTERIAL COMMUNITY DATA ##

bacteria <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/clean_NEON_16S_COMP.csv")

# extract location data into a separate table
table_location_bac <- bacteria %>%
  # keep only the columns listed below
  select(siteID, 
         plotID,
         domainID) %>%
  # keep rows with unique combinations of values, i.e., no duplicate records
  distinct()


# create a taxon table, which describes each completeTaxonomy that appears in the data set
# making the choice to use this because there doesn't seem to be an equivalent to the 
# taxonID column for the plant community. There are no species codes for bacteria, and the 
# other columns don't fully encapsulate the taxonomic classification
table_taxon_bac <- bacteria %>%
  # keep only the columns listed below
  select(taxonSequence, completeTaxonomy,
         phylum, class, order, family, genus,
         scientificName) %>%
  # remove rows with duplicate information
  distinct()

# Make the observation table.
# check for repeated taxa within an dnaSampleID that need to be added together
taxonomyProcessed_summed_bac <- bacteria %>% 
  select(dnaSampleID,
         completeTaxonomy,
         individualCount) %>%
  group_by(dnaSampleID, completeTaxonomy) %>%
  summarize(
    across(c(individualCount), ~sum(.x, na.rm = TRUE)))


# join summed taxon counts back with sample and field data
table_observation_bac <- taxonomyProcessed_summed_bac %>%
  # Join relevant sample info back in by dnaSampleID
  left_join(bacteria %>% 
              select(dnaSampleID,
                     domainID,
                     siteID,
                     plotID,
                     year, month, day) %>%
              distinct())


# check for duplicate records, should return a table with 0 rows
table_observation_bac %>% 
  group_by(dnaSampleID, completeTaxonomy) %>% 
  summarize(n_obs = length(dnaSampleID)) %>%
  filter(n_obs > 1)

# no duplicates. This is organized where each unique taxonomic assignment in each 
# unique site-month of sampling has it's own row, so there shouldn't be any duplicates

# 1,394,677 unique observations - that's a lotta bacteria! 

## Create a site x species table 

# select only site by species density info and remove duplicate records
table_sample_by_taxon_density_long_bac <- table_observation_bac %>%
  select(dnaSampleID, completeTaxonomy, individualCount) %>%
  distinct() %>%
  # filter out NA's in any of the columns 
  # pivoting wider down below can't be done if there are missing cases in the 
  # eventID column 
  filter(!is.na(dnaSampleID)) %>%
  filter(!is.na(completeTaxonomy)) %>%
  filter(!is.na(individualCount))

# this results in 1,394,677 rows, where each row is a complete observation of 
# a taxon with abundance count for a given sampling event at a site 


# pivot to wide format, sum multiple counts per dnaSampleID
table_sample_by_taxon_density_wide_bac <- table_sample_by_taxon_density_long_bac %>%
  tidyr::pivot_wider(id_cols = dnaSampleID, 
                     names_from = completeTaxonomy,
                     values_from = individualCount,
                     values_fill = list(individualCount = 0),
                     values_fn = list(individualCount = sum)) %>%
  column_to_rownames(var = "dnaSampleID") 


# check col and row sums -- mins should all be > 0
colSums(table_sample_by_taxon_density_wide_bac) %>% min()
## [1] 1

rowSums(table_sample_by_taxon_density_wide_bac) %>% min()
## [1] 4,002

# This results in a species x site matrix for all of the NEON sites. 
# the rows are the dnaSampleID, which contains site and year information

######################################################################
## Key components to save out of that process would be:
# The wide species x site matrix for the bacterial communities 
# "table_sample_by_taxon_density_wide_bac"

write.csv(table_sample_by_taxon_density_wide_bac, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_16S_COMP_species_x_site.csv")

# The table of observation data for each dnaSampleID and completeTaxonomy
# "table_observation_bac"

write.csv(table_observation_bac, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_16S_COMP_observation_details.csv")

# The table of all of the location information for each subplot 
# sampled at each of the NEON sites 
# "table_location_bac"

write.csv(table_location_bac, "~/Dropbox/WSU/LTER_SPARC/NEON/NEON_16S_COMP_location_details.csv")


