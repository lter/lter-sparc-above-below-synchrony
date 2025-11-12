# -----------------------------------------------------------------------------#
# Filter and investigate management and events NEON data
# Original Author: L. McKinley Nevins 
# November 7, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     ecocomDP v 1.3.2
#                     lubridate v 1.9.4
#                     stringr v 1.5.1
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(ecocomDP); packageVersion("ecocomDP")
library(lubridate); packageVersion("lubridate")
library(stringr); packageVersion("stringr")

#################################################################################
#                               Main Workflow                                   #
#  Inspect the management and events data from NEON and match to the plots at   #
#  focus for the paired above-belowground data.                                 #
#                                                                               #
#################################################################################

###################
## 1. READ IN DATA 
###################

# Load in events data 
events <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/site_events_NEON.csv")

# Load in list of microbial plots 
microbe_plots <- read.csv("~/Dropbox/WSU/LTER_SPARC/NEON/microbe_plots_by_year.csv")

## Set year column as numeric 
microbe_plots <- microbe_plots %>%
  mutate(year = as.numeric(year))


######################
## 2. MANIPULATE DATA 
######################

# Summarize events data broadly by sites and types of events 
site_summary <- events %>%
  group_by(siteID, eventType) %>%
  summarize(
    n_events = n(),
    first_event = min(ymd(startDate), na.rm = TRUE),
    last_event = max(ymd(endDate), na.rm = TRUE)) %>%
  arrange(siteID, eventType)


# Trying to figure out the spatial scale of the disturbance and events 
# The locationID code is the location where the sample as collected, and it 
# specifies either the base plots, the mammal plot, etc. We are only interested 
# in the things that would be affecting the base plots 
location_summary <- events %>%
  filter(!is.na(locationID)) %>%
  group_by(siteID, locationID, eventType) %>%
  summarize(
    n_events = n(),
    first_event = min(ymd(startDate), na.rm = TRUE),
    last_event = max(ymd(endDate), na.rm = TRUE))

# Separate out the individual plots affected for matching 
events_long <- events %>%
  separate_rows(locationID, sep = ",\\s*") 

# Filter to just rows for locationID that are basePlots 
events_base <- events_long %>%
  filter(str_detect(locationID, "\\.basePlot.all")) #8,702 events at individual plot level 

# Rename and reformat plotID
events_base <- events_base %>%
  mutate(plotID = str_remove(locationID, "\\.basePlot.all"))

# Get the microbial data and pull out just plotID's to merge with 
m_plots <- microbe_plots %>% select(plotID)

m_plots$plotID <- as.factor(m_plots$plotID)

# Merge by plotID to get just events for matching plots 
events_matched <- m_plots %>%
  merge(events_base, by = "plotID") # plotID factor here as 217 levels, so only 217 plots 
                                    # have something recorded for them? 

# Event start and end dates are in year-month-day format, and we want to think about when 
# the disturbance occurred in relation to the microbial sampling 
# We are also interested in if the events happened between time points, so this is less about 
# getting events that exactly overlap with microbial years, and instead trying to get an understanding 
# of the types of management and events that were occurring to the plots during the whole 
# time period from 2016 - 2024

# methodTypeChoice specifies the type of management and disturbance 

events_matched$methodTypeChoice <- as.factor(events_matched$methodTypeChoice) #21 types 

# want to visualize for each affected plot the type of event and the timing of the event 




