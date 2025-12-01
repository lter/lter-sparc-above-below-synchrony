# SETUP ####

# packages
library(tidyverse)
library(neonUtilities)

# themes, etc.
theme_set(theme_minimal())

# environment
neon_api_token <- "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJ6YWhuLmVjb2xvZ3lAZ21haWwuY29tIiwic2NvcGUiOiJyYXRlOnB1YmxpYyIsImlzcyI6Imh0dHBzOi8vZGF0YS5uZW9uc2NpZW5jZS5vcmcvIiwiZXhwIjoxOTE5OTk2NTA2LCJpYXQiOjE3NjIzMTY1MDYsImVtYWlsIjoiemFobi5lY29sb2d5QGdtYWlsLmNvbSJ9.NNaWiAljMLVQyGO1CXzCzdW3gz8fGCw6Cww5dEsdR7KDeG419gMBWG5tHcl2sDUBY57J4u-aAWIyU7Rhk-Lyhg"
cores <- parallel::detectCores() - 1
product_microbe <- "DP1.10108.001" # soil microbial community
product_plant <- "DP1.10058.001" # Plant presence and percent cover
product_foliar_traits <- "DP1.10026.001" # Plant foliar traits
product_plant_structure <- "DP1.10098.001" # plot-level for finding aboveground biomass

sites <- c("BONA","CLBJ","CPER","GUAN","HARV","KONZ","NIWO","ONAQ","ORNL","OSBS","PUUM",
           "SCBI","SJER","SRER","TALL","TOOL","UNDE","WOOD","WREF","YELL")
startdate <- "2016-01" # '%Y-%m'

# DOWNLOAD ####

## Download data sheets from NEON ####

# get microbiome data sheets
microbe_list <- 
  loadByProduct(dpID=product_microbe, 
                site=sites, 
                startdate=startdate,
                package = "expanded",
                release = "current",
                include.provisional = TRUE,
                nCores = cores,
                check.size = FALSE,
                token = neon_api_token)
saveRDS(microbe_list,"microbe_list.RDS")

# add dates convenience columns
microbe_list$mmg_soilMarkerGeneSequencing_16S$m_y <- 
  paste0(month(microbe_list$mmg_soilMarkerGeneSequencing_16S$collectDate),
         "-",
         year(microbe_list$mmg_soilMarkerGeneSequencing_16S$collectDate))
microbe_list$mmg_soilMarkerGeneSequencing_ITS$m_y <- 
  paste0(month(microbe_list$mmg_soilMarkerGeneSequencing_ITS$collectDate),
         "-",
         year(microbe_list$mmg_soilMarkerGeneSequencing_ITS$collectDate))
microbe_list$mmg_soilMarkerGeneSequencing_16S$year <- 
  year(microbe_list$mmg_soilMarkerGeneSequencing_16S$collectDate)
microbe_list$mmg_soilMarkerGeneSequencing_ITS$year <- 
  year(microbe_list$mmg_soilMarkerGeneSequencing_ITS$collectDate)

microbe_list$mmg_soilMarkerGeneSequencing_ITS$date <- 
paste(
  microbe_list$mmg_soilMarkerGeneSequencing_ITS$collectDate %>% lubridate::year(),
  microbe_list$mmg_soilMarkerGeneSequencing_ITS$collectDate %>% lubridate::month(),
  microbe_list$mmg_soilMarkerGeneSequencing_ITS$collectDate %>% lubridate::day(),
  sep= "-"
) %>% as.Date()

microbe_date_summary <- 
microbe_list$mmg_soilMarkerGeneSequencing_ITS %>% 
  group_by(siteID,year) %>% 
  reframe(dates=unique(date))

# get plant cover data sheets
plant_list <- 
  loadByProduct(dpID=product_plant, 
                site=sites, 
                startdate=startdate,
                package = "expanded",
                release = "current",
                include.provisional = TRUE,
                nCores = cores,
                check.size = FALSE,
                token = neon_api_token)
saveRDS(plant_list,"plant_list.RDS")

# data frame: div_1m2Data34 has the percentcover info

# add %m-%Y column to div_1m2Data34
plant_list$div_1m2Data$m_y <- 
  paste0(month(plant_list$div_1m2Data$endDate),
         "-", 
         year(plant_list$div_1m2Data$endDate))

# for each site, in each year, find the date range for percentcover
plant_list$div_1m2Data$year <- year(plant_list$div_1m2Data$endDate)
plant_date_summary <- 
  plant_list$div_1m2Data %>% 
  group_by(siteID,year) %>% 
  summarize(begin = min(endDate),
            end = max(endDate))


## Select Microbiome samples ####
# Select one microbiome sampling point to match plant %cover
# find microbiome (M) dates that are within the plant sampling window (Wp)
# if one M inside Wp, select it
# if no M inside Wp, select nearest M before Wp
# if no M inside Wp & no M before Wp, select nearest M after Wp
# if multiple M inside Wp, select earliest



chosen <- microbe_date_summary %>%
  inner_join(plant_date_summary, by = c("siteID","year")) %>%
  group_by(siteID, year, begin, end) %>%
  summarise(
    pick_in     = { x <- dates[dates >= begin & dates <= end]; if (length(x)) min(x) else as.Date(NA) },
    pick_before = { x <- dates[dates <  begin];                  if (length(x)) max(x) else as.Date(NA) },
    pick_after  = { x <- dates[dates >  end];                    if (length(x)) min(x) else as.Date(NA) },
    .groups = "drop"
  ) %>%
  mutate(chosen = coalesce(pick_in, pick_before, pick_after)) %>%
  select(siteID, year, chosen) %>%
  arrange(siteID, year)
chosen$m_y <- 
  paste0(month(chosen$chosen), "-", chosen$year)

# subset microbe data frames to chosen sampling points

microbe_list$mmg_soilRawDataFiles$m_y <- 
  paste0(month(microbe_list$mmg_soilRawDataFiles$collectDate),
         "-",
         year(microbe_list$mmg_soilRawDataFiles$collectDate))


raw_fastqs <- microbe_list$mmg_soilRawDataFiles %>%
  mutate(collectDate = as.Date(collectDate),
         year = year(collectDate)) %>%
  inner_join(
    chosen %>% transmute(siteID, year, chosen_month = month(chosen)),
    by = c("siteID","year")
  ) %>%
  filter(!is.na(rawDataFilePath), month(collectDate) == chosen_month)

microbe_date_summary_chosen <- 
raw_fastqs %>% 
  group_by(siteID,year) %>% 
  reframe(dates = unique(m_y)) %>% 
  mutate(m_y = dates,
         dates = as.Date(paste(
           year,
           dates %>% str_split("-") %>% map_chr(1),
           "28",
           sep="-"
         )))


# add helper columns
raw_fastqs <- 
  raw_fastqs %>% 
  mutate(amplicon = case_when(grepl("16S",rawDataFileName) ~ "Bacteria",
                              grepl("ITS",rawDataFileName) ~ "Fungi",
                              TRUE ~ NA),
         read = case_when(grepl("_R1",rawDataFileName) ~ "Forward",
                          grepl("_R2",rawDataFileName) ~ "Reverse"))



# make raw_fastqs match up with plot key for plant data
# key is "plotID" in a given year
raw_fastqs <- 
  raw_fastqs %>% 
  mutate(plotID = str_split(namedLocation,"\\.") %>% map_chr(1))

# add amplicon/primer info

microbe_list$mmg_soilPcrAmplification_16S$dnaSampleID %in% raw_fastqs$dnaSampleID
microbe_list$mmg_soilPcrAmplification_16S %>% 
  dplyr::select(dnaSampleID,forwardPrimer,reversePrimer,targetGene) %>% 
  left_join(microbe_list$mmg_soilPcrAmplification_ITS %>% 
              dplyr::select(dnaSampleID,forwardPrimer,reversePrimer,targetGene))

bacteria <- 
microbe_list$mmg_soilPcrAmplification_16S %>% 
  dplyr::select(targetGene,dnaSampleID,forwardPrimer,reversePrimer) %>% 
  right_join(raw_fastqs %>% filter(amplicon == "Bacteria")) %>% 
  unique.data.frame()

fungi <- 
  microbe_list$mmg_soilPcrAmplification_ITS %>% 
  dplyr::select(targetGene,dnaSampleID,forwardPrimer,reversePrimer) %>% 
  right_join(raw_fastqs %>% filter(amplicon == "Fungi")) %>% 
  unique.data.frame()

full <- bind_rows(bacteria,fungi) %>% 
  dplyr::select(plotID,year,m_y,dnaSampleID,siteID,collectDate,
                forwardPrimer,reversePrimer,sequencerRunID,
                rawDataFileName,rawDataFilePath,rawDataFileDescription,amplicon,read) %>% 
  unique.data.frame()
saveRDS(full,"./full_microbiome_dataset_info.RDS")

# download each seq run to its own directory
for(run in unique(full$sequencerRunID)){
  # subset to seqrun
  seqrun <- full %>% 
    dplyr::filter(sequencerRunID == run)
  seqrun$amplicon
  # make new directory for that run
  dir.create(file.path("./data/raw",run,"Bacteria"),recursive = TRUE)
  dir.create(file.path("./data/raw",run,"Fungi"),recursive = TRUE)
  
  # download the files to that directory
  download.file(url=seqrun[["rawDataFilePath"]],destfile = file.path("./data/raw",run,seqrun[["amplicon"]],seqrun[["rawDataFileName"]]),method="libcurl")
  
}



