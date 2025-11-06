# SETUP ####

# packages
library(tidyverse)
library(neonUtilities)

# themes, etc.
theme_set(theme_minimal())

# environment
readRenviron("~/.Renviron.env") # load secret token
neon_api_token <- Sys.getenv("NEON_TOKEN")
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
microbe_list <- readRDS("~/Desktop/microbe_list.RDS")

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
plant_list <- readRDS("~/Desktop/plant_list.RDS")


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

plant_cov_plot <- 
  plant_date_summary %>% 
  ggplot() +
  geom_segment(linewidth = 3,aes(x=begin,xend=end,y=siteID,yend=siteID)) +
  geom_point(
    data=microbe_date_summary,
    aes(x=dates,y=siteID),
    color="red"
  ) +
  labs(title = "Co-temporal sampling overview",x="Date",
       subtitle = "Bars = Plant %cover sampling coverage\nDots = Microbiome data points\nSite level")
plant_cov_plot
ggsave("~/Desktop/plant_sampling_coverage.png",height = 6,width = 10)

## Select Microbiome samples ####
# Select one microbiome sampling point to match plant %cover
# find microbiome (M) dates that are within the plant sampling window (Wp)
# if one M inside Wp, select it
# if no M inside Wp, select nearest M before Wp
# if no M inside Wp & no M before Wp, select nearest M after Wp
# if multiple M inside Wp, select earliest


plant_date_summary
microbe_date_summary

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

plant_date_summary %>% 
  ggplot() +
  geom_segment(linewidth = 5,aes(x=begin,xend=end,y=siteID,yend=siteID)) +
  geom_point(
    data=microbe_date_summary,
    aes(x=dates,y=siteID),
    color="red",shape=4
  ) +
  geom_point(
    data=microbe_date_summary_chosen,
    aes(x=dates,y=siteID),
    color="green",size=3
  ) +
  labs(title = "Co-temporal sampling overview",x="Date",
       subtitle = "Bars = Plant %cover sampling coverage\nDots = Microbiome data points\nSite level")
ggsave("~/Desktop/selected_microbe_dates.png",height = 6,width = 10)

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
microbe_list$mmg_soilMarkerGeneSequencing_16S$sequencerRunID
microbe_list$mmg_soilMarkerGeneSequencing_ITS$sequencerRunID
microbe_list$mmg_soilPcrAmplification_16S$forwardPrimer
microbe_list$mmg_soilPcrAmplification_16S$reversePrimer
microbe_list$mmg_soilPcrAmplification_ITS$forwardPrimer
microbe_list$mmg_soilPcrAmplification_ITS$reversePrimer

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

for(run in unique(full$sequencerRunID)){
  # subset to seqrun
  seqrun <- full %>% 
    dplyr::filter(sequencerRunID == run)
  seqrun$amplicon
  # make new directory for that run
  dir.create(file.path("./data/raw",run,"Bacteria"),recursive = TRUE)
  dir.create(file.path("./data/raw",run,"Fungi"),recursive = TRUE)
  
  # download the files to that directory
  download.file(url=seqrun[["rawDataFilePath"]],
                destfile = file.path("./data/raw",run,seqrun[["amplicon"]],seqrun[["rawDataFileName"]]),
                method= "libcurl")
  
}



# get a glimpse of sample replication for each plotID
raw_fastqs %>% 
  group_by(plotID,year,amplicon,read) %>% 
  summarize(N=n())
raw_fastqs$plotID %>% table

# number of plots sampled in a given site each year...
table(raw_fastqs$year,raw_fastqs$plotID) %>% 
  as.data.frame() %>% 
  mutate(site = Var2 %>% str_split("_") %>% map_chr(1)) %>%
  dplyr::filter(Freq > 0) %>% 
  group_by(Var2,Var1) %>% 
  summarise(N=n()) %>% 
  mutate(site = Var2 %>% str_split("_") %>% map_chr(1)) %>% 
  ggplot(aes(x=Var1,y=N)) +
  geom_col() +
  facet_wrap(~site) +
  theme(axis.text.x = element_text(angle=90))



# how many files are present in each site/year combo?
raw_fastqs %>% 
  group_by(siteID,namedLocation,year,amplicon,read) %>% 
  summarize(N=n()) %>% 
  ggplot(aes(x=N)) +
  geom_histogram() +
  facet_wrap(~siteID)

# which primer sets were used for which years
table(as.Date(microbe_list$mmg_soilPcrAmplification_16S$collectDate),microbe_list$mmg_soilPcrAmplification_16S$targetSubfragment) %>% 
  as.data.frame() %>% 
  mutate(date = row.names(.)) %>% 
    mutate(year=year(Var1)) %>% 
  ggplot(aes(x=year,y=Freq,fill=Var2)) +
    geom_col() + ggtitle("Bacteria")

table(as.Date(microbe_list$mmg_soilPcrAmplification_ITS$collectDate),microbe_list$mmg_soilPcrAmplification_ITS$targetSubfragment) %>% 
  as.data.frame() %>% 
  mutate(date = row.names(.)) %>% 
  mutate(year=year(Var1)) %>% 
  ggplot(aes(x=year,y=Freq,fill=Var2)) +
  geom_col() + ggtitle("Fungi")






raw_fastqs %>% 
  group_by(siteID,namedLocation,year,amplicon,read) %>% 
  summarize(N=n())

raw_fastqs %>% 
  filter(siteID=="BONA" & year == 2018 & amplicon == "Bacteria" & namedLocation == "BONA_071.basePlot.bgc") %>% 
  View()
raw_fastqs %>% 
  filter(siteID=="BONA" & year == 2019 & amplicon == "Fungi" & namedLocation == "BONA_001.basePlot.bgc") %>% 
  View()



## Download raw fastq files for selected microbiomes ####


dat_list

# get list of sampling times from plants
# use to select only appropriate microbe data timepoints
dat_list$mmg_soilRawDataFiles
# download raw data
microbe_list
download.file()

