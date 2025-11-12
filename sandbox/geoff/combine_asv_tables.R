# SETUP ####

## Packages ####
library(dada2)
library(phyloseq)
library(tidyverse)

## Functions ####
as_otu <- function(x, taxa_are_rows = FALSE) {
  m <- as.matrix(x)
  storage.mode(m) <- "integer"
  otu_table(m, taxa_are_rows = taxa_are_rows)  # ASVs in columns, samples in rows
}

n_dupl <- function(x){
  x %>% 
    map(row.names) %>% 
    unlist %>% 
    duplicated() %>% 
    sum()
}

## Options ####
setwd("~/Desktop/GIT_REPOSITORIES/lter-sparc-above-below-synchrony/sandbox/geoff/")
# LOAD DATA ####

## Find ASV tables ####
asv_fps <- list.files("./output",full.names = TRUE, recursive = TRUE)
bact_asv_fps <- asv_fps[grepl(asv_fps,pattern="bacteria_seqtab")]
fung_asv_fps <- asv_fps[grepl(asv_fps,pattern="fungi_seqtab")]

## Read in tables ####
bact_asv_tables <- map(bact_asv_fps,readRDS)
fung_asv_tables <- map(fung_asv_fps,readRDS)

# Check for any duplicated row names
n_dupl(bact_asv_tables)
n_dupl(fung_asv_tables)

## Convert to otu_table objects (phyloseq)
bact_asvs <- map(bact_asv_tables,as_otu)
fung_asvs <- map(fung_asv_tables,as_otu)

# Clean up 
rm(list = c("bact_asv_tables","fung_asv_tables"))
gc()

# COMBINE TABLES ####

# combine bacteria
bact <- bact_asvs %>% reduce(merge_phyloseq)
# save to output
saveRDS(bact,"./output/bacterial_asv_table_combined.RDS")
# clean up
rm(list=c("bact","bact_asvs"))
gc()

# combine fungi
fung <- fung_asvs %>% reduce(merge_phyloseq)
# save to output
saveRDS(fung,"./output/fungal_asv_table_combined.RDS")
# clean up
rm("fung")
gc()
