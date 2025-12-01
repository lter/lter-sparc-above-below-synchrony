# SETUP ####
# WARNING: this build step needs huge RAM, don't rerun unless absolutely necessary

## Packages ####
library(dada2)
library(phyloseq)
library(tidyverse)
library(Biostrings)

## Functions ####
as_otu <- function(x, taxa_are_rows = FALSE) {
  m <- as.matrix(x, drop = FALSE)
  # storage.mode(m) <- "integer"
  # remove low abundance asvs
  good_asvs <- colSums(m) > 9
  m <- m[,good_asvs, drop = FALSE]
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

# LOAD DATA ####

## Get metadata ####
meta <- readRDS("./full_microbiome_dataset_info.RDS")

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
bact <- bact_asvs %>% purrr::reduce(merge_phyloseq)

# identify singleton taxa/samples and remove
bact_good_samples <- rowSums(bact) > 1
bact_good_taxa <- colSums(bact) > 1
bact <- bact[bact_good_samples,bact_good_taxa]

# indentify metadata rows for bacteria
bact_meta <- 
  meta %>% 
  dplyr::filter(rawDataFileName %in% row.names(bact)) %>% 
  dplyr::select(plotID,year,m_y,dnaSampleID,
                siteID,collectDate,forwardPrimer,
                reversePrimer,sequencerRunID,
                rawDataFileName)

# build metadata for physeq
bact_samdata <- sample_data(bact_meta)
sample_names(bact_samdata) <- bact_samdata$rawDataFileName

# build physeq object
bact_ps <- phyloseq(bact_samdata,bact)
bact_ps@sam_data$rawDataFileName <- NULL

# save to output
saveRDS(bact_ps,"./output/bacterial_physeq_no_taxonomy.RDS")
saveRDS(bact,"./output/bacterial_asv_table_combined.RDS")

# write ASVs to file for external taxonomic assignment
bact_reads <- colnames(otu_table(bact_ps))
names(bact_reads) <- bact_reads
writeXStringSet(DNAStringSet(bact_reads), "./output/bact_asv.fasta", width=200000)

# clean up
rm(list=c("bact","bact_asvs"))
gc()

# combine fungi
fung <- fung_asvs %>% purrr::reduce(merge_phyloseq)

# identify singleton taxa/samples and remove
fung_good_samples <- rowSums(fung) > 1
fung_good_taxa <- colSums(fung) > 1
fung <- fung[fung_good_samples,fung_good_taxa]

# indentify metadata rows for fungi
fung_meta <- 
  meta %>% 
  dplyr::filter(rawDataFileName %in% row.names(fung)) %>% 
  dplyr::select(plotID,year,m_y,dnaSampleID,
                siteID,collectDate,forwardPrimer,
                reversePrimer,sequencerRunID,
                rawDataFileName)

# build metadata for physeq
fung_samdata <- sample_data(fung_meta)
sample_names(fung_samdata) <- fung_samdata$rawDataFileName

# build physeq object
fung_ps <- phyloseq(fung_samdata,fung)
fung_ps@sam_data$rawDataFileName <- NULL

# save to output
saveRDS(fung_ps,"./output/fungal_physeq_no_taxonomy.RDS")
saveRDS(fung,"./output/fungal_asv_table_combined.RDS")

# write ASVs to file for external taxonomic assignment
fung_reads <- colnames(otu_table(fung_ps))
names(fung_reads) <- fung_reads
writeXStringSet(DNAStringSet(fung_reads), "./output/fung_asv.fasta", width=200000)


# clean up
rm("fung")
gc()
