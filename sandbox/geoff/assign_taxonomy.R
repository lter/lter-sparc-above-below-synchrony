# SETUP ####

## packages ####
library(tidyverse)
library(dada2)
library(phyloseq)

fung_db <- "/sciclone/scr10/gzahn/DATABASES/Eukaryome/Eukaryome_General_SSU_v1.8_reformatted.fasta.gz"
bact_db <- "/sciclone/scr10/gzahn/DATABASES/Silva/silva_nr99_v138.2_toSpecies_trainset.fa.gz"
bact_spp_db <- "/sciclone/scr10/gzahn/DATABASES/Silva/silva_v138.2_assignSpecies.fa.gz"

## Load physeq objects and extract reads

bact <- readRDS("./output/bacterial_physeq_no_taxonomy.RDS") 
bact_reads <- bact %>% otu_table() %>% colnames()
fung <- readRDS("./output/fungal_physeq_no_taxonomy.RDS")
fung_reads <- fung %>% otu_table() %>% colnames()


# ASSIGN TAXONOMY ####

## bacteria ####
bact_tax <- 
  assignTaxonomy(bact_reads,
                 refFasta = bact_db,
                 multithread = TRUE,
                 tryRC = FALSE,
                 minBoot = 50)

# run addSpecies algorithm and add to previous tax table
bact_tax_spp <- addSpecies(taxtab = bact_tax,refFasta = bact_spp_db,tryRC = TRUE)
# save progress
saveRDS(bact_tax_spp,"./output/bact_tax_table.RDS")


## fungi ####
fung_tax <- 
  assignTaxonomy(fung_reads,
                 refFasta = euk_db,
                 multithread = TRUE,
                 tryRC = FALSE,
                 minBoot = 50)
# save progress
saveRDS(fung_tax,"./output/fung_tax_table.RDS")


# ADD TO PHYSEQ ####
bact_ps <- phyloseq(bact, tax_table(bact_tax))
fung_ps <- phyloseq(fung, tax_table(fung_tax))

# remove non-fungi/non-bacteria
bact_ps <- 
  bact_ps %>% 
  subset_taxa(Kingdom == "Bacteria" & Family != "Mitochondria" & Order != "Chloroplast")

fung_ps <- 
  fung_ps %>% 
  subset_taxa(Kingdom == "Fungi" & !is.na(Phylum))



# EXPORT ####
saveRDS(bact_ps, "./output/bacterial_physeq.RDS")
saveRDS(fung_ps, "./output/fungal_physeq.RDS")
