# SETUP ####

## Packages ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")

## Options ####
set.seed(666)
nthreads <- parallel::detectCores()

seqrun <- "REPLACEME"


## Load metadata ####
paired <- readRDS("./paired_data.RDS")

# SUBSET TO SEQRUN ####
current_run <- 
  paired %>% 
  dplyr::filter(sequencerRunID == seqrun)

# SUBSET TO TAXONOMIC DOMAIN ####
current_bact <- 
  current_run %>% 
  dplyr::filter(amplicon == "Bacteria")

current_fung <- 
  current_run %>% 
  dplyr::filter(amplicon == "Fungi")


# BACTERIAL PIPELINE ####
if(nrow(current_bact) > 0){
  
  ## Get filepaths
  fwd_fps <- current_bact$clean_fwd_filepath
  rev_fps <- current_bact$clean_rev_filepath
  
  ## Learn errors
  bact_errF <- learnErrors(fwd_fps, 
                           multithread=nthreads, 
                           MAX_CONSIST = 10,
                           verbose = 1,
                           randomize = TRUE)
  bact_errR <- learnErrors(rev_fps, 
                           multithread=nthreads, 
                           MAX_CONSIST = 10,
                           verbose = 1,
                           randomize = TRUE)
  # Dereplicate
  bact_derepF <- derepFastq(fwd_fps, verbose=FALSE)
  bact_derepR <- derepFastq(rev_fps, verbose=FALSE)
  
  # Sample inference
  bact_dadaF <- dada(bact_derepF, err=bact_errF, multithread=TRUE, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
  bact_dadaR <- dada(bact_derepR, err=bact_errR, multithread=TRUE, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
  
  # Merge
  mergers <- mergePairs(bact_dadaF, fwd_fps, bact_dadaR, rev_fps, verbose=TRUE)
  
  # Make seq table
  bact_seqtab <- makeSequenceTable(mergers)
  bact_seqtab.nochim <- removeBimeraDenovo(bact_seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
  
  # save output
  saveRDS(bact_seqtab.nochim,file.path("output",seqrun,"Bacteria",paste0(seqrun,"_Bacteria_seqtab.RDS")))
  
}


# FUNGAL PIPELINE ####
if(nrow(current_fung) > 0){
  
  ## Get filepaths
  fwd_fps <- current_fung$clean_fwd_filepath
  rev_fps <- current_fung$clean_rev_filepath
  
  ## Learn errors
  fung_errF <- learnErrors(fwd_fps, 
                           multithread=nthreads, 
                           MAX_CONSIST = 10,
                           verbose = 1,
                           randomize = TRUE)
  fung_errR <- learnErrors(rev_fps, 
                           multithread=nthreads, 
                           MAX_CONSIST = 10,
                           verbose = 1,
                           randomize = TRUE)
  # Dereplicate
  fung_derepF <- derepFastq(fwd_fps, verbose=FALSE)
  fung_derepR <- derepFastq(rev_fps, verbose=FALSE)
  
  # Sample inference
  fung_dadaF <- dada(fung_derepF, err=fung_errF, multithread=TRUE, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
  fung_dadaR <- dada(fung_derepR, err=fung_errR, multithread=TRUE, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
  
  # Merge
  mergers <- mergePairs(fung_dadaF, fwd_fps, fung_dadaR, rev_fps, verbose=TRUE)
  
  # Make seq table
  fung_seqtab <- makeSequenceTable(mergers)
  fung_seqtab.nochim <- removeBimeraDenovo(fung_seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
  
  # save output
  outdir <- file.path("output",seqrun,"fungi")
  if(!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}
  saveRDS(fung_seqtab.nochim,file.path(outdir,paste0(seqrun,"_fungi_seqtab.RDS")))
  
}
