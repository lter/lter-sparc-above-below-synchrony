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

# SUBSET TO SEQRUN (and year) ####
current_run <- 
  paired %>% 
  dplyr::filter(sequencerRunID == seqrun & year == 2020)

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
  
  fwd_fps <- fwd_fps[file.exists(fwd_fps) & file.exists(rev_fps)]
  rev_fps <- rev_fps[file.exists(fwd_fps) & file.exists(rev_fps)]
  
  # Filter and trim
  fwd_filt <- fwd_fps %>% str_replace("/clean/","/filtered/")
  rev_filt <- rev_fps %>% str_replace("/clean/","/filtered/")
  for(i in unique(dirname(fwd_filt))){dir.create(i,recursive = TRUE,showWarnings = FALSE)}
  

  dada2::filterAndTrim(fwd = fwd_fps, filt = fwd_filt,
                       rev = rev_fps, filt.rev = rev_filt,
                       compress = TRUE)
  remaining_files <- file.exists(fwd_filt) & file.exists(rev_filt)
  
  fwd_filt <- fwd_filt[remaining_files]
  rev_filt <- rev_filt[remaining_files]
  
  ## Learn errors
  bact_errF <- learnErrors(fwd_filt, 
                           multithread=nthreads, 
                           MAX_CONSIST = 10,
                           verbose = 1,
                           randomize = TRUE)
  bact_errR <- learnErrors(rev_filt, 
                           multithread=nthreads, 
                           MAX_CONSIST = 10,
                           verbose = 1,
                           randomize = TRUE)
  # Dereplicate
  bact_derepF <- derepFastq(fwd_filt, verbose=FALSE)
  bact_derepR <- derepFastq(rev_filt, verbose=FALSE)
  
  # Sample inference
  bact_dadaF <- dada(bact_derepF, err=bact_errF, multithread=nthreads, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
  bact_dadaR <- dada(bact_derepR, err=bact_errR, multithread=nthreads, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
  
  # Merge
  mergers <- mergePairs(bact_dadaF, fwd_filt, bact_dadaR, rev_filt, verbose=TRUE)
  
  # Make seq table
  bact_seqtab <- makeSequenceTable(mergers)
  bact_seqtab.nochim <- removeBimeraDenovo(bact_seqtab, method="consensus", multithread=nthreads, verbose=FALSE)
  
  # save output
  outdir <- file.path("output",seqrun,"bacteria")
  if(!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}
  saveRDS(bact_seqtab.nochim,file.path(outdir,paste0(seqrun,'_2020',"_bacteria_seqtab.RDS")))
}


# FUNGAL PIPELINE ####
# only using forward reads for fungi
if(nrow(current_fung) > 0){
  
  ## Get filepaths
  fwd_fps <- current_fung$clean_fwd_filepath

  fwd_fps <- fwd_fps[file.exists(fwd_fps)]
  
  # Filter and trim
  fwd_filt <- fwd_fps %>% str_replace("/clean/","/filtered/")
  for(i in unique(dirname(fwd_filt))){dir.create(i,recursive = TRUE,showWarnings = FALSE)}
  

  dada2::filterAndTrim(fwd = fwd_fps, filt = fwd_filt,
                       compress = TRUE)
  fwd_filt <- fwd_filt[file.exists(fwd_filt)]
  

  ## Learn errors
  fung_errF <- learnErrors(fwd_filt, 
                           multithread=nthreads, 
                           MAX_CONSIST = 10,
                           verbose = 1,
                           randomize = TRUE)

  # Dereplicate
  fung_derepF <- derepFastq(fwd_filt, verbose=FALSE)
  
  # Sample inference
  fung_dadaF <- dada(fung_derepF, err=fung_errF, multithread=nthreads, selfConsist = FALSE, verbose=TRUE, pool = "pseudo")
  
  # Merge
  
  # Make seq table
  fung_seqtab <- makeSequenceTable(fung_dadaF)
  fung_seqtab.nochim <- removeBimeraDenovo(fung_seqtab, method="consensus", multithread=TRUE, verbose=FALSE)
  
  # save output
  outdir <- file.path("output",seqrun,"fungi")
  if(!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}
  saveRDS(fung_seqtab.nochim,file.path(outdir,paste0(seqrun,'_2020',"_fungi_seqtab.RDS")))
  
}
