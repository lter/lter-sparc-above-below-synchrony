# clean up raw neon reads

# SETUP ####

## Packages ####
library(tidyverse) 
library(dada2) 
library(ShortRead) 
library(Biostrings) 
library(glue)

## Environment ####
nthreads <- max(1, parallel::detectCores())

## Functions ####

# build all orientations of dna seqs
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}


## Load metadata ####
full <- readRDS("./full_microbiome_dataset_info.RDS")

# get directory structure
dirs <- list.dirs("./data",recursive = TRUE)
bact_dirs <- dirs[grepl("Bacteria",dirs)]
fung_dirs <- dirs[grepl("Fungi",dirs)]


# each individual paired sample needs to have its own cutadapt command
# they could have different primers, for example

# remove any that have missing samples
duplicated(full$dnaSampleID)

full %>% 
  filter(dnaSampleID == "HARV_033-O-33.5-20.5-20170424-GEN-DNA1")

full$dnaSampleID %>% table %>% as.data.frame %>% 
  arrange(Freq) %>% 
  filter(Freq %in% c(1,3,5,7,9))


# find bacterial read pair filepaths
# x <- full[which(full$dnaSampleID == "HARV_033-O-33.5-20.5-20170424-GEN-DNA1"),] 
x <- full
# change to filepath to match actual relative location paths
x$rawDataFilePath <- file.path("./data/raw",x$sequencerRunID,x$amplicon,x$rawDataFileName)

# Normalize read direction using both the 'read' column and filename hints
df_clean <- x %>%
  mutate(
    read_dir = case_when(
      str_to_lower(read) %in% c("forward","r1") ~ "Forward",
      str_to_lower(read) %in% c("reverse","r2") ~ "Reverse",
      str_detect(rawDataFileName, "(^|_)R1(\\.|_|$)") ~ "Forward",
      str_detect(rawDataFileName, "(^|_)R2(\\.|_|$)") ~ "Reverse",
      TRUE ~ NA_character_
    ),
    # Some runs happen multiple times for the same dnaSampleID; keep sequencerRunID to separate them
    run_key = paste(dnaSampleID, amplicon, sequencerRunID, sep = "||")
  )

# If there are accidental duplicates per read within the same run, keep a single "best" row
#    (Here we keep the most recent-looking path)
df_dedup <- df_clean %>%
  arrange(desc(rawDataFilePath)) %>%  # crude recency proxy if years embedded in URL
  group_by(run_key, read_dir) %>%
  slice_head(n = 1) %>%
  ungroup()

# Choose which metadata to retain
meta_cols <- c(
  "plotID","year","m_y","dnaSampleID","siteID","collectDate",
  "forwardPrimer","reversePrimer","amplicon","sequencerRunID"
)

# 4) Pivot to forward/reverse paths
paired <- df_dedup %>%
  select(all_of(meta_cols), rawDataFilePath, read_dir, run_key) %>%
  pivot_wider(
    id_cols  = all_of(meta_cols),
    names_from  = read_dir,
    values_from = rawDataFilePath,
    values_fn = list  # in case there are still dupes; keeps first by default
  ) %>%
  # If lists were created, unlist one value
  mutate(
    forward_path = if (is.list(Forward)) sapply(Forward, `[`, 1) else Forward,
    reverse_path = if (is.list(Reverse)) sapply(Reverse, `[`, 1) else Reverse
  ) %>%
  select(all_of(meta_cols), forward_path, reverse_path)
paired

# clean up samples with a missing data (either fwd or rev)
missing_files <- which(paired$forward_path %>% map_lgl(is.null) | paired$reverse_path %>% map_lgl(is.null))
paired <- paired[-missing_files,]
missing_primers <- which(is.na(paired$forwardPrimer) | is.na(paired$reversePrimer))
paired <- paired[-missing_primers,]


# clean up list cols
paired$forward_path <- unlist(paired$forward_path)
paired$reverse_path <- unlist(paired$reverse_path)

# now we have paired files for cutadapt

# for each ROW (DNA reads sample), build the cutadapt call based on
# forwardPrimer, reversePrimer, sequencerRunID, forward_path, reverse_path


# for each bacterial run (bact_dirs), run cutadapt
# put output clean files into "./data/clean/RUNID/Bacteria"

for(samp in seq_along(paired$forward_path)){
  
  # get primer info
  fwd_primer <- paired$forwardPrimer[samp]
  rev_primer <- paired$reversePrimer[samp]
  fwd_primer_rc <- rc(fwd_primer)
  rev_primer_rc <- rc(rev_primer)
  
  # get raw fp info
  fwd_path <- paired$forward_path[samp]
  rev_path <- paired$reverse_path[samp]
  
  # build clean fp info
  fwd_path_out <- str_replace(fwd_path,"data/raw","data/clean")
  rev_path_out <- str_replace(rev_path,"data/raw","data/clean")
  
  # create output directory if doesn't exist
  if(!dir.exists(dirname(fwd_path_out))){
    dir.create(dirname(fwd_path_out),recursive = TRUE)
  }
  
  # add new paths of cleaned files to paired dataframe
  paired$clean_fwd_filepath <- NA
  paired$clean_rev_filepath <- NA
  paired$clean_fwd_filepath[samp] <- fwd_path_out
  paired$clean_rev_filepath[samp] <- rev_path_out
  
  # cutadapt flags
  R1.flags <- paste("-g", fwd_primer, "-a", rev_primer_rc)
  R2.flags <- paste("-G", rev_primer, "-A", fwd_primer_rc) 
  
  # if the clean file exists already, skip
  if(file.exists(fwd_path_out) & file.exists(rev_path_out)){
    
    system2("cutadapt", args = c(R1.flags, R2.flags, "-n", 2, "--minimum-length 100", # -n 2 required to remove FWD and REV from reads
                                 "-o", fwd_path_out, "-p", rev_path_out, # output files
                                 fwd_path, rev_path)) # input files
  }
}

# export new seq metadata file
saveRDS(paired,"./paired_data.RDS")


