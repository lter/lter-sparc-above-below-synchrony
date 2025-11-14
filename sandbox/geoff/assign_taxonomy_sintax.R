# SETUP ####

## Packages ####
library(tidyverse)
library(phyloseq)

## Functions ####
parse_sintax_vsearch <- function(path, asv_order,
                                 ranks = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
                                 cutoffs = c(Kingdom=.5, Phylum=.6, Class=.7, Order=.7, Family=.8, Genus=.85, Species=.9)) {
  # Map any rank code (lower/upper) to our rank names
  rank_map <- c(
    d = "Kingdom", k = "Kingdom",
    p = "Phylum",
    c = "Class",
    o = "Order",
    f = "Family",
    g = "Genus",
    s = "Species",
    D = "Kingdom", K = "Kingdom",
    P = "Phylum",
    C = "Class",
    O = "Order",
    F = "Family",
    G = "Genus",
    S = "Species"
  )
  
  # Read TSV with unknown number of columns, no names
  raw <- read_tsv(path, col_names = FALSE, show_col_types = FALSE, comment = "#")
  stopifnot(ncol(raw) >= 2)
  
  # Heuristically pick the column that contains rank tokens (k:/d:/p:/...)
  has_rank_tokens <- function(x) any(str_detect(x, "(^|[,\\s])(tax=)?[dkpcfgsoDKPCFGSO]:"))
  cand_cols <- which(vapply(raw, has_rank_tokens, logical(1)))
  stopifnot(length(cand_cols) >= 1)
  sintax_col <- cand_cols[1]             # first column with tokens
  qseqid_col <- 1                        # first column is your header (the ASV sequence)
  
  df <- tibble(
    qseqid = raw[[qseqid_col]],
    sintax = raw[[sintax_col]]
  )
  
  # Normalize to a common token list:
  # - remove optional 'tax=' prefix
  # - split on commas or whitespace between tokens
  # - tokens look like: 'd:Name(0.97)' OR 'k:Name,0.97'
  tokens <- df %>%
    mutate(s = str_remove(sintax, "^tax=")) %>%
    separate_rows(s, sep = "[,\\s]+") %>%              # split on comma or whitespace
    filter(str_detect(s, "^[dkpcfgsoDKPCFGSO]:")) %>% # keep rank tokens only
    mutate(code = str_sub(s, 1, 1),
           # taxon: after "code:" up to first '(' or ',' (whichever comes first)
           taxon = str_remove(str_extract(s, ":[^\\(,]+"), "^:"),
           # confidence: either inside (...) or after a trailing comma
           conf  = coalesce(
             suppressWarnings(as.numeric(str_match(s, "\\(([0-9.]+)\\)")[,2])),
             suppressWarnings(as.numeric(str_match(s, ",([0-9.]+)$")[,2]))
           ),
           rank = unname(rank_map[code])) %>%
    filter(!is.na(rank)) %>%
    transmute(qseqid = df$qseqid[match(qseqid, df$qseqid)], rank, taxon, conf)
  
  # Choose top per-rank (usually one) and apply per-rank cutoffs
  tax_wide <- tokens %>%
    group_by(qseqid, rank) %>%
    slice_max(order_by = conf, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(taxon = ifelse(!is.na(conf) & conf >= cutoffs[rank], taxon, NA_character_)) %>%
    select(qseqid, rank, taxon) %>%
    distinct(qseqid, rank, .keep_all = TRUE) %>%
    pivot_wider(names_from = rank, values_from = taxon)
  
  # Ensure all expected rank columns exist
  for (r in ranks) if (!r %in% names(tax_wide)) tax_wide[[r]] <- NA_character_
  
  # Left-join to full ASV list to preserve order and include missing classifications
  full <- tibble(qseqid = asv_order) %>% left_join(tax_wide, by = "qseqid")
  
  M <- as.matrix(full[, ranks])
  rownames(M) <- full$qseqid
  M
}

# READ TAXONOMY OUTPUT FROM VSEARCH ####
bact <- readRDS("./output/bacterial_physeq_no_taxonomy.RDS")
fung  <- readRDS("./output/fungal_physeq_no_taxonomy.RDS")

bact_reads <- colnames(otu_table(bact))
fung_reads <- colnames(otu_table(fung))

bact_tax <- parse_sintax_vsearch("./output/bact.sintax.tsv", asv_order = bact_reads)
fung_tax <- parse_sintax_vsearch("./output/fung.sintax.tsv", asv_order = fung_reads)

# lock in that the rows exactly match your ASV vectors
stopifnot(identical(rownames(bact_tax), bact_reads))
stopifnot(identical(rownames(fung_tax), fung_reads))

# add to physeq objects
bact_ps <- phyloseq(bact, tax_table(bact_tax))
fung_ps <- phyloseq(fung, tax_table(fung_tax))

# export
saveRDS(bact_ps,"./output/bacterial_physeq_with_taxonomy.RDS")
saveRDS(fung_ps,"./output/fungal_physeq_with_taxonomy.RDS")
