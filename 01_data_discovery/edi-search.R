library(tidyverse)
library(EDIutils)
#library(httr2)
readRenviron(".Renviron")

# Start here for solr query info:
# https://solr.apache.org/guide/solr/latest/query-guide/common-query-parameters.html
# PASTA solr query info: https://pastaplus-core.readthedocs.io/en/latest/doc_tree/pasta_api/data_package_manager_api.html#browse-and-discovery

datadir <- file.path(Sys.getenv("SABRE_ROOT"), "01_DATA")

# Some possible queries
#(above* OR plant* OR vegetation) AND (below* OR soil OR rhizosphere OR microb*)

# Create the queries in a named list
# "subject" is the "full text" search index for EDI, though only inlcudes title, keywords, and abstract
# Here is a broad search
queries <- list( 
  # first a broad search
  #"broad" = "q=subject:(above* OR plant* OR vegetation) AND subject:(below* OR soil OR rhizosphere OR microb*)  AND subject:(paired OR coupled OR linked)",
  # microbial biomass
  "microb*&biomass" = "q=subject:(microb* AND biomass)", 
  # microbial community
  "microb*&comm" = "q=subject:(microb* AND subject:(diversity OR community OR abundance OR composition OR occurrence OR 16S OR ITS OR bacter* OR fung*))",
  # aboveground biomass, npp, etc
  "above*&biomass" = "q=subject:(subject:(plant* OR vegetation OR aboveground OR above-ground OR 'above ground') AND subject:(biomass OR NPP OR producti* OR growth OR exchange))",
  # aboveground cover, community composition, diversity, etc
  "above*&comm" = "q=subject:(subject:(plant* OR vegetation OR aboveground OR above-ground OR 'above ground') AND subject:(community OR cover OR composition OR *diversity))"
)

# Search parameters specific to EDI's solr instance - apply to all searches
params = 'fl=packageid,title,doi&fq=-scope:(ecotrends+lter-landsat*+knb-lter-mcm)'

# Search with EDIutils function
# Note that spaces must be encoded with "%20"
runquery <- function(queryname, querystring){
  result <- EDIutils::search_data_packages(query = paste(gsub(' ', '%20', querystring), params, sep='&')) %>%
  # Format the resulting table
  mutate(landing_url = paste0("https://portal.edirepository.org/nis/mapbrowse?packageid=", packageid),
         metadata_std = "https://eml.ecoinformatics.org/eml-2.2.0",
         # Reformat doi to a URL
         doi = gsub('doi:', 'https://doi.org/', doi),
         repository = 'EDI', query_name=queryname) %>%
  select(repository, repository_id=packageid, landing_url, title, metadata_std, doi, query_name)
  return(result)
}
# Run the queries and get a list of dataframe results
results <- list()
for(q in names(queries)){
  print(q)
  results[[q]] <- runquery(q, queries[q])
  readr::write_csv(results[[q]], file.path(datadir, "01_Raw_data",
      "raw_search_results", paste0("edi-results-", q ,format(Sys.Date(),
                               "_%Y%m%d"), ".csv")))
}

# Create a dataframe with deduplicated search results and columns for 
# an AI agent to fill in with metadata extracted from EML files (see
# AGENTS.md and PLAN.md files)
results_dedup <- bind_rows(results) |> filter(!duplicated(doi)) |>
  mutate(
    metadata_url = paste0("https://portal.edirepository.org/nis/metadataviewer?packageid=",
      repository_id, "&contentType=application/xml"),
    study_type = NA,
    obs_or_exp = NA,
    site_type = NA,
    site_name = NA,
    ecosystem_type = NA,
    aquatic_bool = NA,
    belowground_vars = NA,
    aboveground_vars = NA,
    soil_depth = NA,
    bb_west = NA,
    bb_east = NA,
    bb_south = NA,
    bb_north = NA,
    year_start = NA,
    year_end = NA,
    sampling_freq = NA,
    sub_1_year_bool = NA,
    license = NA,
    claude_eval_date = NA,
    claude_uncertainty = NA
)
# Write deduplicated file with empty columns
readr::write_csv(results_dedup, file.path(datadir, "01_Raw_data",
    paste0('combined_edi_results_', format(Sys.Date(), "%Y%m%d"), "_aiready.csv")))

# Try with httr2 - for the record, this gives the same result as EDIutils
# rq <- httr2::request(paste0("https://pasta.lternet.edu/package/search/eml?",
#                             paste(gsub(' ', '%20', querystring), params, sep='&')))
# xml <- httr2::resp_body_xml(httr2::req_perform(rq))
# xml
