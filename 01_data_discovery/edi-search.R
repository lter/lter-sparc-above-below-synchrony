library(tidyverse)
library(EDIutils)
library(httr2)

# Start here for solr query info:
# https://solr.apache.org/guide/solr/latest/query-guide/common-query-parameters.html

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Some possible queries
#(above* OR plant* OR vegetation) AND (below* OR soil OR rhizosphere OR microb*)

# Create the query string
# "subject" is the "full text" search index for EDI, though only inlcudes title, keywords, and abstract
querystring = 'q=subject:(above* OR plant* OR vegetation) AND subject:(below* OR soil OR rhizosphere OR microb*)  AND subject:("paired" OR "coupled" OR "linked")'
# Search parameters specific to EDI's solr instance
params = 'fl=packageid,title,doi&fq=-scope:(ecotrends+lter-landsat*)'

# Search with EDIutils, note that spaces must be encoded with "%20"
result <- EDIutils::search_data_packages(query = paste(gsub(' ', '%20', querystring), params, sep='&')) %>%
  # Format the resulting table
  mutate(landingURL = paste0("https://portal.edirepository.org/nis/mapbrowse?packageid=", packageid),
         metadataFormat = "https://eml.ecoinformatics.org/eml-2.2.0",
         # Reformat doi to a URL
         doi = gsub('doi:', 'https://doi.org/', doi),
         repoName = 'EDI') %>%
  select(repoName, id=packageid, landingURL, title, metadataFormat, doi )

# Write result
readr::write_csv(result, 'edi-result.csv')

# Try with httr2
# rq <- httr2::request(paste0("https://pasta.lternet.edu/package/search/eml?",
#                             paste(gsub(' ', '%20', querystring), params, sep='&')))
# xml <- httr2::resp_body_xml(httr2::req_perform(rq))
# xml
# Same result
