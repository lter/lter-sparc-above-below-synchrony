library(tidyverse)
library(xml2)
library(httr2)
      
# Start here for solr query info:
# https://solr.apache.org/guide/solr/latest/query-guide/common-query-parameters.html

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#(above* OR plant* OR vegetation) AND (below* OR soil OR rhizosphere OR microb*)

# Create the query string
# "text"  is the "full text" search index
querystring = 'q=text:(above* OR plant* OR vegetation) AND text:(below* OR soil OR rhizosphere OR microb*) AND text:("paired" OR "coupled" OR "linked") AND (-obsoletedBy:* AND formatType:METADATA)'

# Parameters specific to DataONE
params <- 'fl=id,title,formatId,seriesId&rows=100000&wt=xml' #&start=9999'
# Colons, commas, and spaces need to be encoded 
querystring <- gsub(':', '%3A', gsub(',', '%2C', gsub(' ', '%20', paste(querystring, params, sep='&'))))
querystring

# This gives xml back, but only 10000 rows. It counts full result tho
result_xml  <- httr2::request(paste0('https://cn.dataone.org/cn/v2/query/solr/?', querystring)) |>
  httr2::req_perform() |>
  httr2::resp_body_string() |>
  read_xml()

# Check results. If there are more than 10000 found, collect additional 
# results using the "start" query parameter.
numfound <- xml_attr(xml_find_all(result_xml, ".//result"), 'numFound')
numreturned <- length(xml_find_all(result_xml, ".//doc"))

# This works if wt=csv, but there are no header results
# and the max return is 1000 rows
# result_csv <- httr2::request(paste0('https://cn.dataone.org/cn/v2/query/solr/?', querystring)) |>
#   httr2::req_perform() |>
#   httr2::resp_body_string() |>
#   readr::read_csv()

# This gives json back, but only 10000 rows.
# result_json <- httr2::request(paste0('https://cn.dataone.org/cn/v2/query/solr/?', querystring)) |>
#   httr2::req_perform() |>
#   httr2::resp_body_json()# |>
  #readr::read_csv()


# If using an xml result, parse to CSV
# First get a list of all 'doc' nodes, then find the content
# by attribute name. XML_find_first eturns NA when needed.
docs <- xml_find_all(result_xml, xpath = "//doc")
id <- xml_find_first(docs, xpath = ".//str[@name='id']") |> xml_text()
landingURL <- paste0('https://search.dataone.org/view/', id)
title <- xml_find_first(docs, xpath = ".//str[@name='title']") |> xml_text()
metadataFormat <- xml_find_first(docs, xpath = ".//str[@name='formatId']") |> xml_text()
# SeriesId doesn't always contain a DOI... there will be empties
doi <- xml_find_first(docs, xpath = ".//str[@name='seriesId']") |> xml_text()
# Format the result into a tibble
result <- tibble(repoName = "DataONE", id = id, landingURL = landingURL, title = title,
             metadataFormat = metadataFormat, doi = doi) |>
  # Reformat doi to a URL
  mutate(doi = gsub('doi:', 'https://doi.org/', doi))

# Write result
readr::write_csv(result, 'dataONE-result.csv')
