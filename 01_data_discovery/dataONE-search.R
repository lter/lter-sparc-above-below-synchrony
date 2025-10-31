library(tidyverse)
library(xml2)
library(httr2)

# Start here for solr query info:
# https://solr.apache.org/guide/solr/latest/query-guide/common-query-parameters.html

#(above* OR plant* OR vegetation) AND (below* OR soil OR rhizosphere OR microb*)

# Create the query string
# "text"  is the "full text" search index
querystring = 'q=text:(above* OR plant* OR vegetation) AND text:(below* OR soil OR rhizosphere OR microb*) AND text:("paired" OR "coupled" OR "linked") AND (-obsoletedBy:* AND formatType:METADATA)'
# Parameters specific to DataONE
params <- 'fl=id,title,formatId&rows=100000&wt=xml'
# Colons, commas, and spaces need to be encoded 
querystring <- gsub(':', '%3A', gsub(',', '%2C', gsub(' ', '%20', paste(querystring, params, sep='&'))))
querystring

# This gives xml back, but only 10000 rows. It counts full result tho
result_xml  <- httr2::request(paste0('https://cn.dataone.org/cn/v2/query/solr/?', querystring)) |>
  httr2::req_perform() |>
  httr2::resp_body_string() |>
  read_xml()

# Check results
numfound <- xml_attr(xml_find_all(result_xml, ".//result"), 'numFound')
numreturned <- length(xml_find_all(result_xml, ".//doc"))

# This works if wt=csv, but there are no header results
# and the max return is 1000 rows
result_csv <- httr2::request(paste0('https://cn.dataone.org/cn/v2/query/solr/?', querystring)) |>
  httr2::req_perform() |>
  httr2::resp_body_string() |>
  readr::read_csv()

# This gives json back, but only 10000 rows.
# result_json <- httr2::request(paste0('https://cn.dataone.org/cn/v2/query/solr/?', querystring)) |>
#   httr2::req_perform() |>
#   httr2::resp_body_json()# |>
  #readr::read_csv()


# If using an xml result, parse to CSV
# First get a list of all 'str' notes within 'doc' nodes,
# then filter by attribute name and extract content
title <- xml_find_all(result_xml, xpath = "//doc//str[@name='title']") |> xml_text()
id <- xml_find_all(result_xml, xpath = "//doc//str[@name='id']") |> xml_text()
formatId <- xml_find_all(result_xml, xpath = "//doc//str[@name='formatId']") |> xml_text()
# Format into tibble
df <- tibble(id = id, title = title, formatId = formatId)
df

