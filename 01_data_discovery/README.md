# SABRE Data Discovery

These are data discovery materials for an ecological synthesis research project called SABRE - Synchrony of Aboveground and Belowground Responses across Ecosystems. The group has searched the EDI repository (https://edirepository.org) for datasets relevant to observing synchrony between above and belowground communities of organisms.

**Materials here:**

* `edi-search.R` is a query script written in R. It uses the EDI repository's [PASTA REST API](https://pastaplus-core.readthedocs.io/en/latest/doc_tree/pasta_api/index.html) to query an Apache Solr index for datasets in the repository matching some search terms relevant to the research goals of the group. There are several queries in it.
* `dataONE-search.R` is another R query script similar to the above but for the DataONE API.
* `explore_metadata.R` has preliminary analyses of the data discovery results
* `wget_eml.py` python script used to download EML files from EDI.
* AI agent prompting - Some work here was done by using AI.
    - The `PLAN.md` document is the agent work plan referred to when prompting Claude.
    - See `ai_transcript.md` for prompting details and chat transcripts.