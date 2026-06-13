# Work plan for SABRE data discovery analyses

Currently we've done four searches in the EDI repository (https://portal.edirepository.org) for data focusing on belowground microbial commmunity composition and belowground ecosystem processes (biomass, fluxes, etc.) and aboveground plant biomass, growth, or community composition. There are 3238 potential datasets to explore from these search results. We have a few more phases of work, including data and metadata aggregation and exploration, additional searches, pairing of appropriate datasets for synchrony analysis, and then harmonization of the paired datasets.

## Phase 1

This phase focuses on exploring and collating the metadata for the current belowground EDI search results. Those results are in the `combined_edi_bg_results_YYYYMMDD.csv` (most recent version), which is a listing of datasets by repository, identifier, title, DOI, and more. We have added a number of empty columns to this table to fill with metadata that comes from the EML metadata that accompanies each dataset at the EDI repository. Those EML files have been downloaded to the `eml/` subdirectory and are named according to the identifier of the dataset (`id` column in the table.) We need to populate the empty columns of the table with metadata taken from the EML for each dataset according to the specifications for each column shown below. For reference, the EML schema is [here](https://eml.ecoinformatics.org/schema/).

1. `repository`: Do not edit. The name of the repository the dataset comes from. (all EDI at this point)
2. `repository_id`: Do not edit. Identifier of the dataset at the repository. Metadata files for each dataset will have filenames matching this value.
3. `landing_url`: Do not edit. URL for the dataset landing page
4. `title`: Do not edit. Title of the dataset. Can be found in the dataset title element of the metadata.
5. `metadata_std`: Do not edit. Name or URL of the metadata standard for the dataset.
6. `doi`: Do not edit. DOI link for the dataset.
7. `query`: Do not edit. Name of the query used to find the data in the repository.
8. `metadata_url`: Do not edit. URL for the metadata file accompanying the dataset
9. `study_type` : Infer a broad study type based methods metadata and fill in using: field, greenhouse, transplant, mesocosm, laboratory, other.
10. `obs_or_exp` : Infer from methods metadata and fill in using: observational (monitoring or measurements only, no experimental treatments applied), experiment (manipulative or experimental treatments were applied to some or all observational units)
11. `site_type` : Infer from title, abstract, methods, or geographic coverage metadata and fill in using: LTER (for any LTER site), NEON (for any NEON site), or any other well-known research network acronym
12. `site_name` : Fill in with a short name or acronym for the site. For LTER sites use the 29, three-letter acronyms, for NEON sites use the 4-letter acronyms. If there are more than one site names separate them with semi-colons
13. `ecosystem_type` : Fill in with an ecosystem type using the closest match of: forest, desert, grassland, ocean, freshwater, tundra, wetlands, savannah, mountain
14. `aquatic_bool` : Fill in with true or false based on whether this is an aquatic ecosystem
15. `belowground_vars` : Fill in name of any biotic variables measured belowground (or originating belowground) that are of interest for synchrony measurements. These may include microbial biomass, abundance, community composition, soil respiration, and more. Separate multiples with semicolons. These come from the data table attributes or the methods metadata.
16. `aboveground_vars` : Fill in name of any biotic variables measured aboveground that are of interest for synchrony measurements. These may include vegetation cover, plant community composition, plant biomass, plant production, and more. Separate multiples with semicolons. These come from the data table attributes or the methods metadata.
17. `soil_depth` : For any belowground variables, fill in the soil depth measured. Leave as NA if nothing is applicable.
18. `bb_north` : Fill in the most northern bounding coordinate from the geographic coverage metadata.
19. `bb_south` : Fill in the most southern bounding coordinate from the geographic coverage metadata.
20. `bb_east` : Fill in the most eastern bounding coordinate from the geographic coverage metadata.
21. `bb_west` : Fill in the most western bounding coordinate from the geographic coverage metadata.
22. `year_start` : Fill in the first year measurements began using the temporal coverage metadata.
23. `year_end` : Fill in the last year measurements were taken using the temporal coverage metadata.
24. `sampling_freq`: Fill in the sampling frequency using: annual (for yearly sampling), sub-annual (for more frequent than annual), other (for anything longer than annual)
25. `sub_1_year_bool` : Fill in with true or false if measurements lasted less than one year.
26. `license` : Use the intellectual rights and license metadata to fill in a license attached to the dataset using: CC-BY, CC0, unlicensed, or other.
27. `claude_eval`: Fill this this with today's date once metadata for the dataset has been evaluated by Claude and values in the row have been filled in.
28. `claude_uncertainty`: Enter a note about any uncertainties in the analysis of metadata that Claude found.


## Phase 2

In this phase we will collate some other search results.

## Phase 3
