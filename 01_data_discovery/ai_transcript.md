**Prompt:** I have now written a work plan in the `plan.md` document. I would like you  
  to complete Phase 1 of the plan. You'll be reading EML metadata in the      
  `eml/` subdirectory and using extracted information to fill in the table in 
   `combined_edi_bg_results_20260422.csv`. Before you begin read the plan and 
   ask me any clarifying questions you have.

**Claude's clarifying questions, and my answers**

* I have a good picture of the plan and the EML structure. Before starting, I
  have a few clarifying questions:

  1. claude_eval_date and claude_uncertainty columns: The CSV has two columns
  not in the plan. Should I populate these? For claude_eval_date I'd use
  today's date (2026-04-22). For claude_uncertainty, what scale or values
  would you like — e.g., low/medium/high to flag how confident I am on the
  inferred fields (study_type, obs_or_exp, ecosystem_type, etc.)?
  2. site_type when not a known network: For datasets that aren't LTER, NEON,
  or another recognizable network (e.g., standalone edi.* packages), should I
  put NA or something else?
  3. site_name for non-LTER/NEON datasets: The plan specifies 3-letter LTER
  acronyms and 4-letter NEON acronyms. For datasets outside those networks,
  should I use a short descriptive name derived from the metadata, or NA?
  4. belowground_vars / aboveground_vars naming: Should I use the actual
  attribute/column names from the data tables (e.g., BMC, BMN), or
  interpretive descriptive names (e.g., microbial biomass carbon)?
  5. Batch review first?: There are 450 datasets. Would you like me to process
   a small sample (e.g., 10–20) for your review before running the full set?

❯ 1. I just saved the plan document again. It has instructions on filling in  
  those two columns. 2. Put an NA if they are not LTER or NEON. 3. Use a      
  short descriptive name from the metadata. 4. Use interpretive descriptive   
  names. We can look at the metadata later to find the actual column name. 5. 
   Yes! Process the first 20 and then let me review.  

**Prompt after Claude did the first 20**: This looks very good so far. Continue to make notes when important things are missing, like geographic coverages. Proceed with the next 430 rows.