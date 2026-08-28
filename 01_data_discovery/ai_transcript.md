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

## 2026-08-21

I did another round of metadata annotation with Claude Sonnet 5. Note that the search terms have been updated to include aboveground searches in EDI, meaning the dataset to annotate is much bigger now (~3470 datasets).

**Prompt:** You are my detail-oriented metadata review and entry assistant. Please read PLAN.md and populate the metadata for the empty columns in combined_edi_results_20260820_aiready.csv with metadata extracted from the EML files in result_eml/. Follow instructions in PLAN.md carefully. Process the first twenty records first, let me review, then proceed with the rest.

After the initial 20 were reviewed, Claude asked what kind of workflow to use and a chose a multi-agent pipeline based on the current PLAN.md. Claude went off the rails for this and consumed all its tokens.

## 2026-08-26

I tried again with Deepseek v4 Pro in Opencode harness. The prompt was the same as above. The initial 20 went quite quickly as did annotating the remaining records. On review I spot-checked a few records and noticed that the `aboveground_vars` and `belowground_vars` were quite innacurate. I re-prompted with these spot checks and coached the model to refine the pattern matching. It took about three rounds of review and re-prompting to create the python script that worked well-enough.

FULL TRANSCRIPT:

**Date:** 2026-08-26

---

## Task

Populate empty metadata columns in `combined_edi_results_20260820_aiready.csv` from EML files in `result_eml/`, following `PLAN.md`. Process first 20 records for review, then all records.

---

## Phase 1: Understanding the task

- Read `PLAN.md`, which specifies 28 columns, of which columns 9–28 are to be filled:
  - `study_type`, `obs_or_exp`, `site_type`, `site_name`, `ecosystem_type`, `aquatic_bool`, `belowground_vars`, `aboveground_vars`, `soil_depth`, `bb_north/south/east/west`, `year_start/end`, `sampling_freq`, `sub_1_year_bool`, `license`, `claude_eval_date`, `claude_uncertainty`
- Read the CSV (3,747 data rows) and a sample EML (`knb-lter-knz.53.10.xml`) to learn XML structure.
- Discovered child elements in these EML files do NOT carry the `eml:` namespace prefix (only the root `eml:eml` does).

## Phase 2: Writing the extraction script

Wrote `extract_metadata.py` using `xml.etree.ElementTree` + regex with keyword heuristics for:
- study_type (field/greenhouse/transplant/mesocosm/laboratory/other)
- obs_or_exp (experiment/observational)
- site_type (LTER/NEON/NA)
- site_name (LTER 3-letter codes, NEON 4-letter codes, geographic description)
- ecosystem_type (forest/desert/grassland/ocean/freshwater/tundra/wetlands/savannah/mountain)
- aquatic_bool
- belowground/aboveground variables (from pattern matching)
- soil_depth, bounding coordinates, temporal coverage, sampling_freq, license

First run failed (namespace issue), fixed by stripping `eml:` prefix matching.

## Phase 3: First 20 records review

Generated `first20_preview.csv`. Identified several inference issues that were fixed:
- **Word boundary matching** (e.g., `"slaked"` → false "lake" match, `"Great Lakes"` → false "lake")
- **KBS LTER** forced to grassland/non-aquatic (agricultural LTER site)
- **NWT** → mountain, **SEV** → desert, **FCE/GCE** → wetlands, coastal LTERs → ocean/marine
- Institution names filtered from ecosystem matching (e.g., "Center for Ocean-Atmospheric Prediction")

## Phase 4: Full run (v1)

Processed all 3,747 records. Results:
- 3,745 success, 0 missing, 2 parse errors (empty EMLs: `edi.2304.1`, `edi.2369.1`)
- Output: `combined_edi_results_20260820_metadatacomplete.csv`

## Phase 5: Spot check #1 — jrn false positives

User flagged row 597/598 (`knb-lter-jrn.210527001.1`) had spurious `Microbial biomass C:N ratio` and `Microbial community composition`.

Root causes:
- `"biomass ratio"` matched `"leaf to stem biomass ratio"` (attribute `L:S`)
- `"its"` matched the English word `"its"` in "...by its dry mass..."

## Phase 6: Full pattern audit

Audited all bg/ag rules across the dataset. Identified critical overly-greedy patterns:

| Pattern | Rule | Problem |
|---------|------|---------|
| `"biomass"` (standalone) | ag: Plant biomass | 1,204 records; matches any biomass |
| `"community composition"` (standalone) | bg AND ag | 419 records got BOTH labels |
| `"its"` (standalone) | bg | common English word |
| `"biomass ratio"` (standalone) | bg | matches plant tissue ratios |
| `"enzyme"` / `"leachate"` / `"invertebrate"` / `"reproductive"` / `"survival"` / `"litter"` / `"water content"` standalone | various | literature citations / ambiguous context |

## Phase 7: Corrected rules (v2)

Applied corrected bg_rules and ag_rules. Re-ran.

**Result:**
- Microbial community composition: 1,649 → 117 (-93%)
- Plant biomass: 1,204 → 550 (-54%)
- Microbial biomass C:N ratio: 80 → 1 (-99%)
- Double-hit (both CC tags): 419 → 6 (-99%)
- Output: `combined_edi_results_20260820_metadatacomplete_v2.csv`

## Phase 8: Spot check #2 — cdr.405.8

User flagged row 61/62 (`knb-lter-cdr.405.8`) — a root-biomass-only dataset had several aboveground vars and no belowground.

Root causes:
1. **Title/entityName excluded** from variable search corpus — "root biomass" only appeared in `<entityName>e249_Root biomass data</entityName>` and the title, neither of which were searched.
2. **Methods text false matches**:
   - `"plant community composition"` from experimental design description
   - `"productivity"` from a literature citation ("Tilman et al. Diversity and productivity...")
   - `"lai"` from the word "laid"
   - `"plant c"` from substring of "plant community"

User advice: **de-emphasize methods element** — it contains context/citations, not measured variables. Agreed on a **primary-first with methods fallback** approach.

## Phase 9: Primary-corpus approach (v3)

Implemented:
1. **Primary corpus** = title + entityName + entityDescription + attributeName + attributeDefinition
2. **Phase 1**: match all patterns against primary corpus only
3. **Phase 2**: only if Phase 1 yields zero variables, fall back to adding methods text
4. Fixed patterns: removed `"lai"` (→ "leaf area index", "specific leaf area", "sla"), removed `"plant c"` / `"plant n"` (→ "plant carbon", "foliar carbon", etc.)

**Result:**
- cdr.405.8: `bg: Root biomass`, `ag: NA` ✅ (was `bg: NA`, `ag: Leaf area; Plant C; Plant CC; Plant production`)
- Plant biomass: 550 → 388
- Plant community composition: 201 → 111
- Plant carbon content: 649 → 29
- Plant nitrogen content: 142 → 31
- Output: `combined_edi_results_20260820_metadatacomplete_v3.csv`

---

## Final deliverables

- `extract_metadata.py` — extraction script
- `combined_edi_results_20260820_metadatacomplete_v3.csv` — final populated CSV (3,745 records processed, 2 empty-EML parse errors)

## Key design decisions

1. LTER site ecosystem/aquatic overrides applied via `LTER_SITE_INFO` map (keyed by 3-letter LTER code)
2. Word-boundary regex to avoid substring false positives
3. Institution names filtered from ecosystem/aquatic inference
4. Two-phase variable extraction: primary corpus first, methods text only as fallback
5. `claude_uncertainty` records notes for records with empty/inferable fields
