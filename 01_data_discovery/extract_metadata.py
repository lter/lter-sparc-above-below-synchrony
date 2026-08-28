#!/usr/bin/env python3
"""Extract metadata from all EML files and populate CSV. Batch production run."""

import csv
import os
import re
import xml.etree.ElementTree as ET
from datetime import date

CSV_PATH = "/home/ubuntu/agent-sandbox-io/sabre/combined_edi_results_20260820_aiready.csv"
EML_DIR = "/home/ubuntu/agent-sandbox-io/sabre/result_eml"
OUTPUT_PATH = "/home/ubuntu/agent-sandbox-io/sabre/combined_edi_results_20260820_metadatacomplete_v3.csv"

TODAY = "2026-08-26"

# LTER site info: code -> (ecosystem, is_aquatic)
LTER_SITE_INFO = {
    "knz": ("grassland", False),   # Konza Prairie
    "arc": ("tundra", False),       # Arctic
    "kbs": ("grassland", False),    # Kellogg Biological Station (agricultural)
    "fce": ("wetlands", True),      # Florida Coastal Everglades
    "hfr": ("forest", False),       # Harvard Forest
    "hbr": ("forest", False),       # Hubbard Brook
    "nwt": ("mountain", False),     # Niwot Ridge alpine
    "sev": ("desert", False),       # Sevilleta arid/semi-arid
    "and": ("forest", False),       # Andrews Forest
    "bes": ("forest", False),       # Baltimore Ecosystem Study (urban)
    "ble": ("tundra", False),       # Beaufort Lagoon
    "bnz": ("forest", False),       # Bonanza Creek boreal forest
    "cap": ("tundra", False),       # Central Alaska (permafrost/tundra)
    "cce": ("ocean", True),         # California Current Ecosystem
    "cdr": ("grassland", False),    # Cedar Creek grassland
    "cwt": ("forest", False),       # Coweeta forest
    "gce": ("wetlands", True),      # Georgia Coastal Ecosystems
    "jrn": ("desert", False),       # Jornada Basin desert
    "luq": ("forest", False),       # Luquillo tropical forest
    "mcm": ("desert", False),       # McMurdo Dry Valleys
    "mcr": ("ocean", True),         # Moorea Coral Reef
    "ntl": ("freshwater", True),    # North Temperate Lakes
    "pal": ("ocean", True),         # Palmer Antarctic
    "pie": ("wetlands", True),      # Plum Island Ecosystems
    "sbc": ("ocean", True),         # Santa Barbara Coastal
    "sgs": ("grassland", False),    # Shortgrass Steppe
    "vcr": ("wetlands", True),      # Virginia Coast Reserve
}

LTER_CODE_NAMES = {k: k.upper() for k in LTER_SITE_INFO}

# Patterns that should NOT trigger ecosystem inference (institution names, etc.)
FALSE_ECOSYSTEM_PATTERNS = [
    r'\bcenter\s+for\s+ocean', r'\bdepartment\s+of\s+marine',
    r'\bschool\s+of\s+marine', r'\binstitute\s+of\s+marine',
    r'\buniversity\s+of.*marine', r'\bgreat\s+lakes\b',
    r'\b36th\s+st\b.*marsh', r'\bmarsh\s+rd', r'\bmarsh\s+road',
    r'\bmarsh\s+lane', r'\bmarsh\s+street',
    r'\bmarsh\s+drive', r'\bmarsh\s+avenue',
]


def get_all_text(elem):
    if elem is None:
        return ""
    return "".join(elem.itertext()).strip()


def find_text(elem, tag):
    child = elem.find(tag)
    return child.text.strip() if child is not None and child.text else ""


def clean_text_for_ecosystem(text):
    """Remove institution names and addresses that falsely trigger ecosystem types."""
    cleaned = text
    for pat in FALSE_ECOSYSTEM_PATTERNS:
        cleaned = re.sub(pat, ' ', cleaned, flags=re.IGNORECASE)
    return cleaned


def infer_study_type(title, abstract, methods_text):
    """Infer study type. Prioritize title over methods."""
    # Check title/abstract first (higher confidence)
    title_abs = f"{title} {abstract}".upper()
    if re.search(r'GREENHOUSE\s+EXPERIMENT|GROWTH\s+CHAMBER|GREENHOUSE\s+PLOT|'
                 r'GROWN\s+IN\s+(?:POTS|GREENHOUSE)|GLASSHOUSE',
                 title_abs):
        return "greenhouse"
    if re.search(r'TRANSPLANT\s+EXPERIMENT|COMMON\s+GARDEN', title_abs):
        return "transplant"
    if re.search(r'\bMESOCOSM\b|\bMICROCOSM\b', title_abs):
        return "mesocosm"
    if re.search(r'LABORATORY\s+EXPERIMENT|LABORATORY\s+INCUBATION|'
                 r'LABORATORY\s+STUDY', title_abs):
        return "laboratory"

    # Then check methods
    mt = methods_text.upper()
    if re.search(r'GREENHOUSE\s+EXPERIMENT|GROWTH\s+CHAMBER|POTTED\s+IN|'
                 r'GROWN\s+IN\s+POTS|GREENHOUSE\s+TREATMENT', mt):
        return "greenhouse"
    if re.search(r'TRANSPLANT\s+EXPERIMENT|COMMON\s+GARDEN', mt):
        return "transplant"
    if re.search(r'\bMESOCOSM\b|\bMICROCOSM\b', mt):
        return "mesocosm"
    if re.search(r'LABORATORY\s+EXPERIMENT|LABORATORY\s+MICROCOSM', mt):
        return "laboratory"

    return "field"


def infer_obs_or_exp(all_text):
    exp_pattern = re.compile(
        r'\b(?:treatment|experiment|fertiliz(?:ed|ation)|'
        r'nutrient\s+addition|manipulation|control\s+plot|'
        r'herbivore\s+exclusion|burn\s+treatment|mowing\s+treatment|'
        r'nitrogen\s+addition|phosphorus\s+addition|experimental|'
        r'warming\s+treatment|transplant\s+experiment|'
        r'soil\s+amendment|removal\s+experiment|factorial|'
        r'greenhouse\s+treatment|heating|experimental\s+plots?|'
        r'treatment\s+effects|nutrient\s+enrichment|'
        r'climate\s+change\s+across\s+seasons)\b', re.IGNORECASE)
    return "experiment" if exp_pattern.search(all_text) else "observational"


def infer_ecosystem_type(repo_id, all_text, project_title):
    """Infer ecosystem type. Check LTER site info first."""
    repo_lower = repo_id.lower()

    # If it's a known LTER site, use that info
    for code, (eco, _) in LTER_SITE_INFO.items():
        if f"lter-{code}" in repo_lower or (
            repo_lower.startswith("knb-lter-") and code in repo_lower):
            return eco

    # Clean institution names from text
    cleaned = clean_text_for_ecosystem(all_text)

    if re.search(r'\btundra\b', cleaned, re.IGNORECASE):
        return "tundra"
    if re.search(r'\b(?:wetlands?|swamps?|marshes|bogs?|everglades?|fens?|sloughs?)\b',
                 cleaned, re.IGNORECASE):
        return "wetlands"
    if re.search(r'\b(?:alpine|subalpine|mountain)\b', cleaned, re.IGNORECASE):
        return "mountain"
    if re.search(r'\b(?:desert|arid|chihuahuan)\b', cleaned, re.IGNORECASE):
        return "desert"
    if re.search(r'\b(?:grassland|prairie|tallgrass|konza|pasture|meadow|steppe)\b',
                 cleaned, re.IGNORECASE):
        return "grassland"
    if re.search(r'\b(?:forest|woodland|oak\s+stand|hardwood|boreal|'
                 r'hubbard\s+brook|harvard\s+forest)\b', cleaned, re.IGNORECASE):
        return "forest"
    if re.search(r'\b(?:savannah|savanna)\b', cleaned, re.IGNORECASE):
        return "savannah"
    if re.search(r'\b(?:freshwater|lakes?|rivers?|streams?|ponds?|reservoir)\b',
                 cleaned, re.IGNORECASE):
        return "freshwater"
    if re.search(r'\b(?:ocean|marine|coastal|estuar|sea\b)\b', cleaned, re.IGNORECASE):
        return "ocean"
    if re.search(r'\b(?:agricultural|cropland|farmland|row\s+crop)\b',
                 cleaned, re.IGNORECASE):
        return "grassland"
    return "NA"


def infer_aquatic(repo_id, all_text):
    """Infer aquatic bool. Check LTER site info first."""
    repo_lower = repo_id.lower()

    # If it's a known LTER site, use that info
    for code, (_, is_aquatic) in LTER_SITE_INFO.items():
        if f"lter-{code}" in repo_lower or (
            repo_lower.startswith("knb-lter-") and code in repo_lower):
            if is_aquatic is not None:
                return "true" if is_aquatic else "false"

    # Clean text
    cleaned = clean_text_for_ecosystem(all_text)

    aquatic_re = re.compile(
        r'\b(?:aquatic|marine|ocean|freshwater|lakes?|rivers?|streams?|ponds?|'
        r'estuar(?:y|ine)|coastal|coral|reefs?|lagoon|'
        r'everglades?|sloughs?|salt\s+marsh|tidal|wetlands?|'
        r'swamps?|marshes|bogs?)\b', re.IGNORECASE)
    return "true" if aquatic_re.search(cleaned) else "false"


def extract_vars(dataset, title, methods_text):
    """Extract belowground and aboveground variables.
    Phase 1: match against primary corpus (title + entity descriptions + attributes).
    Phase 2: only if Phase 1 yields no results, fall back to methods text.
    """
    # Build primary corpus from high-confidence sources
    primary_parts = [title]
    for dt in dataset.findall("dataTable"):
        en = get_all_text(dt.find("entityName"))
        ed = get_all_text(dt.find("entityDescription"))
        if en:
            primary_parts.append(en)
        if ed:
            primary_parts.append(ed)
        al = dt.find("attributeList")
        if al is None:
            continue
        for attr in al.findall("attribute"):
            primary_parts.append(get_all_text(attr.find("attributeName")))
            primary_parts.append(get_all_text(attr.find("attributeDefinition")))

    primary_text = " ".join(primary_parts).lower()

    bg_vars = set()
    ag_vars = set()

    bg_rules = [
        ("Microbial biomass carbon",
         ["microbial biomass c", "microbial biomass carbon", "bmc", "microbial c",
          "biomass c", "microbial carbon"]),
        ("Microbial biomass nitrogen",
         ["microbial biomass n", "microbial biomass nitrogen", "bmn", "microbial n",
          "biomass n", "microbial nitrogen"]),
        ("Microbial biomass",
         ["microbial biomass"]),
        ("Soil respiration",
         ["soil respiration", "co2 flux", "co2 efflux", "soil co2"]),
        ("Soil inorganic nitrogen",
         ["inorganic n", "inorganic nitrogen", "soil inorganic nitrogen",
          "soil nitrate", "soil ammonium", "soil no3", "soil nh4",
          "extractable nitrogen", "extractable n"]),
        ("Soil water content",
         ["soil water content", "soil moisture",
          "gravimetric water", "soil water"]),
        ("Soil organic carbon",
         ["soil organic carbon", "soc", "total organic carbon",
          "soil organic matter", "som"]),
        ("Soil total nitrogen",
         ["total nitrogen", "total n", "soil nitrogen",
          "total kjeldahl", "total soil n"]),
        ("Microbial community composition",
         ["microbial community composition", "bacterial community",
          "fungal community", "plfa", "phospholipid fatty acid",
          "16s", "its1", "its2", "its region", "its gene",
          "internal transcribed spacer", "metagenom", "amplicon"]),
        ("Soil pH",
         ["soil ph", "ph of soil"]),
        ("Nematode abundance",
         ["nematode", "nematod"]),
        ("Root biomass",
         ["root biomass", "root mass", "belowground biomass",
          "fine root"]),
        ("Soil bulk density",
         ["bulk density", "soil bulk"]),
        ("Soil carbon flush",
         ["carbon flush", "cflush", "co2 flush"]),
        ("Soil nitrogen flush",
         ["nitrogen flush", "nflush"]),
        ("Microbial biomass C:N ratio",
         ["microbial biomass c:n", "microbial c:n", "c:n biomass",
          "microbial biomass ratio", "microbial c n"]),
        ("Mycorrhizal colonization",
         ["mycorrhizal", "mycorrhiza", "amf", "arbuscular"]),
        ("Soil enzyme activity",
         ["enzyme activity", "extracellular enzyme",
          "phosphatase", "glucosidase", "cellulase"]),
        ("Soil carbon stock",
         ["carbon stock", "c stock", "soil c stock"]),
        ("Dissolved organic carbon",
         ["dissolved organic carbon", "dissolved organic c", "doc "]),
        ("Soil food web biomass",
         ["food web", "foodweb", "microarthropod", "protozoa"]),
        ("Nitrogen mineralization",
         ["nitrogen mineralization", "n mineralization",
          "net mineralization"]),
        ("Microbial activity",
         ["microbial activity"]),
        ("Soil total carbon",
         ["soil total carbon", "total c", "soil carbon",
          "total soil carbon", "soil c"]),
        ("Microbial mat biomass",
         ["microbial mat", "microalgal"]),
    ]

    ag_rules = [
        ("Plant biomass",
         ["aboveground biomass", "plant biomass", "shoot biomass",
          "standing biomass"]),
        ("Vegetation cover",
         ["vegetation cover", "plant cover", "percent cover",
          "cover %", "canopy cover"]),
        ("Plant community composition",
         ["plant community composition", "plant species composition",
          "vegetation community"]),
        ("Plant production",
         ["productivity", "npp", "gpp", "net primary", "plant production",
          "aboveground production"]),
        ("Plant survival",
         ["plant survival", "seedling survival", "mortality rate"]),
        ("Plant height",
         ["plant height", "canopy height", "stem height"]),
        ("Plant density",
         ["plant density", "stem density", "seedling density"]),
        ("Leaf area index",
         ["leaf area index", "leaf area", "specific leaf area", "sla"]),
        ("Litter biomass",
         ["litterfall", "litter mass", "litter biomass", "leaf litter"]),
        ("Arthropod abundance",
         ["arthropod", "sweepnet", "aboveground arthropod"]),
        ("Plant carbon content",
         ["plant carbon", "foliar carbon", "leaf carbon"]),
        ("Plant nitrogen content",
         ["plant nitrogen", "foliar nitrogen", "leaf nitrogen"]),
        ("Seed germination",
         ["germination", "seedling emergence"]),
        ("Macrophyte biomass",
         ["macrophyte", "aquatic plant", "submerged vegetation"]),
        ("Reproductive output",
         ["reproductive output", "flower count", "fruit count",
          "seed production"]),
    ]

    # Phase 1: match against primary corpus
    for var_name, patterns in bg_rules:
        if any(p in primary_text for p in patterns):
            bg_vars.add(var_name)
    for var_name, patterns in ag_rules:
        if any(p in primary_text for p in patterns):
            ag_vars.add(var_name)

    # Phase 2: only if Phase 1 found nothing, fall back to methods + primary
    if not bg_vars and not ag_vars:
        fallback_text = f"{primary_text} {methods_text}".lower()
        for var_name, patterns in bg_rules:
            if any(p in fallback_text for p in patterns):
                bg_vars.add(var_name)
        for var_name, patterns in ag_rules:
            if any(p in fallback_text for p in patterns):
                ag_vars.add(var_name)

    return ("; ".join(sorted(bg_vars)) if bg_vars else "NA",
            "; ".join(sorted(ag_vars)) if ag_vars else "NA")


def extract_soil_depth(dataset, methods_text):
    """Extract soil depth."""
    combined = methods_text.lower()
    for dt in dataset.findall("dataTable"):
        al = dt.find("attributeList")
        if al is None:
            continue
        for attr in al.findall("attribute"):
            an = get_all_text(attr.find("attributeName"))
            ad = get_all_text(attr.find("attributeDefinition"))
            combined += f" {an} {ad}".lower()

    m = re.search(r"(\d+(?:\.\d+)?)\s*(?:-|—|to)\s*(\d+(?:\.\d+)?)\s*cm",
                   combined, re.IGNORECASE)
    if m:
        return f"{m.group(1)}-{m.group(2)} cm"

    m = re.search(r"(?:depth|sampled at)\s*(?:of\s*)?(\d+(?:\.\d+)?)\s*cm",
                   combined, re.IGNORECASE)
    if m:
        return f"0-{m.group(1)} cm"

    for dt in dataset.findall("dataTable"):
        al = dt.find("attributeList")
        if al is None:
            continue
        for attr in al.findall("attribute"):
            if get_all_text(attr.find("attributeName")).lower() == "depth":
                return "See attribute: depth"

    return "NA"


def extract_sampling_freq(all_text, year_start, year_end):
    """Infer sampling frequency."""
    try:
        dur = int(year_end) - int(year_start)
    except (ValueError, TypeError):
        dur = 0

    mt = all_text.lower()
    sub_annual_kw = ["weekly", "monthly", "biweekly", "bi-weekly",
                     "sampled every", "collected every",
                     "multiple times per year", "three times a year",
                     "twice a year", "quarterly", "every month",
                     "every 2 weeks", "bi-monthly", "bimonthly",
                     "several times", "times per year",
                     "sampling dates", "sampled in june and august",
                     "sampled june and", "sampled in spring and",
                     "sampled in summer and", "sampled in fall and",
                     "pre-burn", "post-burn"]

    if any(kw in mt for kw in sub_annual_kw):
        return "sub-annual"
    if dur >= 1:
        return "annual"
    return "other"


def extract_site_name(repo_id, dataset):
    """Extract site name. Check LTER codes first, then geog description."""
    repo_lower = repo_id.lower()
    for code, name in LTER_CODE_NAMES.items():
        if f"lter-{code}" in repo_lower or (
            repo_lower.startswith("knb-lter-") and code in repo_lower):
            return name

    for cov in dataset.findall("coverage"):
        for geo in cov.findall("geographicCoverage"):
            gd = get_all_text(geo.find("geographicDescription"))
            if gd:
                short = gd.split(":")[0].split(",")[0].strip()
                if len(short) < 60 and short and not short.startswith("http"):
                    return short
    return "NA"


def extract_license(dataset):
    """Extract license."""
    ir = dataset.find("intellectualRights")
    if ir is None:
        return "unlicensed"
    text = get_all_text(ir).lower()
    if "cc0" in text or "public domain" in text:
        return "CC0"
    if "cc-by" in text or "cc by" in text or "creative commons attribution" in text:
        return "CC-BY"
    return "other"


def parse_eml(filepath):
    """Parse EML XML file and extract all metadata."""
    try:
        tree = ET.parse(filepath)
        root = tree.getroot()
        dataset = root.find("dataset")
        if dataset is None:
            return None
    except Exception:
        return None

    result = {}

    repo_id = os.path.basename(filepath).replace(".xml", "")
    title = get_all_text(dataset.find("title"))
    abstract = get_all_text(dataset.find("abstract"))

    # Keywords
    kw_text = ""
    for ks in dataset.findall("keywordSet"):
        for kw in ks.findall("keyword"):
            kw_text += " " + get_all_text(kw)

    # Methods
    methods_text = ""
    for ms in dataset.findall("methods"):
        for mstep in ms.findall("methodStep"):
            methods_text += " " + get_all_text(mstep)

    # Project
    project_title = get_all_text(dataset.find("project/title"))

    result["study_type"] = infer_study_type(title, abstract, methods_text)
    result["obs_or_exp"] = infer_obs_or_exp(
        f"{title} {abstract} {methods_text} {project_title}".lower())

    # Geographic coverage
    geo_list = []
    for cov in dataset.findall("coverage"):
        for g in cov.findall("geographicCoverage"):
            geo_list.append(g)
    if geo_list:
        g = geo_list[0]
        result["bb_west"] = find_text(g, "boundingCoordinates/westBoundingCoordinate") or "NA"
        result["bb_east"] = find_text(g, "boundingCoordinates/eastBoundingCoordinate") or "NA"
        result["bb_north"] = find_text(g, "boundingCoordinates/northBoundingCoordinate") or "NA"
        result["bb_south"] = find_text(g, "boundingCoordinates/southBoundingCoordinate") or "NA"
    else:
        result["bb_west"] = result["bb_east"] = result["bb_north"] = result["bb_south"] = "NA"

    # Temporal coverage
    for cov in dataset.findall("coverage"):
        tc = cov.find("temporalCoverage")
        if tc is None:
            continue
        rod = tc.find("rangeOfDates")
        if rod is not None:
            begin = find_text(rod, "beginDate/calendarDate")
            end = find_text(rod, "endDate/calendarDate")
            result["year_start"] = begin[:4] if len(begin) >= 4 else (begin or "NA")
            result["year_end"] = end[:4] if len(end) >= 4 else (end or "NA")
            break
        sdt = tc.find("singleDateTime/calendarDate")
        if sdt is not None and sdt.text:
            y = sdt.text.strip()[:4]
            result["year_start"] = y
            result["year_end"] = y
            break
    else:
        result["year_start"] = result["year_end"] = "NA"

    # All text for inference
    all_text = f"{title} {abstract} {kw_text} {methods_text} {project_title}"

    # Site type
    if "knb-lter-" in repo_id.lower():
        result["site_type"] = "LTER"
    elif re.search(r'\b(?:LTER|Long-Term Ecological Research)\b',
                   f"{title} {project_title}", re.IGNORECASE):
        result["site_type"] = "LTER"
    elif re.search(r'\bNEON\b', f"{title} {project_title}", re.IGNORECASE):
        result["site_type"] = "NEON"
    else:
        result["site_type"] = "NA"

    result["site_name"] = extract_site_name(repo_id, dataset)

    # Ecosystem and aquatic - use all text for non-LTER, LTER overrides
    result["ecosystem_type"] = infer_ecosystem_type(repo_id, all_text, project_title)
    result["aquatic_bool"] = infer_aquatic(repo_id, all_text)

    # Variables
    bg, ag = extract_vars(dataset, title, methods_text)
    result["belowground_vars"] = bg
    result["aboveground_vars"] = ag
    result["soil_depth"] = extract_soil_depth(dataset, methods_text)

    # Sampling
    result["sampling_freq"] = extract_sampling_freq(all_text,
                                                     result["year_start"],
                                                     result["year_end"])
    try:
        dur = int(result["year_end"]) - int(result["year_start"])
    except (ValueError, KeyError):
        dur = 0
    result["sub_1_year_bool"] = "true" if dur < 1 and result["year_start"] != "NA" else "false"

    result["license"] = extract_license(dataset)

    # Uncertainty notes
    uncertainties = []
    if result["ecosystem_type"] == "NA":
        uncertainties.append("ecosystem_type uncertain")
    if result["study_type"] == "other":
        uncertainties.append("study_type uncertain")
    if result["belowground_vars"] == "NA" and result["aboveground_vars"] == "NA":
        uncertainties.append("no variables extracted")
    if result["site_name"] == "NA":
        uncertainties.append("site_name uncertain")
    if result["soil_depth"] == "NA":
        uncertainties.append("soil_depth uncertain")
    result["claude_uncertainty"] = "; ".join(uncertainties)

    return result


def main():
    # Read CSV
    rows = []
    with open(CSV_PATH, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            rows.append(row)

    total = len(rows)
    print(f"Processing {total} records...")

    stats = {
        "processed": 0,
        "missing_eml": 0,
        "parse_error": 0,
        "study_type": {},
        "ecosystem_type": {},
        "aquatic_true": 0,
        "aquatic_false": 0,
        "obs_or_exp": {},
        "site_type": {},
        "with_uncertainty": 0,
    }

    updated = []
    for i, row in enumerate(rows):
        repo_id = row["repository_id"]
        eml_path = os.path.join(EML_DIR, f"{repo_id}.xml")

        if not os.path.exists(eml_path):
            row["claude_uncertainty"] = "EML file not found"
            stats["missing_eml"] += 1
            updated.append(row)
            continue

        meta = parse_eml(eml_path)
        if meta is None:
            row["claude_uncertainty"] = "EML parse error"
            stats["parse_error"] += 1
            updated.append(row)
            continue

        # Fill row from meta
        for key, value in meta.items():
            if key in row:
                if row[key] == "NA" or not row[key]:
                    row[key] = value

        row["claude_eval_date"] = TODAY
        updated.append(row)
        stats["processed"] += 1

        # Aggregate stats
        for key in ["study_type", "ecosystem_type", "obs_or_exp", "site_type"]:
            val = meta.get(key, "NA")
            stats[key][val] = stats[key].get(val, 0) + 1

        if meta.get("aquatic_bool") == "true":
            stats["aquatic_true"] += 1
        else:
            stats["aquatic_false"] += 1

        if row.get("claude_uncertainty"):
            stats["with_uncertainty"] += 1

        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{total} records processed...")

    # Write output
    with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(updated)

    # Print summary
    print(f"\n{'='*60}")
    print(f"COMPLETE: {total} records")
    print(f"  Successfully processed: {stats['processed']}")
    print(f"  Missing EML files:     {stats['missing_eml']}")
    print(f"  Parse errors:          {stats['parse_error']}")
    print(f"  Records with uncertainty notes: {stats['with_uncertainty']}")
    print(f"\nstudy_type distribution:")
    for k, v in sorted(stats["study_type"].items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")
    print(f"\necosystem_type distribution:")
    for k, v in sorted(stats["ecosystem_type"].items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")
    print(f"\nobs_or_exp distribution:")
    for k, v in sorted(stats["obs_or_exp"].items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")
    print(f"\nsite_type distribution:")
    for k, v in sorted(stats["site_type"].items(), key=lambda x: -x[1]):
        print(f"  {k}: {v}")
    print(f"\naquatic_bool: true={stats['aquatic_true']} false={stats['aquatic_false']}")
    print(f"\nOutput: {OUTPUT_PATH}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
