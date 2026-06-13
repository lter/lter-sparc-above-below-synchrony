"""
Download EML XML files from the PASTA REST API.

Reads the id column from combined_edi_result.csv and parses each value into
scope, identifier, and revision by splitting on the last two dots
(e.g. knb-lter-mcm.79.9 -> scope=knb-lter-mcm, identifier=79, revision=9).

Endpoint: GET https://pasta.lternet.edu/package/metadata/eml/{scope}/{identifier}/{revision}

Files are saved to eml/<id>.xml. Already-downloaded non-empty files are skipped.
Rate-limited to 0.25s between requests to stay under EDI's 5 req/sec limit.
"""
import csv
import os
import subprocess
import time

CSV_FILE = os.path.join("/Users/gmaurer/Library/CloudStorage/GoogleDrive-gmaurer@nmsu.edu/Shared drives/LTER-SPARC_Above-Below-Synchrony/01_DATA/01_Raw_data", "combined_edi_results_20260613_aiready.csv")
EML_DIR = os.path.join("/Users/gmaurer/Library/CloudStorage/GoogleDrive-gmaurer@nmsu.edu/Shared drives/LTER-SPARC_Above-Below-Synchrony/01_DATA/01_Raw_data", "result_eml")
BASE_URL = 'https://pasta.lternet.edu/package/metadata/eml'

os.makedirs(EML_DIR, exist_ok=True)

with open(CSV_FILE) as f:
    rows = list(csv.DictReader(f))

total = len(rows)
downloaded = 0
skipped = 0

for i, r in enumerate(rows, 1):
    pkg_id = r['repository_id']
    outfile = os.path.join(EML_DIR, f'{pkg_id}.xml')

    if os.path.isfile(outfile) and os.path.getsize(outfile) > 0:
        skipped += 1
        continue

    # Parse scope.identifier.revision (scope may contain dots, identifier and
    # revision are always the final two numeric segments)
    parts = pkg_id.rsplit('.', 2)
    if len(parts) != 3:
        print(f'[{i}/{total}] WARNING: cannot parse id "{pkg_id}", skipping')
        continue
    scope, identifier, revision = parts
    url = f'{BASE_URL}/{scope}/{identifier}/{revision}'

    print(f'[{i}/{total}] Downloading {pkg_id}...')
    subprocess.run(['wget', '-q', '--timeout=60', '--tries=2', '-O', outfile, url])
    time.sleep(0.25)
    downloaded += 1

print(f'\nDone. Downloaded: {downloaded}, Skipped (already present): {skipped}')
