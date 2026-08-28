"""
Download EML XML files from the PASTA REST API.

Reads the id column from combined_edi_result.csv and parses each value into
scope, identifier, and revision by splitting on the last two dots
(e.g. knb-lter-mcm.79.9 -> scope=knb-lter-mcm, identifier=79, revision=9).

Endpoint: GET https://pasta.lternet.edu/package/metadata/eml/{scope}/{identifier}/{revision}

Requests are authenticated via the EDI_KEY API key (appended as a ?key= query
parameter, per the ropensci/EDIutils add_api_key() mechanism).  The key is read
from the EDI_KEY environment variable, falling back to the project .Renviron file.

Files are saved to <SABRE_ROOT>/01_DATA/01_Raw_data/result_eml/<id>.xml.
Already-downloaded non-empty files are skipped.
Rate-limited to 0.25s between requests to stay under EDI's 5 req/sec limit.
"""
import csv
import os
import subprocess
import sys
import time
import urllib.parse

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RENVIRON_FILE = os.path.join(SCRIPT_DIR, '.Renviron')


def load_renviron():
    """Parse .Renviron and set variables in os.environ (preserves existing)."""
    if not os.path.isfile(RENVIRON_FILE):
        return
    with open(RENVIRON_FILE) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '=' not in line:
                continue
            name, _, value = line.partition('=')
            name = name.strip()
            if not name:
                continue
            value = value.strip().strip('"').strip("'")
            if name not in os.environ:
                os.environ[name] = value


load_renviron()
EDI_KEY = os.environ.get('EDI_KEY')
if not EDI_KEY:
    sys.exit('EDI_KEY not found. Set it in .Renviron or the environment.')

SABRE_ROOT = os.environ.get('SABRE_ROOT')
if not SABRE_ROOT:
    sys.exit('SABRE_ROOT not found. Set it in .Renviron or the environment.')

RAW_DATA_DIR = os.path.join(SABRE_ROOT, '01_DATA', '01_Raw_data')
CSV_FILE = os.path.join(RAW_DATA_DIR, 'combined_edi_results_20260820_aiready.csv')
EML_DIR = os.path.join(RAW_DATA_DIR, 'result_eml')
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
    url = f'{BASE_URL}/{scope}/{identifier}/{revision}?key={urllib.parse.quote(EDI_KEY, safe="")}'

    print(f'[{i}/{total}] Downloading {pkg_id}...')
    subprocess.run(['wget', '-q', '--timeout=60', '--tries=2', '-O', outfile, url])
    time.sleep(0.25)
    downloaded += 1

print(f'\nDone. Downloaded: {downloaded}, Skipped (already present): {skipped}')
