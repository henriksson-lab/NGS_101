# Source manifest

Sources live outside the checkout at `$CHEM_DATA/sources/t-scan__10.1016+j.cell.2019.07.009/` (default `/husky/henriksson/NGS101_papers/sources/`).

Re-obtain with `python3 tools/get_sources.py --slugs t-scan__10.1016+j.cell.2019.07.009`, then `python3 tools/scrape_primers.py "$CHEM_DATA/sources/t-scan__10.1016+j.cell.2019.07.009/"`. The bundle MANIFEST.tsv records exact URLs, validation and actual gated responses. Read scope and version qualifications are in the main note and `catalogue/source_reviews/t-scan__10.1016+j.cell.2019.07.009.json`. Original source material and extracted text are not committed.

Published versions take precedence; associated author manuscripts are explicitly qualified in the notes. Superseded preprints are excluded from current chemistry. Unreviewed files remain unreviewed and gated originals do not establish confirmed source absence.
