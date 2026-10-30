# Source manifest

Sources live outside the checkout at `$CHEM_DATA/sources/tcr-map__10.1038+s41587-024-02248-6/` (default `/husky/henriksson/NGS101_papers/sources/`).

Re-obtain with `python3 tools/get_sources.py --slugs tcr-map__10.1038+s41587-024-02248-6`, then `python3 tools/scrape_primers.py "$CHEM_DATA/sources/tcr-map__10.1038+s41587-024-02248-6/"`. The bundle MANIFEST.tsv records exact URLs, validation and actual gated responses. Read scope and version qualifications are in the main note and `catalogue/source_reviews/tcr-map__10.1038+s41587-024-02248-6.json`. Original source material and extracted text are not committed.

Published versions take precedence; associated author manuscripts are explicitly qualified in the notes. Superseded preprints are excluded from current chemistry. Unreviewed files remain unreviewed and gated originals do not establish confirmed source absence.
