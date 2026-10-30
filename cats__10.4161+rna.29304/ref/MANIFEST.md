# Source manifest

Source originals and extracted text live under `$CHEM_DATA/sources/cats__10.4161+rna.29304/` (default `/husky/henriksson/NGS101_papers/sources/`). Nothing is mirrored into the checkout.

Re-obtain with `python3 tools/get_sources.py --slugs cats__10.4161+rna.29304`, followed by `python3 tools/scrape_primers.py "$CHEM_DATA/sources/cats__10.4161+rna.29304/"`. The external MANIFEST.tsv preserves exact URLs, file identities and actual download/size limits. Main notes and source-review JSON qualify inspected versions and capped scope. Missing local originals make source-dependent checks skip; neither a failed fetch nor zero scraper hits establishes missing disclosure.
