# Source manifest

Source originals and extracted text live under `$CHEM_DATA/sources/spacecat__10.1038+s41467-021-25279-y/` (default `/husky/henriksson/NGS101_papers/sources/`). Nothing is mirrored into the checkout.

Re-obtain with `python3 tools/get_sources.py --slugs spacecat__10.1038+s41467-021-25279-y`, followed by `python3 tools/scrape_primers.py "$CHEM_DATA/sources/spacecat__10.1038+s41467-021-25279-y/"`. The external MANIFEST.tsv preserves exact URLs, file identities and actual download/size limits. Main notes and source-review JSON qualify inspected versions and capped scope. Missing local originals make source-dependent checks skip; neither a failed fetch nor zero scraper hits establishes missing disclosure.
