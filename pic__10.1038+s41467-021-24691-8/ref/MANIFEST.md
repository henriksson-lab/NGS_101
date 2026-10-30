# Source manifest

Originals and extracted text belong under `$CHEM_DATA/sources/pic__10.1038+s41467-021-24691-8/` (default `/husky/henriksson/NGS101_papers/sources/`), outside the checkout.

Run `python3 tools/get_sources.py --slugs pic__10.1038+s41467-021-24691-8` then `python3 tools/scrape_primers.py "$CHEM_DATA/sources/pic__10.1038+s41467-021-24691-8/"`. The external MANIFEST.tsv preserves actual URLs, duplicate identities and failure reasons. Main-note/source-review records distinguish browser-read primary originals from files successfully cached by the standard tool. Missing local material makes source-dependent checks skip; this is not evidence of missing disclosure. No fetched papers, images or text twins are committed.
