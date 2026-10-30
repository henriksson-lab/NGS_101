# Source manifest

Originals and extracted text belong under `$CHEM_DATA/sources/hypr-seq__10.1073+pnas.2010738117/` (default `/husky/henriksson/NGS101_papers/sources/`), outside the checkout.

Run `python3 tools/get_sources.py --slugs hypr-seq__10.1073+pnas.2010738117` then `python3 tools/scrape_primers.py "$CHEM_DATA/sources/hypr-seq__10.1073+pnas.2010738117/"`. The external MANIFEST.tsv preserves actual URLs, duplicate identities and failure reasons. Main-note/source-review records distinguish browser-read primary originals from files successfully cached by the standard tool. Missing local material makes source-dependent checks skip; this is not evidence of missing disclosure. No fetched papers, images or text twins are committed.
