# Source retrieval

Run `python3 tools/get_sources.py --slugs carlin__10.1016+j.cell.2020.04.048`, then the standard `tools/scrape_primers.py` on its external source directory.

Selected chemistry source is accepted-manuscript PMC7529102. Final publisher HTML returned HTTP 403; actual manuscript-era supplementary originals including Table S1 are gated. Exact URLs/availability are in the external `MANIFEST.tsv`. Primers in the fetched preprint SI are excluded, not a substitute final-version table. No sources or converted text belong in this checkout.
