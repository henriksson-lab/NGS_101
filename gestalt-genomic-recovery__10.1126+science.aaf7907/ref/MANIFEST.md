# Source retrieval

Run `python3 tools/get_sources.py --slugs gestalt-genomic-recovery__10.1126+science.aaf7907`, then the standard `tools/scrape_primers.py` on its external source directory.

External `MANIFEST.tsv` records exact original identities and failures. Selected original is accepted-manuscript PMC4967023. Both named manuscript supplements are gated; the declared institutional published-main PDF returned HTML rather than PDF. Superseded preprint originals are excluded. No source files belong in this checkout.
