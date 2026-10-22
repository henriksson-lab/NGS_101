# MoClo reference originals

Sources are stored outside the checkout in
`$CHEM_DATA/sources/moclo__10.1371+journal.pone.0016765/` (default storage root:
`/husky/henriksson/NGS101_papers`). The external `MANIFEST.tsv` records URLs, hashes,
download outcomes and duplicate originals. No third-party source files are committed.

- Defining article: [Weber et al. 2011](https://doi.org/10.1371/journal.pone.0016765),
  [PMC3041749](https://pmc.ncbi.nlm.nih.gov/articles/PMC3041749/).
- The `pmccloud:PMC3041749` row in `catalogue/extra_sources.tsv` retrieves the
  original article PDF/JATS and its published attachments through the existing
  PMC Cloud fetcher. Verify the article's own DOI before using its attachments.
- Table S1 is the original module inventory, `pone.0016765.s001.doc`.

After catalogue registration, retrieve with:

```sh
python3 tools/get_sources.py --slugs moclo__10.1371+journal.pone.0016765
python3 tools/scrape_primers.py "$CHEM_DATA/sources/moclo__10.1371+journal.pone.0016765/"
```

Retrieval is not source review. Inspected sections and deferred targets are recorded
in `01_moclo.md` and the pooled-cloning curation ledger.
