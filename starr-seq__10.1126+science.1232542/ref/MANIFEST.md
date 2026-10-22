# STARR-seq reference originals

Defining published source: [10.1126/science.1232542](https://doi.org/10.1126/science.1232542).
Originals and text twins live under `$CHEM_DATA/sources/starr-seq__10.1126+science.1232542/`.
The external `MANIFEST.tsv` records original URLs, hashes, duplicate copies and
failed access routes. No third-party source files are committed here.

Retrieve with the existing source tools:

```sh
python3 tools/get_sources.py --slugs starr-seq__10.1126+science.1232542
python3 tools/scrape_primers.py "$CHEM_DATA/sources/starr-seq__10.1126+science.1232542/"
```

Protocol-specific additional originals and author-copy routes are declared in
`catalogue/extra_sources.tsv`. Downloading a source does not establish that it
was read. Consult `../01_starr_seq.md` and the pooled-cloning curation ledger for
inspected sections, version exclusions and actual deferred original documents.
