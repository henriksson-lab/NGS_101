# CHyMErA reference originals

Defining published source: [10.1038/s41587-020-0437-z](https://doi.org/10.1038/s41587-020-0437-z).
Originals and text twins live under `$CHEM_DATA/sources/chymera__10.1038+s41587-020-0437-z/`.
The external `MANIFEST.tsv` records original URLs, hashes, duplicate copies and
failed access routes. No third-party source files are committed here.

Retrieve with the existing source tools:

```sh
python3 tools/get_sources.py --slugs chymera__10.1038+s41587-020-0437-z
python3 tools/scrape_primers.py "$CHEM_DATA/sources/chymera__10.1038+s41587-020-0437-z/"
```

Protocol-specific additional originals and author-copy routes are declared in
`catalogue/extra_sources.tsv`. Downloading a source does not establish that it
was read. Consult `../01_chymera.md` and the pooled-cloning curation ledger for
inspected sections, version exclusions and actual deferred original documents.
