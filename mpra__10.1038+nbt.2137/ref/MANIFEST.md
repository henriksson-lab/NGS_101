# MPRA reference originals

Defining published source: [10.1038/nbt.2137](https://doi.org/10.1038/nbt.2137).
Originals and text twins live under `$CHEM_DATA/sources/mpra__10.1038+nbt.2137/`.
The external `MANIFEST.tsv` records original URLs, hashes, duplicate copies and
failed access routes. No third-party source files are committed here.

Retrieve with the existing source tools:

```sh
python3 tools/get_sources.py --slugs mpra__10.1038+nbt.2137
python3 tools/scrape_primers.py "$CHEM_DATA/sources/mpra__10.1038+nbt.2137/"
```

Protocol-specific additional originals and author-copy routes are declared in
`catalogue/extra_sources.tsv`. Downloading a source does not establish that it
was read. Consult `../01_mpra.md` and the pooled-cloning curation ledger for
inspected sections, version exclusions and actual deferred original documents.
