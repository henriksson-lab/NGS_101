# PhIP-Seq originals

Third-party originals and extracted text belong outside the checkout at
`$CHEM_DATA/sources/phip-seq__10.1038+nbt.1856/` (default storage root:
`/husky/henriksson/NGS101_papers/`).

Re-obtain them with:

```sh
python3 tools/get_sources.py --slugs phip-seq__10.1038+nbt.1856
python3 tools/scrape_primers.py "$CHEM_DATA/sources/phip-seq__10.1038+nbt.1856/" --all
```

The external `MANIFEST.tsv` identifies each URL, checksum and provenance.
PMC4169279 is an associated author manuscript. Its NIHMS284271 supplement is
longer than the inspected13-page publisher supplement; neither filename nor the
shared DOI proves full revision equivalence. See the source-review record and
main note for exact inspected sections and remaining published-procedure checks.
