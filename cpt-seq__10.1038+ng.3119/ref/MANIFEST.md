# CPT-seq originals

Third-party originals and text twins belong in
`$CHEM_DATA/sources/cpt-seq__10.1038+ng.3119/`, normally under
`/husky/henriksson/NGS101_papers/`.

```sh
python3 tools/get_sources.py --slugs cpt-seq__10.1038+ng.3119
python3 tools/scrape_primers.py "$CHEM_DATA/sources/cpt-seq__10.1038+ng.3119/" --all
```

The external `MANIFEST.tsv` identifies source URLs/checksums. PMC4409979 is an
associated author manuscript. Publisher MOESM1 PDF and MOESM2 workbook carry
the20May2015 correction notice restoring the first five Table4 rows. The NIHMS
Word supplement remains gated; no complete equivalence to the publisher files
has been assumed. See the main note and source-review JSON for inspected scope.
