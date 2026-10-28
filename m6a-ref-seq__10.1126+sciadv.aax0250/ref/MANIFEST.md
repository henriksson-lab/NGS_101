# m6A-REF-seq sources

Originals and text twins stay in `$CHEM_DATA/sources/m6a-ref-seq__10.1126+sciadv.aax0250/`, outside this checkout.

```sh
python3 tools/get_sources.py --slugs m6a-ref-seq__10.1126+sciadv.aax0250
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/m6a-ref-seq__10.1126+sciadv.aax0250/
```

The manifest identifies published PMC6609220 PDF/XML (same original) and `aax0250_SM.pdf`. The note and source review record exactly inspected sections. The superseded preprint is excluded. Applicable historical VAHTS NR801 documentation is not identified; no current kit revision is substituted. No downloaded originals are checked in.
