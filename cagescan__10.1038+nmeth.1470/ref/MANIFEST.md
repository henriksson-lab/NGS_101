# Source manifest

Retrieve originals without committing them:

```sh
python3 tools/get_sources.py --slugs cagescan__10.1038+nmeth.1470
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/cagescan__10.1038+nmeth.1470
```

- Defining article: https://doi.org/10.1038/nmeth.1470 (PMC2906222).
- Original Online Methods: `nmeth.1470_supp_NIHMS207200-supplement-supplement_1.pdf`.
- Original eight-oligo CAGEscan table: MOESM315 PDF p.26, Supplementary Table 8.
- Other original NIHMS supplements contain figures/result tables; their whole linked-source
  scope is not inferred from the inspected methods/table.

Downloaded originals, extracted text and figure derivatives stay under `$CHEM_DATA`.
Exact independent reading locations, source-unit notice and remaining obligations are
recorded in `catalogue/source_reviews/`.
