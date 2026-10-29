# NOMe-seq sources

All originals/text twins remain outside the checkout under `$CHEM_DATA/sources/nome-seq-bulk__10.1101+gr.143008.112/`.

```sh
python3 tools/get_sources.py --slugs nome-seq-bulk__10.1101+gr.143008.112
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/nome-seq-bulk__10.1101+gr.143008.112/
```

Published main PDF/XML are two representations of PMC3514679, one original. The Europe PMC supplementary zip contains main-article figure images. The actual publisher supplemental index is declared in `catalogue/extra_sources.tsv`; it returned403 and is a pending source, not evidence of omission. Exact reading sections and linked preparative targets are recorded in the note and source-review JSON. No source originals or extracted text are checked in.
