# Source manifest

Re-obtain originals without committing them:

```sh
python3 tools/get_sources.py --slugs chro-seq__10.1038+s41588-018-0244-3
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/chro-seq__10.1038+s41588-018-0244-3
```

- Defining article: https://doi.org/10.1038/s41588-018-0244-3 (PMC6204104).
- Original Springer MOESM1 supplement, MOESM2 reporting summary and MOESM3 workbook;
  **Supplementary Table 1 is the workbook's first sheet**, not the reporting summary.
- Actual library predecessor: https://doi.org/10.1038/nprot.2016.086 (PMC5502525).
- Actual chromatin-isolation predecessor: Wuarin/Schibler 1994 (PMC359256).
- Original study GEO metadata GSE117832 and article-linked Danko-Lab PRO-seq pipeline.

Downloads and extracted text stay under `$CHEM_DATA/sources/`. Independent readings and
remaining source identities/link obligations are recorded in `catalogue/source_reviews/`.
