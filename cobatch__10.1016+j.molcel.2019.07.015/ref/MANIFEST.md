# Source manifest

- Defining article: <https://doi.org/10.1016/j.molcel.2019.07.015>
- Re-obtain article and supplementary oligo workbook with `python3 tools/get_sources.py --doi 10.1016/j.molcel.2019.07.015`.

No source documents are committed.

- Original publisher attachments: `https://ars.els-cdn.com/content/image/1-s2.0-S1097276519305453-mmc1.pdf` (Figures S1–S4), `mmc2.xlsx` (Table S1 oligos), `mmc3.xlsx` (Table S2 sequencing metadata), `mmc4.xlsx` (Table S3 comparison), `mmc5.xlsx` (Table S4 tissue digestion). All are declared in `catalogue/extra_sources.tsv`. The published article-plus-SI PDF (`mmc7.pdf`) and Methods S1 bench protocol (`mmc6.pdf`) at the same `ars.els-cdn.com/content/image/1-s2.0-S1097276519305453-` prefix were hand-fetched on 2026-10-17 and are declared in `catalogue/extra_sources.tsv` (the cell.com PDF route returns HTTP 403).
- Author-deposited GEO records (fetched by hand, no `extra_sources.tsv` row yet): `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129327&targ=self&form=text&view=quick` and `...acc=GSM3705255&targ=self&form=text&view=quick` (sample processing text with Read 1/Read 2 split cycles).
