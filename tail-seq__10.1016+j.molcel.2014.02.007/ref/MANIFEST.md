# Source manifest

- Defining article: <https://doi.org/10.1016/j.molcel.2014.02.007>
- Re-obtain with `python3 tools/get_sources.py --doi 10.1016/j.molcel.2014.02.007`; publisher full text may require manual retrieval.

No source documents are committed.

## Round 2026-10-17-deferred-followup-02

- GEO sample text read in place of the gated article: `geo_GSM1242324_self_full.txt`, `geo_GSM1242325_self_full.txt` (from <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1242324&targ=self&form=text&view=full>, same for GSM1242325) and the series record `geo_GSE51299_self_brief.txt`, stored with MANIFEST.tsv rows in `$CHEM_DATA/sources/tail-seq__10.1016+j.molcel.2014.02.007/`.
- Still manual: article PDF <https://www.cell.com/article/S109727651400121X/pdf>, supplement mmc1 (cell.com) and the mirror <https://phmpt.org/wp-content/uploads/2024/06/125752_S2_M3_33_chang-2014.pdf> (all HTTP 403). SNU repository submitted version <https://hdl.handle.net/10371/171907> timed out.

## Round deferred-followup-03 (2026-10-17)

- Article PDF `j.molcel.2014.02.007.pdf` from <https://www.cell.com/article/S109727651400121X/pdf>: hand-download on an institutional network (fetcher gets HTTP 403).
- Supplements mmc1–mmc3 from the Elsevier CDN, e.g. <https://ars.els-cdn.com/content/image/1-s2.0-S109727651400121X-mmc1.pdf>. Hand-fetched on an institutional network.
- MANIFEST.tsv rows in the source bundle record the md5 of each file.
