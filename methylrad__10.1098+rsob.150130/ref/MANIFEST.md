# Original sources

Run `python3 tools/get_sources.py --slugs methylrad__10.1098+rsob.150130`. All originals/text conversions belong in the external source store under that slug, not the checkout. The external `MANIFEST.tsv` records exact URLs, checksums and validation.

- [Published main, PMC4680569](https://pmc.ncbi.nlm.nih.gov/articles/PMC4680569/), DOI10.1098/rsob.150130. PMC Cloud `.1` JATS/PDF are the published article, not an author manuscript. Selected overview §2.1/Fig1 and Methods §§5.1–5.3.
- [Published supplementary PDF](https://pmc-oa-opendata.s3.amazonaws.com/PMC4680569.1/rsob150130supp1.pdf): Figures S1–S3 and Tables S1/S3/S5. Table S5's whole oligo rows, amine notation and index highlighting were checked in the rendered original, p6.
- Published `rsob150130supp2.xlsx` and `rsob150130supp3.xlsx` are fetched result workbooks. Scraper-only sequence pointers from the latter were not adopted as preparation facts; original workbook contents were not reviewed for chemistry.

Source-linked predecessor and applicable historical vendor originals remain ordinary follow-up targets, not failed-download or confirmed-absence claims. See the protocol note and `catalogue/source_reviews/` record for the exact bounded scope.
