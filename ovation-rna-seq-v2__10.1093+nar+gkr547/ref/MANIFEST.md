# Source manifest

- Tariq et al.: <https://doi.org/10.1093/nar/gkr547>
- NuGEN Ovation RNA-Seq System V2 user guide (commercial authority; locate by document number M01206)

Retrieve the paper with `python3 tools/get_sources.py --doi 10.1093/nar/gkr547`; source files are not committed.

The Supplementary Material is obtained from the Europe PMC supplementary zip
(`nar_gkr547_PMC3185437_supplementary.zip`, fetched by `get_sources.py`). Its
`supp_gkr547_Sup_Mat_Transcriptome_analysis_R2_SPI_clean.docx` member is a legacy
Word .doc despite the extension; extract it with `catdoc`. The Tecan route for the
V2 user guide returned HTTP 400 on 2026-10-17.
