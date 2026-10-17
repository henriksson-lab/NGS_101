# Reference manifest

- Defining protocol: https://doi.org/10.1016/j.ymeth.2019.07.014
- Author-hosted full text: https://www.hubrecht.eu/app/uploads/2017/11/4C-seq-from-beginning-to-end_-A-detailed-protocol-for-sample-preparation-and-data-analysis.pdf
- Re-fetch pointer: `python3 tools/get_sources.py --doi 10.1016/j.ymeth.2019.07.014`

- Registered current pooling comparator: https://support-docs.illumina.com/SHARE/IndexAdaptersPooling/Content/TruSeqPoolingGuidelines.htm (v17; not assigned as the paper's historical revision).
- Appendix A supplementary data (inspected 2026-10-17; analysis figures only): https://ars.els-cdn.com/content/image/1-s2.0-S1046202318304742-mmc2.pdf (Supplemental figure 1) and https://ars.els-cdn.com/content/image/1-s2.0-S1046202318304742-mmc3.pdf (Supplemental figure 2). Not reached by `get_sources.py`; fetched by hand into the source directory.
- Pending source obligations: https://github.com/deLaatLab/pipe4C and https://github.com/deWitLab/peakC/. These are identified targets, not completed original reads.

Table 1 is raster content in PDF page 4 / printed page 20; inspect the original image, not only extracted text. Fetched files and conversion derivatives stay in `/husky/henriksson/NGS101_papers/sources/` and `_data/tmp/`. No third-party source file is committed.
