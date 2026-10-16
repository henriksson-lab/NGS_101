# Source manifest

- Llorens-Rico et al. 2025 and supplements: `python3 tools/get_sources.py --doi 10.1038/s41592-025-02700-8`.
- The original publisher/PMC-cloud `MOESM3` XLSX contains Supplementary Table 1a–d: oligo descriptions, general modified oligos, and separate 96-member ATAC and RNA sample-index families. Run `python3 tools/doctext.py` on that workbook to create its text twin; it is large and the automatic fetch does not convert it.
- `MOESM1` is the supplementary table legend and figures; `MOESM2` is an image-only reporting summary. `MOESM4` is an ATAC marker TSV ZIP, while `MOESM5` and `MOESM6` are result workbooks. Source-review details are in `catalogue/source_reviews/sum-seq__10.1038+s41592-025-02700-8.json`.
