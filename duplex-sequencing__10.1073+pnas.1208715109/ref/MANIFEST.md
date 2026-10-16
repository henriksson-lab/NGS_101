# Source manifest

- Defining paper: https://doi.org/10.1073/pnas.1208715109, PMCID PMC3437896.
- Recover with `python3 tools/get_sources.py --slugs duplex-sequencing__10.1073+pnas.1208715109` and run the full oligo scraper before reading.
- The registered public original-PDF route is the USPTO exhibit URL in `catalogue/extra_sources.tsv`; its content is the original article plus SI, not patent commentary. It saves as `source_pdf_837ee6137323e8b633b0.pdf` with a doctext twin. Original SI begins on PDF page 7 and explicitly says **Supporting Information Corrected March 27, 2013**. The model is scoped to that version; pre-correction equivalence is unresolved.
- The PMC/PNAS direct SI route `1208715109_pnas.201208715SI.pdf` was gated in earlier retrieval attempts; that history does not make its chemistry absent now that an original corrected copy has been recovered.
- Later author protocol https://doi.org/10.1038/nprot.2014.170 (PMC4271547) publishes redesigned adapters with reversed tail polarity. It is a comparator, not a source for the corrected-original adapter bases.
- Exact original-source passes, source identity and remaining historical kit/version scope are recorded in `catalogue/source_reviews/duplex-sequencing__10.1073+pnas.1208715109.json`.

Third-party originals, text twins and rendered inspection images stay under ignored `_data/`.
