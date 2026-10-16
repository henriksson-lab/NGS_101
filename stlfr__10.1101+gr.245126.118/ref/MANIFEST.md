# Source manifest

- Defining article: https://doi.org/10.1101/gr.245126.118 (PMC6499310).
- Recover with `python3 tools/get_sources.py --slugs stlfr__10.1101+gr.245126.118`, then run `python3 tools/scrape_primers.py _data/sources/stlfr__10.1101+gr.245126.118/` before reading.
- Detailed original: https://doi.org/10.1038/protex.2018.116; registered protocols.io import `x54v9rnqpv3e/v1`, printed integer 122070. The imported PDF identifies its original manuscript as attachment 5; attachment 6 is Supplement_0.pdf and attachment 7 Supplement_0.xlsx. Read the original attachment rather than assuming the import's internal references are unchanged.
- The oligo workbook is also recovered as bioRxiv `324392_324392-2.xlsx` (same content as protocol attachment 7).
- Actual linked precursors and author decoding repository have registered rows in `catalogue/extra_sources.tsv`: doi:10.1093/dnares/dsy037, doi:10.1126/science.1181498, https://github.com/stLFR/stLFR_read_demux.
- Downloads, extracted text and rendered inspection images remain in ignored `_data/`. Exact reading passes and remaining source scope are recorded in `catalogue/source_reviews/stlfr__10.1101+gr.245126.118.json`.
