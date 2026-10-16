# Source manifest

- Defining article: https://doi.org/10.1016/j.cell.2022.04.003.
- Published original author copy: https://www.repository.cam.ac.uk/bitstreams/1d4e7717-3e50-4448-8bda-a35871dfbfef/download?format=pdf. The 38-page PDF includes Supplemental Figures S1–S7 on pp. 25–38; it does not make the separately cited Table S1 an inspected source.
- Original repository entry: https://www.repository.cam.ac.uk/items/bee9af47-0faa-4c0a-826e-fdabac4f63ef.
- Article-linked detailed procedure and oligo table: https://db.cngb.org/stomics/assets/html/stereo.seq.html.
- Preprint: https://doi.org/10.1101/2021.01.17.427004, v3; keep its TTGNAA precursor variant distinct from final Cell TTGAA.
- Recover via `python3 tools/get_sources.py --slugs stereo-seq__10.1016+j.cell.2022.04.003`; then run the full source oligo scraper before reading. Special routes are registered in `catalogue/extra_sources.tsv`; failed candidate routes are retrieval history, not proof of an absent original.
- Exact readings and unresolved Table S1/historical MGI manual identity are recorded in `catalogue/source_reviews/stereo-seq__10.1016+j.cell.2022.04.003.json`.

Third-party originals, extracted text and inspection images remain under ignored `_data/`; none are committed.
