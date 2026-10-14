# CapSeq reference material

Third-party papers, supplements and extracted text belong in
`$CHEM_DATA/sources/capseq__10.1126+science.ady7227/` and are not committed here.

- Defining paper: https://doi.org/10.1126/science.ady7227
- Author preprint: https://doi.org/10.1101/2025.03.14.642805
- Preprint PMC mirror and supplements: https://pmc.ncbi.nlm.nih.gov/articles/PMC11957016/
- Authors' analysis code: https://github.com/mazutislab/capseq
- Sequencing data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291189
- Submitter's library construction and read processing: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8830586
- Later inventors' patent (Example 3, Tables 1–5): https://eureka.patsnap.com/patent/WO2026190526A1
- Later CapMux read-layout paper: https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1846065/full

`python3 tools/get_sources.py --slugs capseq__10.1126+science.ady7227`
retrieved the bioRxiv PDF, its two supplements and the PMC mirror. The Science
publisher full text remains a manual row in the source manifest. Tables S21–S25
are cited by the preprint methods but absent from those two archived supplements;
revalidate a legitimately obtained Science oligo table in the same source directory.
