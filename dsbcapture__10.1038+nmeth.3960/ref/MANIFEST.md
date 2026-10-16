# Source manifest

Run `python3 tools/get_sources.py dsbcapture`. The source routes include the defining
paper (doi:10.1038/nmeth.3960, PMCID PMC5045719), its PMC Cloud supplements, and the
original authors' detailed Protocol Exchange method (doi:10.1038/protex.2016.52).
Do not commit downloaded source files.

The article's Methods refer to Protocol Exchange with an unfilled “Insert URL”
placeholder. The recovered original method is `supplement_0.pdf` in the source bundle;
its extracted text prints the modified P5/P7 duplexes, the non-biotinylated control,
PCR F and the AD-indexed PCR R family on PDF pp. 8–9 (`supplement_0.pdf.txt`,
lines 608–653). It specifies a six-base TruSeq single index on PDF p. 7, line 574.
`tools/selftest_source_oligos.py dsbcapture__10.1038+nmeth.3960` checks exact graph
oligo states against this text when the source bundle is available.

The detailed method also cites Illumina's historical *TruSeq Library Prep Pooling
Guide* #15042173-01 (2015) for which barcodes to combine. Its availability and
relevance are tracked in `catalogue/source_reviews/`; the authors' own method
already prints the barcode sequences.
