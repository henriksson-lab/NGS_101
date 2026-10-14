# CRISPR-sciATAC — research notes

## Assay target

The sequenced material is human chromatin accessibility together with each cell's integrated perturbation guide cassette. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (barcoded TnY transposome ends, barcode plates and guide reverse-transcription and amplification primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 Liscovitch-Brauer et al. 2021, DOI `10.1038/s41587-021-00902-x`, preprint/full methods and supplements mirrored by the source fetcher.
- 🟢 Teichlab CRISPR-sciATAC page; CROP-seq and semi-suppressive-PCR papers were read only for inherited components.

## Construct path

🟢 Cells receive a lentiviral guide cassette, then accessible DNA is indexed by sci-ATAC-style combinatorial indexing. 🟢 A separate targeted amplification recovers the integrated guide identity and associates it with the same cellular index combination. The two outputs must remain separate libraries because the ATAC insert and guide-cassette amplicon have different molecular origins.
