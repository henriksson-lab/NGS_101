# wellDA-seq

**Defining source:** doi:10.1101/2024.09.06.611519.

## Molecular boundary

🟢 Indexed Tn5 accessibility products and a whole-genome branch are recovered from the same well-isolated cell.

## Evidence boundary

🟢 The open preprint and PMCID PMC11418942 establish the named method; the page does not infer one immutable commercial adapter sequence.

## Oligo source audit

🟢 The author preprint supplement (`2024.09.06.611519_media-1.pdf.txt`, Extended Data Fig. 1) prints the ATAC and DNA constructs at tagmentation, gap filling, co-enrichment PCR and modality enrichment PCR. Its ATAC branch uses the printed first-stage primer `AATGATACGGCGACCACCGAGATCTACACTGA-ATAC_CB1-GCCTCCCTCGCGCCAT` and second-stage primer `CAAGCAGAAGACGGCATACGAGAT-ATAC_CB2-CTTGCCAGCCCGCTCAG`; `ATAC_CB1` and `ATAC_CB2` are selected from Supplementary Table 5.

🟢 The DNA branch uses `AATGATACGGCGACCACCGAGATCTACAC-DNA_CB1-TCGTCGGCAGCGTC` and `CTGAGTCGGAGACACGCA-DNA_CB2-GTCTCGTGGGCTCGG` for co-enrichment, then introduces `DNA_CB3` in its modality-specific enrichment primer. Supplementary Table 6 gives the DNA barcode sequences. The figure prints the transposon and finished-library adapter handles in both strands, so the broad claim that the primer sequence is implementation-dependent is too strong for these demonstrated libraries.

🟡 The indexed windows vary by well. An exact page may choose one published barcode set or represent the windows symbolically while retaining the published flanking bases. The cited supplement is the preprint source; any later protocol variant should be kept distinct.
