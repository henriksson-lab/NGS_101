# spatial-ATAC–RNA-seq

**Defining paper:** Zhang et al., *Nature* (2023), doi:10.1038/s41586-023-05795-1. The later step-by-step protocol is doi:10.1038/s41596-025-01145-9.

## Molecular boundary

🟢 Formaldehyde-fixed, permeabilized tissue is tagmented with Tn5 for accessible DNA while a biotinylated poly(dT) primer initiates in-situ reverse transcription.

🟢 Two orthogonal microfluidic channel sets deliver barcode A and barcode B. Their ordered pair is a two-dimensional tissue-pixel address shared by recovered ATAC and RNA molecules; the RNA construct includes a UMI.

## Evidence boundary

🟢 The defining article's extended-data chemistry diagrams distinguish the ATAC and RNA branches. Supplementary Table 1 of the detailed protocol (Li et al., *Nature Protocols* 2025, doi:10.1038/s41596-025-01145-9) prints the complete 100-member A and B plates plus RT primer, Tn5ME-A/B/Rev, both ligation linkers and PCR/index primers.

🟢 The RT primer and Tn5ME-B share the 15-nt `CATCGGCGTACGACT` handle. Linker 1 is its reverse complement followed by the reverse complement of the barcode-A ligating arm. Linker 2 analogously bridges the barcode-A outgoing arm to the barcode-B ligating arm. The schematic derives both splinted nicks from those exact sequences and shows A-1/B-1 as representative members of the published 100 × 100 spatial address.

🟢 The same table gives complete N501 and N701 primers. Their 3′ ends land on Tn5ME-A and the barcode-B PCR arm, respectively, so the ATAC branch's indexed P5-to-P7 shell and both PCR binding sites are derived rather than typed as generic adapter roles. The final representative panel is deliberately the ATAC molecule: the RNA product undergoes Nextera tagmentation after pre-amplification, so its insert-side cut position varies from molecule to molecule.
