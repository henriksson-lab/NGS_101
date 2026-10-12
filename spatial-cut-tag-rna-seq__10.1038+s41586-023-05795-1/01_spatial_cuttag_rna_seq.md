# spatial-CUT&Tag–RNA-seq

**Defining paper:** Zhang et al., *Nature* (2023), doi:10.1038/s41586-023-05795-1. The later step-by-step protocol is doi:10.1038/s41596-025-01145-9.

## Molecular boundary

🟢 An antibody selects a histone modification in fixed, permeabilized tissue and protein-A–Tn5 tags nearby DNA. A biotinylated poly(dT) primer independently captures RNA by in-situ reverse transcription.

🟢 Orthogonal microfluidic barcode-A and barcode-B ligations give both branches a shared tissue-pixel address before CUT&Tag and RNA products are recovered and separated.

## Evidence boundary

🟢 This remains separate from spatial-ATAC–RNA-seq: targeted pA–Tn5 replaces untargeted accessibility tagmentation even though the paper and downstream spatial barcoding are shared.

🟢 Supplementary Table 1 of the detailed protocol (Li et al., *Nature Protocols* 2025, doi:10.1038/s41596-025-01145-9) prints the complete 100-member A and B plates plus RT primer, Tn5ME-A/B/Rev, both ligation linkers and PCR/index primers. Tn5ME-B and the RT primer carry the same 15-nt `CATCGGCGTACGACT` ligation handle, allowing one exact barcode-A splint geometry on both branches. The second exact splint joins barcode B to the 5′-phosphorylated barcode-A arm.

🟢 Complete N501 and N701 sequences close the CUT&Tag branch without a role-level adapter placeholder: their 3′ arms bind Tn5ME-A and the barcode-B PCR arm, and their 5′ tails supply P5/i5 and P7/i7. The final representative panel is deliberately the CUT&Tag molecule. The separately recovered RNA product is Nextera-tagmented after pre-amplification, so no single source-defined cDNA cut position can stand for every RNA library molecule.
