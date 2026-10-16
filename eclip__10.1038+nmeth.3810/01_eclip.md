# eCLIP

## 1. What it is

Enhanced CLIP maps direct RNA contacts of an immunoprecipitated RNA-binding protein and retains cDNAs that terminate at the crosslink.

## 2. Sources and evidence

- 🟢 Van Nostrand EL et al. *Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP).* Nature Methods (2016). [doi:10.1038/nmeth.3810](https://doi.org/10.1038/nmeth.3810).
- 🟢 ENCODE eCLIP SOP v1.P: https://www.encodeproject.org/documents/842f7424-5396-424a-a1a3-3f18707c3222/@@download/attachment/eCLIP_SOP_v1.P_110915.pdf
- 🟢 The defining paper's own Supplementary Protocol 1 (PMC4887338, `NIHMS765516-supplement-2.pdf`, “Synthetic oligonucleotides”) independently prints X1A/X1B/X2A/X2B, RiL19, AR17, rand103Tr3, D501/D502 and D701/D702.
- 🟢 The source workflow UV-crosslinks RBP to RNA, partially digests RNA, immunoprecipitates the RNP, ligates a 3′ RNA adapter, reverse-transcribes, and ligates a second adapter to the cDNA 3′ end.
- 🟢 The cDNA-adapter ligation permits amplification of both RT-truncated and read-through molecules.
- 🟢 SOP v1.P publishes the X1A RNA adapter with its seven-base inline barcode and five random bases, rand103Tr3 with ten random bases, AR17, and the complete D5/D7 PCR-primer series. The model uses the X1A, D501 (`TATAGCCT`) and D701 (`CGAGTAAT`) examples.
- 🟢 The defining paper sequences 50-nt paired-end libraries. The PCR primers establish standard TruSeq Read 1, Index 1, Index 2 and Read 2 landing sites: this is dual-index, not single-end.
- 🟢 The paper's cited Supplementary Protocol 2 is in the second half of the same NIHMS PDF as Protocol 1 (PDF text lines 1289–1485). Its Cutadapt commands list adapter-trimming motifs and staggered substrings. Those strings document read processing; they are not additional ordered wet-lab oligos.
- 🟢 The [authors' eCLIP pipeline barcode reference](https://github.com/YeoLab/eclip/blob/master/example/inputs/yeolabbarcodes_20170101.fasta) lists X1A/X1B/X2A/X2B read strings as five random bases followed by the reverse complement of each seven-base inline barcode printed in Supplementary Protocol 1. This confirms read orientation for those tags. The lookup file does not identify which D5/D7 index pair was used for each paper sample.
- 🟢 The original X1A adapter ends in a 22-base `AGAUCGGAAGAGCGUCGUGUAG` tail (Supplementary Protocol 1, lines 329–333). `RiL19` is a separate 19-base adapter. The earlier model substituted an incorrect 21-base RiL19-like tail in the X1A ligation state; this was a repository transcription error, not a published discrepancy.
- 🟢 The original Supplementary Protocol prints the RNA adapters with `/5Phos/`, individual ribonucleotide notation and `/3SpC3/`; rand103Tr3 has `/5Phos/` and `/3SpC3/`. Those end modifications are reagent properties even though the DNA-only final library cannot display them.

## 3. Molecular path

The crucial eCLIP innovation is adapter ligation after reverse transcription. An RT stop at the residual crosslinked peptide is therefore retained instead of being lost for lack of a second priming site.

Read 1 first reports X1A's five-base randomer followed by `CCTATAT`, the reverse complement of its seven-base RNA-adapter barcode `ATATAGG`. Read 2 first reports rand103Tr3's ten-base randomer, then approaches the RNA 5′ end or crosslink-stop boundary.
