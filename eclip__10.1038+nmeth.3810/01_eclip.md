# eCLIP

## 1. What it is

Enhanced CLIP maps direct RNA contacts of an immunoprecipitated RNA-binding protein and retains cDNAs that terminate at the crosslink.

## 2. Sources and evidence

- 🟢 Van Nostrand EL et al. *Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP).* Nature Methods (2016). [doi:10.1038/nmeth.3810](https://doi.org/10.1038/nmeth.3810).
- 🟢 ENCODE eCLIP SOP v1.P: https://www.encodeproject.org/documents/842f7424-5396-424a-a1a3-3f18707c3222/@@download/attachment/eCLIP_SOP_v1.P_110915.pdf
- 🟢 The source workflow UV-crosslinks RBP to RNA, partially digests RNA, immunoprecipitates the RNP, ligates a 3′ RNA adapter, reverse-transcribes, and ligates a second adapter to the cDNA 3′ end.
- 🟢 The cDNA-adapter ligation permits amplification of both RT-truncated and read-through molecules.
- 🟢 SOP v1.P publishes the X1A RNA adapter with its seven-base inline barcode and five random bases, rand103Tr3 with ten random bases, AR17, and the complete D5/D7 PCR-primer series. The model uses the X1A, D501 (`TATAGCCT`) and D701 (`CGAGTAAT`) examples.
- 🟢 The defining paper sequences 50-nt paired-end libraries. The PCR primers establish standard TruSeq Read 1, Index 1, Index 2 and Read 2 landing sites: this is dual-index, not single-end.

## 3. Molecular path

The crucial eCLIP innovation is adapter ligation after reverse transcription. An RT stop at the residual crosslinked peptide is therefore retained instead of being lost for lack of a second priming site.

Read 1 first reports X1A's seven-base barcode and five-base randomer. Read 2 first reports rand103Tr3's ten-base randomer, then approaches the RNA 5′ end or crosslink-stop boundary.
