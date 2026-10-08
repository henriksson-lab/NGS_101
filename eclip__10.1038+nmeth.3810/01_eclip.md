# eCLIP

## 1. What it is

Enhanced CLIP maps direct RNA contacts of an immunoprecipitated RNA-binding protein and retains cDNAs that terminate at the crosslink.

## 2. Sources and evidence

- 🟢 Van Nostrand EL et al. *Robust transcriptome-wide discovery of RNA-binding protein binding sites with enhanced CLIP (eCLIP).* Nature Methods (2016). [doi:10.1038/nmeth.3810](https://doi.org/10.1038/nmeth.3810).
- 🟢 ENCODE eCLIP experimental protocol: https://www.encodeproject.org/documents/fa2a3246-6039-46ba-b960-17ef06e7876a/
- 🟢 The source workflow UV-crosslinks RBP to RNA, partially digests RNA, immunoprecipitates the RNP, ligates a 3′ RNA adapter, reverse-transcribes, and ligates a second adapter to the cDNA 3′ end.
- 🟢 The cDNA-adapter ligation permits amplification of both RT-truncated and read-through molecules.
- 🟡 Barcode identities are experiment-specific; the page shows their roles and inferred canonical run geometry.

## 3. Molecular path

The crucial eCLIP innovation is adapter ligation after reverse transcription. An RT stop at the residual crosslinked peptide is therefore retained instead of being lost for lack of a second priming site.
