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
- 🟢 The defining paper sequences paired-end 50-nt (N5 cDNA adapter) or 55-nt (N10) libraries (Methods). The PCR primers establish standard TruSeq Read 1, Index 1, Index 2 and Read 2 landing sites: this is dual-index, not single-end.
- 🟢 The paper's cited Supplementary Protocol 2 is in the second half of the same NIHMS PDF as Protocol 1 (PDF text lines 1289–1485). Its Cutadapt commands list adapter-trimming motifs and staggered substrings. Those strings document read processing; they are not additional ordered wet-lab oligos.
- 🟢 The [authors' eCLIP pipeline barcode reference](https://github.com/YeoLab/eclip/blob/master/example/inputs/yeolabbarcodes_20170101.fasta) lists X1A/X1B/X2A/X2B read strings as five random bases followed by the reverse complement of each seven-base inline barcode printed in Supplementary Protocol 1. This confirms read orientation for those tags. The lookup file does not identify which D5/D7 index pair was used for each paper sample.
- 🟢 The original X1A adapter ends in a 22-base `AGAUCGGAAGAGCGUCGUGUAG` tail (Supplementary Protocol 1, lines 329–333). `RiL19` is a separate 19-base adapter. The earlier model substituted an incorrect 21-base RiL19-like tail in the X1A ligation state; this was a repository transcription error, not a published discrepancy.
- 🟢 The original Supplementary Protocol prints the RNA adapters with `/5Phos/`, individual ribonucleotide notation and `/3SpC3/`; rand103Tr3 has `/5Phos/` and `/3SpC3/`. Those end modifications are reagent properties even though the DNA-only final library cannot display them.

## 3. Molecular path

The crucial eCLIP innovation is adapter ligation after reverse transcription. An RT stop at the residual crosslinked peptide is therefore retained instead of being lost for lack of a second priming site.

Read 1 first reports X1A's five-base randomer followed by `CCTATAT`, the reverse complement of its seven-base RNA-adapter barcode `ATATAGG`. Read 2 first reports rand103Tr3's ten-base randomer, then approaches the RNA 5′ end or crosslink-stop boundary.

## 4. Reaction order and conditions (2026-10-17 bounded pass)

Read: paper Online Methods (PMC4887338 JATS text line 123) and Figure 1a legend (line 363); Supplementary Protocol 1 (NIHMS765516 supplement-2 text lines 308–369 oligos, 575–620 on-bead 3′ ligation, 800–880 gel/membrane/proteinase K, 1046–1230 RT through PCR).

- 🟢 After stringent washes, RNA ends are treated on beads with FastAP and T4 PNK (pH 6.5 buffer), then the barcoded RNA adapter is ligated on beads with T4 RNA ligase 1 (high concentration) in PEG 8000/DMSO for 75 min at room temperature. Each sample receives two different barcoded adapters (SOP pairs such as X1A + X1B); the drawn panels use X1A.
- 🟢 Samples are run on a protein gel and transferred to nitrocellulose. The lane from the RBP band to 75 kDa above it is cut and proteinase K releases the RNA.
- 🟢 AR17 is annealed (65 °C) and AffinityScript reverse-transcribes at 55 °C for 45 min. ExoSAP-IT removes excess oligonucleotide, NaOH at 70 °C hydrolyses RNA (HCl neutralises), and MyONE Silane beads clean the cDNA.
- 🟢 rand103Tr3 (`/5Phos/N10…/3SpC3/`) is ligated by its 5′ phosphate to the cDNA 3′ end (Methods) with T4 RNA ligase 1 high concentration, PEG 8000 and DMSO overnight. The model previously drew the adapter at the cDNA 5′ end; this was a repository drawing error and is corrected.
- 🟢 qPCR on a 1:10 dilution sets the cycle number; Q5 PCR uses 6 cycles at 68 °C annealing plus a Ct-derived number of two-step cycles, followed by agarose gel size selection (Methods). The excised size range is not given in the inspected passages.
- 🟢 Sequencing was HiSeq 2500 or 4000 paired-end, 50 nt with the N5 cDNA adapter or 55 nt with the N10 adapter (Methods). The model draws rand103Tr3 N10 and therefore records 2 × 55; index-read cycle counts are not stated.

## 5. Published discrepancies

- The Online Methods give the 75 kDa membrane window as “~150 nt of RNA”; the Figure 1a legend gives “~220 nt of RNA”. The model uses only the 75 kDa cut, which both passages and the SOP state.
- Supplementary Protocol 1 lists the cDNA adapter as `rand103Tr3` (lines 352–353) but the ligation step is headed “5’ linker ligate cDNA” and adds “rand3Tr3 adapter” (lines 1126–1128). The model uses rand103Tr3 at the cDNA 3′ end as the Methods state; identity of the two names is inferred.

## 6. Open questions

- Which D5/D7 index pairs and which RNA adapter pairs were used for each published library (the author barcode lookup does not assign them).
- The N5 version of the cDNA adapter cited in the Methods is not printed in Supplementary Protocol 1; not modelled.
- Agarose size-selection range and acquired index-read cycles.
- Predecessor iCLIP (König 2010, Huppertz 2014) and RNAtag-Seq (Shishkin 2015) methods were not re-read in this pass.
