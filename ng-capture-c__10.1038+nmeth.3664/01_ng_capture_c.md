# NG Capture-C — research notes

## Citation

Davies JOJ et al. *Multiplexed analysis of chromosome conformation at vastly improved sensitivity.* Nature Methods 13, 74–80 (2016). [doi:10.1038/nmeth.3664](https://doi.org/10.1038/nmeth.3664)

## Sources read

- 🟢 The defining article and Methods in [PMC4724891](https://pmc.ncbi.nlm.nih.gov/articles/PMC4724891/).
- 🟢 The later modular Capture-C protocol, [PMC7613269](https://pmc.ncbi.nlm.nih.gov/articles/PMC7613269/), used only to clarify the three modules and the distinction from later NuTi Capture-C.

## 1. What it is

🟢 NG Capture-C is a many-viewpoint-versus-all assay. A 3C proximity-ligation library is converted to an ordinary indexed Illumina library, then selected twice with biotinylated DNA oligos directed to the chosen restriction fragments.

## 2. Construct-changing steps

1. 🟢 Fix, restriction-digest and proximity-ligate chromatin to make 3C DNA.
2. 🟢 Sonicate the 3C library to about 200 bp. Random endpoints permit PCR duplicates to be recognized.
3. 🟢 Add paired-end Illumina sequencing adapters and amplify/index the library.
4. 🟢 Hybridize 120-nt biotinylated DNA capture oligos that overlap the restriction sites at both ends of each target fragment.
5. 🟢 Capture on streptavidin and repeat the capture. The defining paper reports that the sequential second capture is the key NG Capture-C enrichment change.

## 3. Final library and reads

🟢 Paired-end reads identify a viewpoint-containing sonication fragment and the reporter across its 3C ligation junction. The defining workflow does not use Hi-C-style biotin fill-in to select all ligation junctions.

🟢 The original Supplementary Data workbook (`nmeth3664_MOESM218_ESM.xlsx`, ‘Oligonucleotide sequences’ sheet, row 94) prints a 120-base Hba-1 viewpoint probe as 5′-`/5Biosg/GAT CCA TGG TAG CAC AGG GCA GCT AAG ATG CAA GTC TGA AGG AGG AGT CTG GCG AGC TGC TCC TGC AGT TCC CTG GAC CCA GAA GGA TGA GCT AGC AGA TTC ACT TGA GCC AAA GGA TTC`-3′. The workbook lists further viewpoint-specific probes and positions. The public page shows this one source-written example; the graph retains the capture-oligo role because each viewpoint uses a different sequence.

🟡 The schematic uses canonical Illumina-compatible primer sites.

## 4. Boundary with related protocols

🟢 Later Nuclear, Low-Input, Titrated and NuTi Capture-C variants are not silently folded into this page.
