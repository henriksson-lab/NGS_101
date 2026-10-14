# FIND-seq — research notes

Primary protocol: Shin *et al.*, *Nature Protocols* (2024), [10.1038/s41596-024-01021-y](https://doi.org/10.1038/s41596-024-01021-y) ([open author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC11537836/)). Defining study: Clark *et al.*, *Nature* (2023), [10.1038/s41586-022-05613-0](https://doi.org/10.1038/s41586-022-05613-0).

## 1. What it is

Focused interrogation of cells by nucleic-acid detection and sequencing (FIND-seq) physically sorts rare cells using RNA- or DNA-specific TaqMan fluorescence, while preserving their polyadenylated transcriptomes for downstream bulk-bin or single-cell RNA sequencing.

## Sources read

- 🟢 The protocol's development overview, three-device workflow, bead chemistry, TaqMan selection and transcriptome-recovery sections were read.
- 🟢 The defining Nature paper's FIND-seq Methods and workflow establish the original RNA-marker implementation.
- 🟢 The 2024 protocol's Materials → Oligonucleotides list gives the complete Smart-seq2 and Smart-seq3 capture, template-switching and amplification oligos, including their 5′ and terminal modifications.
- The source fetch commands are `python3 tools/get_sources.py --doi 10.1038/s41596-024-01021-y` and `python3 tools/get_sources.py --doi 10.1038/s41586-022-05613-0`.

## Construct-changing and selection steps

1. 🟢 Acrydite-modified oligo(dT) is covalently linked to allyl agarose with APS/TEMED.
2. 🟢 A bubble-triggered microfluidic device co-encapsulates single cells, molten agarose and strong Proteinase K/LiDS lysis buffer. Cooling forms a genome-capture hydrogel bead.
3. 🟢 Polyadenylated RNA hybridizes to the covalent oligo(dT) matrix and is reverse-transcribed; genomic DNA remains physically entrapped for marker interrogation.
4. 🟢 Beads are re-injected with TaqMan PCR reagents. Thermocycling generates fluorescence in droplets containing the selected RNA or DNA marker.
5. 🟢 Fluorescent droplets are sorted. The marker amplicon is a selection signal, not the transcriptome library.
6. 🟢 Selected matrix-bound cDNA is recovered and subjected to Smart-seq2 or, in the updated protocol, UMI-aware Smart-seq3 whole-transcriptome amplification and library preparation.

## Published oligos

- 🟢 **Smart-seq2 branch:** 5′-acrydite-TTTTTTAAGCAGTGGTATCAACGCAGAGTAC(T)₃₀VN-3′ capture primer; 5′-AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG-3′ template-switch oligo; 5′-AAGCAGTGGTATCAACGCAGAGT-3′ Smart PCR primer.
- 🟢 **Smart-seq3 branch:** 5′-biotin-AGAGACAGATTGCGCAATGNNNNNNNNrGrGrG-3′ template-switch oligo; 5′-acrydite-TTTTTTACGAGCATCAGCAGCATACGA(T)₃₀V-3′ capture primer; 5′-TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGATTGCGCAA*T*G-3′ forward primer and 5′-ACGAGCATCAGCAGCATAC*G*-3′ reverse primer. The source marks `*` as a phosphorothioate bond. The eight random bases in the Smart-seq3 TSO supply a UMI.
- 🟢 The same Materials list gives exact env and Ψ TaqMan primer/probe sets for the demonstrated viral-marker selection. They assay the bead before transcriptome library preparation and do not enter that library; other marker choices require their own primers and probes.

## Final library and read layout

- 🟢 The final transcriptome is a tagmented Smart-seq library; the exact UMI/read layout depends on whether Smart-seq2 or Smart-seq3 and bulk-bin or single-bead output is selected.
- 🟢 The public page uses the standard Nextera primer geometry for the displayed recovered-cDNA branch. Selection-target TaqMan primers do not become part of that library.

## Boundary

FIND-seq here means nucleic-acid cytometry for rare-cell transcriptomics, not the unrelated “full interrogation of nuclease DSBs” acronym.
