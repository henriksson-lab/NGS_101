# FLASH-seq oligo source review

Hahaut et al., “Fast and highly sensitive full-length single-cell RNA sequencing using
FLASH-seq,” *Nature Biotechnology* (2022),
[doi:10.1038/s41587-022-01312-3](https://doi.org/10.1038/s41587-022-01312-3).

Evidence: 🟢 printed in the defining paper's Supplementary Note 2 · 🟡 model-derived ·
🔴 unresolved.

## Sources read

- 🟢 `s41587-022-01312-3_supp_41587_2022_1312_MOESM1_ESM.pdf.txt`, Supplementary
  Note 2, “List of all the oligonucleotides used in this study” (under
  `$CHEM_DATA/sources/flash-seq__10.1038+s41587-022-01312-3/`).
- 🟢 The same supplement's FLASH-seq and FLASH-seq-with-UMI protocols distinguish the
  main chemistry from the tested optional UMI branch.

## Main FLASH-seq branch

🟢 The defining supplement prints these ordered oligos; `/5Biosg/` is a 5′ biotin and
`rG` denotes a ribonucleotide:

```text
FS dT30VN  /5Biosg/AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN
FS TSO    /5Biosg/AAGCAGTGGTATCAACGCAGAGTACrGrGrG
ISPCR     /5Biosg/AAGCAGTGGTATCAACGCAGAGT
```

The oligo-dT, TSO, and PCR primer share the same SMART handle. The TSO has an `AC`
linker and no UMI or SMART-seq3 tag. Its common-handle PCR and Nextera tagmentation
lead to a full-length, plate-indexed library. 🟡 The final Nextera read-primer sites
are computed from the common library construct.

## Separate FS-UMI branch

🟢 Supplementary Note 2 prints several tested UMI TSOs. The recommended `TSO_CTAAC`
is `/5Biosg/AAGCAGTGGTATCAACGCAGAGTNNNNNNNNCTAACrGrGrG`; it uses the SMART handle,
not the SMART-seq3 `AGAGACAGATTGCGCAATG` handle. The supplement also lists a
`Tn5_ISPCR_F` enrichment primer for this branch. Its distinct enrichment and read
layout need a separate construct; the public FLASH-seq page depicts the main branch.

## Open scope

- 🔴 Exact per-well i5/i7 sample-index assignments are not established here. The page
  shows their source-known Nextera primer architecture with index windows.
- 🔴 The optional FS-UMI branch is source-identified but has no separate page or complete
  primer-placement model yet.
