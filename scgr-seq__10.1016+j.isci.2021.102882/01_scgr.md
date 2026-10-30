# scGR-seq

Defining original: **Integrated analysis of glycan and RNA in single cells**, DOI[10.1016/j.isci.2021.102882](https://doi.org/10.1016/j.isci.2021.102882).

Capped source/model pass,2026-10-29. Private and incomplete;no accepted source-absence claim. Structural graph checks are separate from chemistry/source completeness.

## Selected originals

- `j.isci.2021.102882_PMC8349903.xml`: Published Results88/Fig1 and Methods330–350: conjugation/quench/desalt/click/affinity/dialysis, cell stain/washes/pick/UV, tag supernatant versus RNA pellet, tag PCR/AMPure/MiSeq and GenNextRamDA/NovaSeq branches.
- `j.isci.2021.102882_supp_mmc2_TableS1.xlsx`: Published TableS1 rPVL row: full reporter DNA and22-base lectin identity. Original file is a member of the published mmc2 TablesS1–S15 archive, not an independently guessed download.
- `j.isci.2021.102882_supp_mmc2_TableS15.xlsx`: Published TableS15 printed primer header, i5/i7 families and cell barcode description. Exported worksheet name TableS11 does not override printed TableS15/file identity.

## Chemistry and limits

🟢 Source Results specifies photocleavable DBCO-NHS attachment to lectin amino groups followed by coupling to5′-azide DNA;Methods abbreviates the reagent name as DBCO-NHS. Exact linker catalogue/structure and photo-product formula are unassembled;the model does not claim ordinary DBCO itself is photocleaved. Lectin activation20°C1h is followed by Tris quenching20°C15min and SephadexG25 removal of excess reagent. AzideDNA coupling at4°C precedes sugar–Sepharose affinity purification,dialysis and10kDa concentration. Quenching alone is not removal.

🟢 The representative rPVL reporter is `CGACGCTCTTCCGATCTCTGTCGCCCTGAACAACGGCAATGGCAGTCAGATCGGAAGAGCACAC`;its TableS1 lectin identity is `TCGCCCTGAACAACGGCAATGG`. The input DNA carries5′azide. After click attachment the original azide is not presented as an unreacted terminal group;after UV the local5′remnant is explicitly unassembled,not assigned phosphate or hydroxyl by assumption. This label is a lectin identity,not a UMI or cell identity.

🟢 Cells receive41DNA-barcoded proteins0.5µg/ml4°C1h,three PBS/BSA washes,and visual single-cell picking after counting. Source365nm15W15min exposure releases reporterDNA. Source Fig1/Results explicitly recovers DNA-bearing supernatant after centrifugation and retains the same-cell pellet forRNA. The two fractions can carry soluble or cell-associated remnants;neither is asserted chemically pure. No DNA barcode is transferred into RNA as part of this separation.

🟢 Tag supernatant uses NEBNextUltraIIQ5 and i5/i7 well-index primers for20cycles98/65,then AMPureXP cleanup and MultiNA assessment. TableS15 calls the8X windows cell barcodes;they are sequencing index windows that identify cell wells,not molecular UMIs. The model carries those meanings without assigning particular experimental cells. Source reports MiSeqv2 50cycle kit and26bp paired-end reads;both printed statements remain explicit without a silent25bp substitution. Exact acquired index cycles/read primers are unassembled.

🟢 Cell-pellet RNA continues through GenNextRamDA-seqSingleCellKit as directed by the source,then PowerUpSYBR/CFXlibraryqPCR and NovaSeq6000S4,151bp paired-end. Whole historical kit/cDNA/library geometry and exact manufacturer operations remain unread. Existing RamDA-seq code is not imported solely because the product name matches.

Component ancestry: published Methods350 explicitly directs the RNA branch to the GenNextRamDA kit and Hayashi etal.2018 source. That establishes a branch-specific source referral,not full recipe equivalence with an existing module. scGlycan-seq-only,bulkGlycan-seq and jointscGR-seq remain branches of this one defining technology;bulk-only cell number is not silently assigned to the single-cell graph.

## Deferred

- Published mmc1 figures including detailed FigS1 linker depiction; its cached existence does not mean inspected
- Direct GenNextRamDA kit historical guide/revision and Hayashi2018 predecessor; complete RNA preparation/library steps
- Exact photocleavable linker/released DNA end chemistry and full PCR/read geometry
- Remaining linked originals and independent source review

## Validation

Local source transcription and physical branch/state tests are recorded in the source-review JSON. No finished library or full source review is claimed.
