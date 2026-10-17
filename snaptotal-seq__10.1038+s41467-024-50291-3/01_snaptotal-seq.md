# snapTotal-seq — source intake

Defining paper: Niu, Luo & Zong, "Single-cell total-RNA profiling unveils regulatory
hubs of transcription factors", *Nature Communications* 15:5941 (2024), DOI
[10.1038/s41467-024-50291-3](https://doi.org/10.1038/s41467-024-50291-3); read as the
PMC11251146 JATS/PDF with the publisher supplements. The bioRxiv preprint
(10.1101/2023.12.14.571592) is superseded and not used for active claims.

Sources read (2026-10-17 bounded pass): JATS Methods "Cell lysis, reverse transcription,
and amplification by snapTotal-seq" and "Sequencing library construction" (text lines
246–270 of the `.xml.txt` twin), and Supplementary Data 12 (MOESM15 workbook, sheet
Table S12). Supplementary Fig. 1/2 (MOESM1) were not inspected in this pass.

## Workflow (Methods)

1. 🟢 Sorted single cells are lysed (42 °C, 3.5 min) in a buffer containing the primer
   mix (1.5 µM GAT27-N5-T12, 5 µM GAT27-N5-T3, 5 µM GAT27-N5-G3).
2. 🟢 Maxima H Minus RT runs 10 cycles of 8/15/20/30/42/50 °C followed by 50 °C 15 min;
   template switching uses the existing MALBAC primers as template-switching oligos
   (Results). RT is inactivated at 74 °C for 25 min.
3. 🟢 GoTaq Flexi PCR with GAT27, 23–26 cycles; 1.2X AMPure XP.
4. 🟢 Per-cell double-strand conversion (DSC) with the cell barcode primer, Deep Vent
   (exo-), 15 cycles 63/72 °C, EDTA stop. The authors add the barcode after PCR so that
   it lands on one cDNA end only.
5. 🟢 Pool 4 µL per cell, 1.1X AMPure, sonicate ~350 ng to 150–250 bp (Covaris S220),
   1.8X AMPure.
6. 🟢 Second DSC with Adapter-DSC (6 cycles 56/72 °C) "to enrich the DNA fragments with
   UMI and cell barcode"; 1.4X AMPure. 🟡 This is a primer-extension step, not a physical
   separation: untagged fragments are carried along and fall behind only in later PCRs.
7. 🟢 dA tailing with Klenow (3′→5′ exo-) and dATP, RT 30 min; 1.4X AMPure. No separate
   end-repair step is described (the earlier page drew one; removed).
8. 🟢 Quick Ligase Y-adapter ligation, RT 20 min, EDTA, 1.2X AMPure.
9. 🟢 PCR with Adapter-PCR + Illumina-P5, 10 cycles; 1.2X AMPure.
10. 🟢 DSN treatment (Evrogen) "to remove ribosomal reads": 100 ng, 95 °C 30 s, 80 °C
    3 h, DSN 80 °C 15 min, EDTA; 1.2X AMPure.
11. 🟢 Enrichment PCR with Adapter-DSC + Illumina-P5, 5 cycles; 1.2X AMPure.
12. 🟢 Final PCR, 5 cycles: Adapter-P7 + Illumina-P5 (custom sequencing), or Truseq-P7
    + Illumina-P5 (standard Illumina sequencing, 55 °C anneal).

## Published oligos (Supplementary Data 12, verbatim, spaces removed)

GAT27-5N-3G: `GTGAGTGATGGTTGAGGATGTGTGGAGNNNNNGGG`, GAT27-5N-3T: `GTGAGTGATGGTTGAGGATGTGTGGAGNNNNNTTT`, GAT27-5N-dT: `GTGAGTGATGGTTGAGGATGTGTGGAGNNNNNTTTTTTTTTTTT` 🟢 (NNNNN: UMI).

- GAT27 `GTGAGTGATGGTTGAGGATGTGTGGAG`
- cell barcode primer `GTGTGCTCTTCCGATCTNNNNNNNNAGGAGAGTGTGAGTGATGGTTGAGGATGTGTGGAG` (NNNNNNNN: cell barcode 1)
- Adapter-DSC `GCACGACATCTGCTAACGCAGTAGTGTGCTCTTCCGATCT`
- Truseq-P5-index primer `AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT` (NNNNNNNN: "cell barcode 2
  (library barcode)")
- Y-shape complementary sequence `P*GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTG` (`P*` as printed)
- Adapter-PCR `GCACGACATCTGCTAACGCAGTA`; Illumina-P5 `AATGATACGGCGACCACCGAGA`
- Adapter-P7 `CAAGCAGAAGACGGCATACGAGATGCACGACATCTGCTAACGCAGTA`
- Truseq-P7 `CAAGCAGAAGACGGCATACGAGATCGAGTAATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` (CGAGTAAT printed as an example index)
- custom Read 2 `AGAGGTGAGTGAGTGATGGTTGAGGATGTGTGGAG`; custom Index i7 `CTCCACACATCCTCAACCATCACTCACTCACCTCT`;
  Read 1 and Index i5 are the standard TruSeq primers.

## Final library and reads

🟡 The Y adapter ligates to both ends of every dA-tailed fragment; Adapter-PCR anneals
inside the Adapter-DSC handle, so the amplified molecule runs P5–i5–Read 1–insert–UMI–
spacer–cell barcode–handle. In the standard option, Truseq-P7 anneals to the cell-barcode
primer's 5′ `GTGTGCTCTTCCGATCT`, so the TruSeq Read 2 primer reads the 8-nt cell barcode,
the 35-nt primer-derived spacer and the 5-nt UMI. 🟢 This matches
the Methods: "Read 2 sequenced the cell barcode (bases 1–8) and UMI (bases 44–48)".
The earlier model placed the UMI next to the Read 2 primer; corrected.

🟢 In the custom run (NextSeq 500), Read 1 reads the captured RNA and Read 2 reads the
UMI; the analysis takes the UMI as the first five bases of Read 2.

## Open questions

- 🔴/🟡 Published discrepancy (shown on the page): the custom Read 2 primer's 5′ 8 nt
  (`AGAGGTGA`) differ from the cell barcode primer's `AGGAGAGT` before GAT27; the custom
  i7 primer (its reverse complement) therefore has 3′ mismatches against the drawn
  template. The source does not say which sequence is intended or how the custom index
  read primes; the custom-primer layout is not drawn.
- Acquired read lengths and the cell-barcode list are not given in the inspected passages.
- Supplementary Fig. 1/2 and the MALBAC/VASA-seq predecessor methods were not inspected.
