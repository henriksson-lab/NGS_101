# CELLO-seq — research notes

## Sources read

- 🟢 Berrens et al. 2022, DOI `10.1038/s41587-021-01093-1`: Europe PMC author manuscript PMC7614850 (JATS `PMC7614850.xml.txt` Methods 'CELLO-seq' lines 176–181, 'ONT sequencing' line 207, sarlacc preprocessing line 219; PDF p35 Fig. 1A rendered) and Supplementary Table 5 (`EMS134943-supplement-Supplementary_Table_5.xlsx`, md5-identical to Springer MOESM6). The Springer Supplementary Information PDF (MOESM1) is image-only and was not read.
- 🟢 Berrens et al. 2025, DOI `10.1038/s41596-025-01203-2`, PMC detailed protocol and supplementary materials.

## Construct path

🟢 Poly(A)-primed RT (SuperScript IV with betaine and ET SSB; 10 min 57 °C, 120 min 42 °C) uses a 3′-amino-blocked TSO. Exonuclease I (37 °C 1 h) digests residual oligos, then 95 °C 3 min denatures enzyme and cDNA. Top and bottom splints are pre-annealed in hybridisation buffer and ligated with HiFi Taq DNA ligase plus NAD+ at 55 °C for 1 h, followed by Proteinase K, pooling and 0.5× AMPure XP. 🟢 PCR uses KAPA HiFi Uracil+ with betaine and ET SSB, a proxy qPCR to set 25–29 cycles (10 min extension) and 0.5× AMPure XP; 200 fmol are prepared with ONT SQK-LSK109 (hiPSC, 2-cell blastomeres) or SQK-LSK110 (mESC) and sequenced on MinION R9.4.1. It adds the cellular barcode and a 22-base UMI designed as repeated RYN classes to avoid long homopolymers. 🟢 High PCR duplicate counts are deliberate: noisy ONT reads sharing one long UMI are grouped to form a consensus suitable for locus-specific transposable-element assignment.

## Original oligos and version boundary

🟢 Supplementary Table 5 of the defining paper prints the original phosphorylated oligo(dT), amino-blocked TSO, top splint, and all 96 bottom splints. The page models the first listed bottom splint (`splintlig_spc_idx1`) explicitly. It carries, in order, the Nextera R2-derived tail, a 24-base cell barcode, `AGTGGTATC`, a 22-base patterned UMI (`NRY` repeated seven times plus `N`), the reverse complement of the top splint, and a terminal `T15` tract. The top splint has a 3′ inverted-dT block. The oligo(dT) is `Phos-` + 40-nt handle + `T25VN`.

🟡 Model correction (2026-10-17): the earlier model transcribed the oligo(dT) as `T30VN` and the bottom-splint tail as `T14`; Table S5 prints `T25VN` (67 nt after `Phos-`) and `T15` (162-nt idx1). These were repository transcription errors, not source discrepancies. Some bottom splints carry a 25-nt rather than 24-nt barcode (for example idx53, 163 nt); only idx1 is modelled.

🟢 The later 2025 Nature Protocols article changes this design: its bottom splint uses a 50-base UMI and the article explicitly describes that as replacing the original 22-base UMI. Those oligos define a revised implementation, not missing detail that can be copied into the 2022 construct. The public schematic therefore stays with the defining paper's 22-base design and links the revised protocol separately.

🟡 The published oligos determine the internal splint duplex and identifier-bearing product. The terminal Oxford Nanopore motor-loaded adapter is named by the workflow, but its bases and end geometry are not disclosed in the sources read, so that terminal part remains inferred rather than being invented.

## Open questions

- 🟡 Ligation geometry: the Methods say only that the annealed splint is ligated to first-strand cDNA, and Fig. 1A draws it at the 5′-phosphorylated oligo(dT) end. The bottom splint ends in the full oligo(dT) handle plus `T15` and the top splint is the handle's reverse complement plus an M13-like arm with a 3′ inverted dT; which strand templates the nick between the bottom-splint 3′ end and the cDNA 5′ phosphate is not described. The drawn junction is a model interpretation.
- 🟡 PCR primers: the Methods name only "forward and reverse PCR primers (IDT)" with "All oligo sequences ... in Table S5"; the NexteR1/NexteR2 pair matching the TSO and bottom-splint handles is inferred.
- Deferred: Springer Supplementary Information PDF (image-only), ONT SQK-LSK109/LSK110 guides (end-prep/dA-tailing and adapter geometry), the 2025 Nature Protocols detailed procedure beyond the version boundary.
