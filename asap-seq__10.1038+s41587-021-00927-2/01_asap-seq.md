# ASAP-seq — research notes

Defining source: Mimitou et al., *Nat. Biotechnol.* 39:1246–1258 (2021),
[doi:10.1038/s41587-021-00927-2](https://doi.org/10.1038/s41587-021-00927-2) (PMC8763625,
author manuscript). Published version used; the bioRxiv preprint 10.1101/2020.09.08.286914
is superseded and not used as authority.

## Sources read (bounded pass 2026-10-17, round breadth-first-06)

- 🟢 PMC8763625 JATS (`s41587-021-00927-2_PMC8763625.xml.txt`) Results line 122 (bridge
  principle, UBI), 124 (0.5 µmol BOA, pre-/post-SPRI tag input), 132 (BOA+BOB co-use);
  Methods lines 225–242 ("Cell fixation and permeabilization for ASAP-seq", "Transposition
  and barcoding for ASAP-seq"); Extended Data Fig. 1b/2a legends (lines 332, 338).
- 🟢 Publisher supplementary protocol MOESM1 (`…MOESM1_ESM.pdf.txt`), ASAP-seq protocol
  (Version September 2020): staining, fixation/lysis, 10x modifications, tag PCR, buffers, Oligos.
- 🟢 Historical 10x CG000209 Rev A (v1.1) guide, Appendix "Oligonucleotide Sequences" p48
  (`.pdf.txt` lines 3520–3625): Read 1N/2N, gel-bead primer PN-2000210, linear-amplification
  product, SI-PCR Primer B, Sample Index N primer and final ATAC library. Fetched in this pass.
- 🟢 MOESM3 Supplementary Table 1: TotalSeq antibody barcodes (all 15 nt in the barcode columns).
- Not inspected: the cited CG000209 **Rev D** and CG000168 Rev D (v1); the PMC author-manuscript
  supplements `ASAP_DOGMA_protocols.pdf`/`ASAP_DOGMA_supp_tables.xlsx` (PMC download gate);
  MOESM2 (no extractable text); BioLegend TotalSeq-A/-B conjugate documentation.

## Workflow (source order)

1. 🟢 Stain the live single-cell suspension: Fc block, antibody-oligo pool 30 min at 4 °C,
   three washes (MOESM1 §1). Staining precedes fixation.
2. 🟢 Fix in 1% formaldehyde 10 min RT, glycine quench, PBS washes; permeabilize on ice with
   LLL (mtDNA retained) or OMNI buffer; optional intracellular staining after lysis (MOESM1 §2).
3. 🟢 10x scATAC workflow (Rev D cited): transposition, then barcoding mix (step 2.1) with 0.5 µl
   1 µM bridge oligo. GEM incubation (step 2.5) gets a 40 °C 5 min annealing step before
   72 °C 5 min, 98 °C 30 s and 12 cycles (98 °C 10 s, 59 °C 30 s, 72 °C 1 min).
4. 🟢 The 3′-blocked bridge templates extension of the antibody tag in the first cycle, adding
   the Read 1N handle (and, for BOA, the N9V UBI); the bead primer then barcodes the extended
   tag with the chromatin fragments (Results line 122, Fig. 1b legend).
5. 🟢 Silane elution 43.5 µl, 3 µl kept for tag PCR; 1.2× SPRI-bound fraction continues as the ATAC
   library; supernatant raised to 2.0× SPRI (32 µl beads) and eluted as tag input.
6. 🟢 Tag PCR: KAPA 2×, P5 + RPI-x (TSA) or D7xx_s (TSB/hashtag), 95 °C 3 min, 14–16 cycles,
   expected product 190 bp (MOESM1). The P5 primer is named but not printed.

## Source-printed oligos

| Role | 5′→3′ as printed | Source |
|---|---|---|
| BOA, TotalSeq-A bridge | `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGNNNNNNNNNVTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT/3InvdT/` | Methods; MOESM1 |
| BOB, TotalSeq-B bridge | `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTTGCTAGGACCGGCCTTAAAGC/3InvdT/` | Methods; MOESM1 |
| RPI-x (TSA) | `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` | Methods |
| D7xx_s (TSB, TSA hashing) | `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCAGACGTGTGC` | Methods |
| Gel-bead primer PN-2000210 | `AATGATACGGCGACCACCGAGATCTACAC-N16-TCGTCGGCAGCGTC` | CG000209 Rev A p48 |
| SI-PCR Primer B | `AATGATACGGCGACCACCGAGA` | CG000209 Rev A p48 |
| Sample Index N | `CAAGCAGAAGACGGCATACGAGAT-N8-GTCTCGTGGGCTCGG` | CG000209 Rev A p48 |

🟢 Methods: the RPI-x 3′ part is the TruSeq Small RNA handle present in TSA tags; `x` is a
user-chosen sample index. 🟡 The BOA 5′ 33 nt equal Nextera Read 1 (computed in the model).

## Final libraries

- 🟢 ATAC: `P5–N16 cell barcode–Read 1N–insert–Read 2N rc–N8 i7–P7 rc` exactly as CG000209 Rev A
  prints it. The cell barcode is between P5 and Read 1N and is acquired by the i5 read; there is
  no i5 sample index. (Model correction: the earlier page placed an 8-nt i5 sample index at P5
  and put the cell barcode after the mosaic end.)
- 🟡 TotalSeq-A tag: `P5–N16–Read 1N–N9V UBI–T30–[antibody barcode rc]–RPI-x arm rc–N8–P7 rc`,
  derived from the bead primer, BOA and RPI-x. The barcode adjacency to the BOA poly(dT) is
  inferred because the conjugate sequence (poly(A) length, any linker) is not printed in the
  inspected sources. (Model correction: the earlier page drew a generic TruSeq Read 1 library
  with a 10-nt tag UMI.)
- TotalSeq-B (BOB) and hashtag libraries are not drawn.

## Read layout

🔴 Not stated in the inspected article sections: acquired read cycles, and the Read 2 / Index 1
primers for the tag library's small-RNA arm. The model assigns no cycle counts. (Model
correction: the earlier 50/50 and 28/50 lengths had no source.) CG000209 Rev A p42 recommends
ATAC cycles 50/8/16/50 (lines 3132–3146); that is a vendor recommendation, not the authors' run.

## Open questions

- Compare CG000209 Rev D and CG000168 Rev D appendices with Rev A.
- Read the PMC author-manuscript protocol supplement for any differences from MOESM1.
- BioLegend TotalSeq-A/-B conjugate structure; the exact P5 primer used for tag PCR.
- Missing-information claims remain blocked (no two-reader gate in this pass).
