# DOGMA-seq — research notes

Defining source: Mimitou et al., *Nat. Biotechnol.* 39:1246–1258 (2021),
[doi:10.1038/s41587-021-00927-2](https://doi.org/10.1038/s41587-021-00927-2) (PMC8763625,
author manuscript). Published version used; the bioRxiv preprint 10.1101/2020.09.08.286914 is
superseded and not used as authority.

## Sources read (bounded pass 2026-10-17, round breadth-first-06)

- 🟢 PMC8763625 JATS (`s41587-021-00927-2_PMC8763625.xml.txt`) Methods lines 244–262 ("Cell
  permeabilization for DOGMA-seq", "Transposition and barcoding for DOGMA-seq").
- 🟢 Publisher supplementary protocol MOESM1, DOGMA-seq protocol (Version January 2021): staining,
  permeabilization 2a/2b, 10x modifications 1–4, notes 1–4, buffers, Oligos.
- 🟢 Cited 10x CG000338 Rev A Multiome guide (`.pdf.txt`, copied in this pass from the identical
  URL/md5 file in the REDeeM source directory): Introduction pp17–18 (lines 1205–1310), sequencing
  p68 (lines 4925–5030), Appendix "Sequences" p79 (lines 5625–5690).
- 🟢 MOESM3 Supplementary Table 1: TotalSeq-A barcodes (15 nt).
- Not inspected: PMC author-manuscript supplements (PMC download gate); BioLegend TotalSeq-A
  conjugate documentation; 10x Dual Index TT / Single Index N primer sequences.

## Boundary

🟢 DOGMA-seq is not an alias for ASAP-seq: it uses the Multiome kit, adds a gene-expression branch,
captures TotalSeq-A tags through the RNA bead's poly(dT) (no bridge oligo) and yields three libraries.

## Workflow (source order)

1. 🟢 Stain with TotalSeq-A antibodies after Fc block; wash; optional FACS of live cells (MOESM1 §1).
2. 🟢 Permeabilize: LLL (0.1% formaldehyde 5 min RT, glycine quench, LLL lysis 3 min on ice) or DIG
   (0.01% digitonin, 5 min on ice) (Methods; MOESM1 §2a/2b).
3. 🟢 Multiome workflow (CG000338 Rev A): transposition; GEMs; incubation gives barcoded transposed
   DNA and barcoded cDNA; silane cleanup; pre-amplification of both (guide pp17–18).
4. 🟢 Pre-amplification (step 4.1) spiked with 1 µl 0.2 µM ADT additive; after SPRI (step 4.3k) elute in
   100 µl EB; 25% to ATAC indexing, 35% to cDNA amplification, 35% to tag indexing.
5. 🟢 Tag PCR: KAPA 2×, SI-PCR + RPI-x, 95 °C 3 min, ~10 cycles (Methods: 10; MOESM1: ~10–12),
   expected product ~190 bp. MOESM1 note 4 offers an optional post-cDNA-amplification supernatant
   route, not drawn.

## Source-printed oligos

| Role | 5′→3′ as printed | Source |
|---|---|---|
| ADT additive | `CCTTGGCACCCGAGAATT*C*C` (* phosphorothioate) | Methods; MOESM1 |
| HTO additive | `GTGACTGGAGTTCAGACGTGTGC*T*C` | MOESM1 |
| SI-PCR | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC` | Methods; MOESM1 |
| RPI-x | `CAAGCAGAAGACGGCATACGAGATxxxxxxxxGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA` | Methods |
| RNA gel-bead primer PN-2000261 | `CTACACGACGCTCTTCCGATCT-N16-N12-T30VN` | CG000338 Rev A p79 |
| ATAC gel-bead primer | `AATGATACGGCGACCACCGAGATCTACAC-N16-CGCGTCTG` | CG000338 Rev A p79 |

🟡 The ADT additive is the 5′ 20 nt of the TSA small-RNA handle named in Methods (computed in the
model). MOESM1's example RPxx/D7xx members (`…GAGATAAGATTACGTGACT…`, `…GAGATCTTCCGGAGTGACT…`)
are examples only.

## Final libraries

- 🟢 RNA: `P5 (25 nt as printed)–i5 N10–TruSeq Read 1–N16–N12–T30VN–cDNA–Read 2 arm–i7 N10–P7 rc`
  (guide p79). Another 10x guide was reported in this repository's TEA-seq audit to print a 29-nt
  P5 prefix; not re-checked in this pass.
- 🟢 ATAC: `P5–N16–CGCGTCTG spacer–Read 1N–insert–Read 2N rc–i7 N8–P7 rc` (guide p79). (Model
  correction: the earlier page put an 8-nt i5 sample index at P5, omitted the spacer and placed the
  cell barcode after the mosaic end.)
- 🟡 TotalSeq-A tag: `P5–TruSeq Read 1–N16–N12–T30VN–[antibody barcode rc]–RPI-x arm rc–i7 N8–P7 rc`,
  derived from the bead primer, SI-PCR and RPI-x; barcode adjacency to poly(dT)VN is inferred.
  (Model correction: earlier 10-nt tag UMI and an invented "antibody-tag capture primer" step;
  the tag is captured by the RNA bead primer in the GEM and enriched by the ADT additive.)

## Read layout

🔴 Acquired read cycles are not stated in the inspected article sections; tag-library Read 2/Index 1
primers are not printed. No cycle counts are assigned. (Model correction: the earlier 28/90, 50/50 and
28/50 lengths had no source.) The guide recommends GEX 28/10/10/90 and ATAC 50/8/16/49 (p68); these
are vendor recommendations, not the authors' run.

## Open questions

- PMC author-manuscript protocol supplement vs MOESM1; BioLegend TotalSeq-A conjugate structure.
- 10x Dual Index TT and ATAC Single Index N primer sequences for exact endpoint geometry.
- Missing-information claims remain blocked (no two-reader gate in this pass).
