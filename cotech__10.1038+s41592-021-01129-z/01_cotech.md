# CoTECH

**Defining paper:** Xiong et al., *Nature Methods* (2021), doi:10.1038/s41592-021-01129-z.

## Sources read

- 🔴 Published article text: not available (publisher page only; manual row in the source
  manifest). Its Methods have **not** been read.
- 🟢 Published Supplementary Table 1 (Springer `MOESM3` XLSX): sheets `Custom_primers`
  (text lines 1–11), `PA-Tn5_barcodes` (12–63), `barcoded_oligo_dT` (64–162),
  `Nextera_index` (163–188), `Truseq_index` (189–).
- Historical only (excluded from current source authority): 🟡 bioRxiv preprint 10.1101/2020.10.15.339226v1, used **provisionally** only for step
  order because the published Methods are unreadable: Results lines 150–172 and the
  Figure 1 legend lines 347–358. It prints no STAR Methods text. It must be replaced by
  the published Methods when they are obtained.

## Historical provisional molecular interpretation (superseded by batch09 below)

🟢 (Table S1) Each PA-T5 and PA-T7 barcode primer is "annealed with TnExtMErev_primer, then
assembled with PA-Tn5"; the common annealing oligo is printed `5Phos/CTGTCTCTTATACACATCT`.
🟡 (preprint lines 159–161) each first-round well receives a unique T5/T7 combination, so
the chromatin first-round identity is the pair of T5 and T7 barcodes.

🟢 (Table S1) RT uses well-barcoded oligo dT (`RT-barcode_01`–`96`) with a `D` UMI window,
and a template-switch oligo `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G`. 🟢 The printed example directly contains an8base RTbarcode and D8 UMI; this molecular width does not establish acquired sequencing cycles.

🟡 (preprint lines 164–172, 349–358) Cells are pooled and redistributed into secondary
96-/384-well plates, both branches are preamplified, purified with AMPure XP and split in
halves. The RNA half is tagmented by Tn5 and its 3′ end is amplified "introducing a second
round of RNA index"; the DNA half is amplified with "a pair of common adaptor primers in
the first round and a pair of well-specific indexing primers in the second round of PCR".

**Model correction (this round):** the previous page attached the second-round index at
redistribution and labelled the library index a sample index. The second combinatorial
coordinate is the secondary well's library-PCR index; it is now a cell-barcode feature
added in `dna_pcr2` and `rna_pcr`. Barcode collisions were previously modelled as a
physically discarded waste product; they are a computational exclusion and were removed
from the physical graph.

## Oligos (verbatim from published Table S1)

🟢 `Mosaic-A_primer (MEA)` `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`; common annealing primer
`5Phos/CTGTCTCTTATACACATCT`; TSO `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G`; IS primer
`AAGCAGTGGTATCAACGCAGAGT`; P2 primer `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC`; connector A
`TCGTCGGCAGCGTCTCCACGC`; connector B `GTCTCGTGGGCTCGGCTGTCC`; TruSeq 1st-round forward
`ACACTCTTTCCCTACACGACGCTCTTCCGATCTTCGTCGGCAGCGTCTCCACGC` and reverse
`GACTGGAGTTCAGACGTGTGCTCTTCCGATCTGTCTCGTGGGCTCGGCTGTCCC`.

🟢 T5-1 `TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG` (24 T5 entries),
T7-1 `GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG` (25 T7 entries),
`RT-barcode_01` `TCAGACGTGTGCTCTTCCGATCTAACGTGATDDDDDDDDttttttttttttttttttttttttt`
(96 entries). Index families: Nextera i5 N501–N526, TruSeq i5 501–516, TruSeq i7 701–708.
The self-test checks the displayed oligos against the workbook text when present.

## Open questions

- 🔴 Published Methods: preamplification primers and cycling (connector A/B and IS/P2 are
  candidates, shown INFERRED), Tn5 loading for the cDNA branch (MEA + common annealing
  primer inferred), which i5/i7 carries the second-round RNA index, read lengths and
  which reads cover the T5/T7 barcodes, RT barcode and UMI.
- Predecessor CoBATCH (Wang et al., 2019) PA-Tn5 assembly not followed.
- Exact segment boundaries in the T5/T7 adaptors (barcode window length) not modelled
  base by base.


## Bounded native migration — 2026-10-18 batch09

Inspected 2 originals: Published SpringerMOESM3 SupplementaryTable1 Custom_primers1–11, PA-Tn5barcode/exampleT5/T7rows12–63, RTbarcode01/D8/polyTexample64onward and indexfamilyheadings; publishedMOESM1 SupplementaryNotes opening/entire shortNotes text comparing assay modalities/countingintronicRNA. Manifest originalidentities checked. OfficialNature definingarticle attempted once; nonretryableidpredirect, not read or counted. Mandatory bundle scraper/hits/modifications reviewed.

Removed multimodal_protocols Definition/ident, chain/role_library and active supersededpreprint order. SourceTable1 exact duplex transfer/common strands, RTprimer withD8 and sourceTSOrGrG+G are local structuredinputs. Conditional nativefirstwell/commonprocessing/precursor boundaries explicitly state publishedMethods/order unread; no preprint enzyme/primer assignment promoted. Mixed DNA/RNA class Scenes, actual incoming source-bounded ss/referencePCRsubstrates and residualoriginal lineage replace false DNA/RNA chains and expected-endpoint templates. Both endpoints preserve knownT5/T7orRT/UMI features; secondarywell width/run unknown. This is a sourcecomponent migration with unverified complete execution sequence, not full published chemistry validation.

Deferred: Published definingMethods including actual operationorder/reagent assignments, preamp/pooling/redistribution/recovery and RNA/DNA precursor formation; Cited CoBATCH/predecessor and native exact first-cycle strand geometry; Native secondarywell indexassignment/member/width and actual read/indexacquisition; sourceavailability unresolved, no accepted absence. This capped pass does not accept source-wide absence or establish complete experimental validation. Earlier dated/provisional interpretations above are historical; this entry describes the current source scope.
