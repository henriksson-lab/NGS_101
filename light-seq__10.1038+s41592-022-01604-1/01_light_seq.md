# Light-Seq — research notes

Kishi et al., *Nature Methods* (2022),
[doi:10.1038/s41592-022-01604-1](https://doi.org/10.1038/s41592-022-01604-1).

## Sources read

- 🟢 Defining open-access article and Methods (PMCID: PMC9636025).
- 🟢 Published Supplementary Information (MOESM1, especially its oligo table) and
  reporting summary (MOESM2, visually inspected because its PDF has no text layer).
- 🟢 Authors' complete Light-Seq protocol and oligo table,
  [doi:10.17504/protocols.io.x54v9jno4g3e/v1](https://doi.org/10.17504/protocols.io.x54v9jno4g3e/v1).
- 🟢 Authors' deposited analysis code, especially `ExtractTwoBarcodes.py` and
  `ExtractThreeBarcodes.py`, for the actual Read 1 parsing boundary.

## In-situ capture and light-directed indexing

1. 🟢 Fix and permeabilize cells or tissue, then reverse-transcribe total RNA with
   RT.5N.3G, `TTTACACGATTGAGTTATNNNNNGGG`. Five degenerate bases followed by three G
   bases prime without requiring poly(A); the 5′ constant becomes a docking site.
2. 🟢 Terminal transferase adds a poly(dA) tract to the free 3′ cDNA end.
3. 🟢 Hybridize one cnvK-bearing barcode strand throughout the sample, illuminate the
   chosen ROI and wash stringently. The photoactivated cnvK makes a covalent interstrand
   crosslink only in the illuminated area. Repeat with orthogonal barcodes for more ROIs.
4. 🟢 RNase H cleaves the RNA in RNA:cDNA hybrids and releases the selected, crosslinked
   barcode–cDNA complexes without destroying the specimen.

The three published barcode strands contain the same 27-base GATE PCR handle and 12
IUPAC-D UMI positions. Their ROI-identifying prefixes are `TATGGA`, `GTTAGG` and
`AGGGTA`; the authors' parser matches those six bases, followed by fixed `TGAGTTAT`.
The protocol preserves the ordered modifications as `[cnvK]`, `[Cy5-3]`, `[Cy3-3]` and
`[Fl-3]`.

## Crossing a covalent junction

🟢 The crosslink prevents ordinary PCR. GATC.20T primes the poly(dA) tract and Bst DNA
polymerase large fragment (NEB M0275L; Methods JATS line 204, protocols.io resource
line 82) extends through cDNA until it stops at the crosslink. (An earlier version of
this note named Bst 2.0 WarmStart; neither inspected source says so — corrected
2026-10-17.) A shared seven-base P
domain—potentially eight after Bst A-tailing—allows the newly synthesized end to invade
the opposite strand and resume synthesis on the barcode template. The uncrosslinked
product therefore carries GATC–transcript–P-domain–ROI–UMI–GATE on one strand. GATC and
GATE primers amplify that product; a small qPCR aliquot determines the endpoint cycle.

🟢 Reaction details from the published Methods (JATS lines 202–208): Maxima H Minus RT
with 1 µM RT primer (30 min 22 °C, 12-cycle 8–42 °C ramp, 42 °C 30 min); formamide
washes; A-tailing with terminal transferase, 1 mM dATP plus 25 µM ddATP (37 °C, 45 min);
RNase H displacement in ThermoPol (37 °C), eluate mixed with GATC.20T and heat
inactivated 75 °C 20 min; Bst LF 37 °C 30 min, 80 °C 20 min. A 5-µl aliquot is run in a
30-cycle SYBR qPCR with KAPA HiFi HotStart and GATE/GATC primers; the remainder is
bulk amplified (20 cycles for retina layers; 24 cycles for cell mixing, line 238).

## Nextera endpoint and reads

🟡 Published discrepancy — i5-side indexing primer. The Methods (JATS line 244) say
1.75 µl of "custom i5 primer (`GATE*.P5*` in Supplementary Table 6)" is added to the
indexing PCR. Supplementary Table 6 (MOESM1 text lines 1199–1331) and the protocols.io
"Sample indexing PCR" step 2c instead use S50X.GATE for that PCR and list `GATE*.P5*`
as the custom i5 *sequencing* primer. `GATE*.P5*` is the reverse complement of the
custom Read 1 primer and carries no P5, so it cannot add the flow-cell end. The model
uses S50X.GATE for PCR; whether the Methods sentence is a naming slip is not stated.

🟡 The authors' protocols.io v1 resource list (PDF text line 105) names Nextera DNA
kit FC-121-1030, but its tagmentation procedure (lines 2580–2582) explicitly uses
Nextera XT FC-131-1096, as does the published paper's Methods (JATS line 244).
The model uses XT because the procedure and defining article agree. Whether the
resource-list entry describes an alternate version remains unresolved.

🟢 After 1.2× AMPure XP cleanup, 2 ng amplicon is tagmented with Nextera XT (TD + ATM,
55 °C 5 min, NT stop); indexing PCR uses NPM for 12 cycles and 0.9× AMPure follows
(Methods line 244; protocols.io lines 2580–2690). Nextera XT tagmentation fragments the Light-Seq amplicon. Indexing PCR is asymmetric:
an ordinary Next.N7XX primer amplifies the tagmented end, while S5XX.GATE replaces an
ordinary Nextera S5 primer and selects the original GATE end. Thus Read 1 is not an
ordinary Nextera insert boundary. The required P5.GATE custom Read 1 primer places cycle
1 at the 12-base UMI, followed by the six-base ROI code, fixed `TGAGTTAT`, the eight-base
P domain and transcript cDNA. The authors extracted up to 40 transcript bases from R1;
R2 approaches the transcript from the random tagmentation boundary.

🟢 A GATE*.P5* custom i5 sequencing primer is published. It worked on HiSeq, but the
protocol reports inefficient i5 index sequencing on NovaSeq and therefore recommends a
unique i7 index for every pooled sample. The page displays the complete S502/N701 example;
the protocol also publishes S503 and S505–S508 plus N702–N706.

## Bounded source pass 2026-10-17 (round breadth-first-12)

Inspected: published JATS Methods lines 200–248 (RT, A-tailing, barcoding, displacement,
cross-junction synthesis, qPCR/bulk PCR, library preparation, sequencing); Supplementary
Table 6 (MOESM1 text lines 1100–1340); protocols.io v1 library-preparation, indexing-PCR
and sequencing sections (PDF text lines 2270–2770). Oligos re-checked against both
tables. Not followed: Primer Exchange Reaction predecessor (ref. 35), cnvK chemistry
references, authors' parsing code beyond earlier notes, lightseq.io updated protocols.
