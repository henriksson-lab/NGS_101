# GUIDE-seq — research notes

## Assay scope

The library captures genomic DNA beside dsODN insertions at nuclease-induced
breaks. The local declaration covers the tag, half-functional adapter, nested
library primers and sequencing readout; genomic targets remain schematic.

## Sources read

- 🟢 Tsai et al. 2015, DOI [10.1038/nbt.3117](https://doi.org/10.1038/nbt.3117),
  PMC4320685 — main text and Online Methods (JATS XML and PDF).
  - Online Methods, "dsODN for GUIDE-Seq": both modified tag strands.
  - Online Methods, "Library preparation": "end-repaired, A-tailed, and ligated to
    half-functional adapters, incorporating a 8-nt random molecular index".
- 🟢 Springer supplementary file MOESM21, **corrected 25 June 2015** (the correction
  notice un-swaps the `Nuclease_off_+_GSP1` / `Nuclease_off_-_GSP1` labels and restores
  the discovery thermocycling conditions):
  - Supplementary Fig. 1 — the library cartoon, which places the molecular index in the
    "Barcoded Y-adapter with Molecular Index (random 8-mer)" at the ligation step, and the
    oligo tag next to the P7 barcode in the final library.
  - Supplementary Table 4 — every common oligo: P7 adapters P701–P708, P5 adapter primers
    P5_1/P5_2, the two custom sequencing primers, Y-adapters A01–A16 with the MiSeq common
    adapter, and the four nested GSP primers.
  - Supplementary Methods, "Protocol for GUIDE-seq Library Preparation" — Y-adapter
    annealing, shearing, end repair/A-tailing, adapter ligation, PCR 1, PCR 2, and the
    MiSeq run.
- 🟢 Springer supplementary file MOESM22 — off-target site tables; contains no oligos.

2026-10-17 bounded re-reading (round 2026-10-17-breadth-first-11): Online Methods
"Isolation and preparation of genomic DNA for GUIDE-Seq" (JATS lines 129–131) and the
corrected MOESM21 Supplementary Methods protocol steps 1–12 (extracted lines 33074–33262),
including the equipment list (line 33026). The four AMPure XP SPRI cleanups, the DNAdvance
genomic-DNA isolation, the on-bead PCR 2 input and the qPCR normalization to 1.2 × 10¹⁰
molecules are now represented in `protocol_graph()` with their removed fractions.

Nothing in this note or in `tools/guide_seq.py` is taken from a secondary schematic.

## Oligos, verbatim with their modifications

🟢 Supplementary Table 4 and the Online Methods, exactly as printed
(`*` = phosphorothioate linkage, `[Phos]`/`P` = 5′ phosphate):

| Name | Sequence 5′→3′ |
|---|---|
| dsODN strand A | `5'-P-G*T*TTAATTGAGTTGTCATATGTTAATAACGGT*A*T-3'` |
| dsODN strand B | `5'-P-A*T*ACCGTTATTAACATATGACAACTCAATTAA*A*C-3'` |
| Y-adapter A01 (long strand) | `AATGATACGGCGACCACCGAGATCTACACTAGATCGCNNWNNWNNACACTCTTTCCCTACACGACGCTCTTCCGATC*T` |
| MiSeq common adapter | `[Phos]GATCGGAAGAGC*C*A` |
| P5_1 | `AATGATACGGCGACCACCGAGATCTA` |
| P5_2 | `AATGATACGGCGACCACCGAGATCTACAC` |
| Nuclease_off_+_GSP1 | `GGATCTCGACGCTCTCCCTATACCGTTATTAACATATGACA` |
| Nuclease_off_+_GSP2 | `CCTCTCTATGGGCAGTCGGTGATACATATGACAACTCAATTAAAC` |
| Nuclease_off_-_GSP1 | `GGATCTCGACGCTCTCCCTGTTTAATTGAGTTGTCATATGTTAATAAC` |
| Nuclease_off_-_GSP2 | `CCTCTCTATGGGCAGTCGGTGATTTGAGTTGTCATATGTTAATAACGGTA` |
| P701 | `CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTGACTGGAGTCCTCTCTATGGGCAGTCGGTGA` |
| custom Index 1 primer | `ATCACCGACTGCCCATAGAGAGGACTCCAGTCAC` |
| custom Read 2 primer | `GTGACTGGAGTCCTCTCTATGGGCAGTCGGTGAT` |

A02–A16 differ from A01 only in the 8-nt sample barcode; P702–P708 differ from P701 only
in the 8-nt P7 barcode. Every oligo above is reproduced in
`tools/guide_seq.py`, and `tools/selftest_source_oligos.py` re-checks the model's copies
against the extracted source text.

## How they interlock (computed in `tools/guide_seq.py`, not eyeballed)

🟢 Every identity below is asserted by `_require(...)` at module import, so the page
cannot be built if one stops holding.

- **A01 is a stock TruSeq P5 arm with two identifiers spliced in:**
  `A01 == illumina.P5 + TAGATCGC + NNWNNWNN + illumina.TRUSEQ_READ1`, i.e.
  {{= len(guide_seq.Y_ADAPTER_A01) }} nt.
- **The adapter is a 12-bp duplex with a 3′ dT overhang.** The common adapter's 5′
  `GATCGGAAGAGC` is `illumina.STEM`, the reverse complement of A01's
  `GCTCTTCCGATC`; A01's terminal `T` is left unpaired as the overhang that pairs with the
  dA-tailed fragment end.
- **The common adapter's 3′ `CA` has no partner.** It sits opposite A01's
  `{{= guide_seq.FORK_BASES }}` (the last two bases of the Read 1 arm before the stem) and
  neither base pairs. There is no second arm at all: no P7, no i7, no Read 2 site. *That*
  is the published "single-tail" / "half-functional" property — the molecule is a Y
  adapter with one arm missing rather than one arm blocked, so a fragment that never took
  up the dsODN carries P5 at both ends and cannot bridge a flow cell.
- **The P7 side is not TruSeq.** The P7-side arm is
  `{{= guide_seq.P7_ARM }}` ({{= len(guide_seq.P7_ARM) }} nt) — the 5′ tail of the GSP2
  primers, completed by P701–P708:
  `P701 == illumina.P7 + TCGCCTTA + P7_ARM[:-1]`. It is an AMP/Nextera-style arm, which is
  exactly why the run needs two spiked-in custom primers: the custom Read 2 primer **is**
  that arm, and the custom Index 1 primer is its exact reverse complement.
- **The GSP bodies are verbatim dsODN stretches.** `+_GSP1` carries strand B bases 1–22,
  `+_GSP2` strand B bases 13–34, `-_GSP1` strand A bases 1–29 and `-_GSP2` strand A bases
  7–33. GSP2 is therefore nested outward of GSP1 on the same strand, and the nested PCR
  replaces the tag's first 12 bases with the P7 arm.
- **P5_1 and P5_2 are nested prefixes of canonical P5**: `P5_2 == illumina.P5` and
  `P5_1 == illumina.P5[:26]`.

## Molecular order

🟢 Supplementary Methods, "Protocol for GUIDE-seq Library Preparation", in order:

1. **Y-adapter preparation.** A01–A16 are annealed to the MiSeq common adapter
   (95 °C → 60 °C → slow ramp to 4 °C). The molecular index enters the workflow here, as
   part of a synthetic oligo — not during PCR.
2. **dsODN capture.** The blunt, 5′-phosphorylated 34-bp dsODN with phosphorothioate
   linkages at both the 5′ and 3′ ends of both strands is nucleofected with the RGN
   components and captured at on- and off-target breaks by NHEJ.
3. **Genomic DNA isolation.** 🟢 Online Methods: solid-phase reversible immobilization
   magnetic beads (Agencourt DNAdvance). The graph now carries this bead isolation as an
   explicit extraction with its unbound fraction, ahead of shearing.
4. **Shearing and cleanup.** 400 ng genomic DNA in 120 µl, sheared to an average length of
   500 bp, then a 1X AMPure XP SPRI cleanup eluted in 15 µl.
   🟡 **Published source conflict:** the Online Methods name a "Covaris S200"; the
   Supplementary Methods protocol and equipment list name the Covaris S2 (S-2 Focused
   Ultra-sonicator). The model declares only the ~500 bp target, not an instrument model.
5. **End repair and A-tailing** — one tube: end-repair mix plus non-hot-start Taq, with
   12 °C / 37 °C / 72 °C steps. The 72 °C step is the dA addition.
6. **Adapter ligation.** Annealed Y adapter plus T4 DNA ligase, 16 °C then 22 °C. The
   adapter's 3′ dT pairs with the fragment's 3′ dA, so both ends of every fragment receive
   the adapter — but only one tail exists, so nothing is amplifiable without the tag.
   A 0.9X SPRI cleanup (eluted in 12 µl) follows, removing unligated adapter.
7. **PCR 1** — GSP1 + P5_1, touchdown programme. For discovery, separate master mixes and
   separate reactions for the `+` (sense) and `−` (antisense) GSP sets. Deep-sequencing
   validation instead uses one master mix with a normalized primer pool and the
   non-touchdown validation programme; that branch is not drawn.
   A 1.2X SPRI cleanup (eluted in 15 µl) follows; 🟢 the PCR 2 recipe takes "DNA sample
   with beads (from Step 8)", so the cleanup beads are carried into the nested reaction.
8. **PCR 2** — nested GSP2 + P5_2 + P7_#, same programme. GSP2 installs the P7 arm;
   P7_# completes it and adds the 8-nt P7 barcode. At least four different P7_# are used
   per run for image registration. A 0.7X SPRI cleanup (eluted in 30 µl) follows.
9. **Sequencing.** qPCR quantification, normalization to 1.2 × 10¹⁰ molecules, MiSeq
   Reagent Kit v2 300-cycle, **3 µl of 100 µM custom Index 1 primer into cartridge
   position 13 and 3 µl of 100 µM custom Read 2 primer into position 14**, cycles
   **151 | 8 | 16 | 151** "with the paired-end Nextera sequencing protocol".

## Final library and read layout

🟢 The library is a **hybrid**: a TruSeq P5/Read 1 side installed by ligation, and an
AMP-style arm on the P7 side installed by PCR.

```
{{oligo: guide_seq.FINAL_LIBRARY}}
```

- The **dA/dT junction exists on the P5 side only** — the ligated end. The P7 end is a
  PCR tail and has no dA at all. The drawn base at the junction is the adapter's 3′ dT;
  the fragment's dA is the facing base on the bottom strand.
- **Genomic DNA lies between the Read 1 arm and the dsODN**, because the adapter is at the
  sheared end while the tag is internal. Read 1 therefore reads genomic DNA from the shear
  site; Read 2 reads {{= len(guide_seq.DSODN_IN_LIBRARY) }} bp of dsODN and then crosses
  the integration junction into genomic DNA at the nucleotide-resolved breakpoint.
- 🟢 **The 16-cycle index-2 read is the molecular-index read.** Reading from P5 it reports
  8 cycles of Y-adapter sample barcode and then 8 cycles of `NNWNNWNN`; the 8-cycle
  index-1 read reports the P7 barcode. This is why the run is `151 | 8 | 16 | 151` rather
  than `151 | 8 | 8 | 151`. Cycle spans on the page are derived with
  `seqprimers.feature_spans()`, never written out by hand.
- **Index 2 mechanism and effective model anchor.** The recovered official
  *Indexed Sequencing Overview Guide* 15057455 **v08, November 2020**, p. 4
  (extracted lines 391–418), says MiSeq uses the surface-grafted P5 primer, then
  **seven chemistry-only cycles without imaging**, before imaging Index 2. The
  local `illumina.P5` reference marks the computational **post-extension read-start
  boundary**, not a claimed exact 29-nt physical grafted oligo. Its 16 imaged cycles
  cover sample barcode8 then molecular index8. The seven dark cycles are additional
  instrument chemistry, not extra reported index bases.
- The guide's actual revision history lists initial Rev. A in July 2014 and Rev. B
  in February 2015. This establishes document lineage, not the authors’ selected
  historical guide or instrument-recipe revision. The contemporary manufacturer
  custom-primer bulletin also states that MiSeq has no custom Index 2 primer
  option because it uses grafted P5. Neither source supplies a selected GUIDE-seq
  sample's adapter pair merely from the instrument name.

## `NNWNNWNN` is a constrained degenerate molecular index

🟢 Supplementary Methods: "The adapters also contain 8-mer NNWNNWNN (N = A, C, T, or G;
W = A or T) molecular indexes." Positions 3 and 6 are restricted to A/T, so the index has
4⁶ × 2² = 16 384 forms, not 4⁸ = 65 536. The molecular feature records
`encoding="degenerate"` to preserve those restrictions. The paper’s broader term
"random" does not imply that all eight positions are unrestricted N bases.

## Recorded source conflict: where the molecular index is installed

Two statements in the same paper disagree:

- 🟢 **Supplementary Methods and Supplementary Fig. 1 (the authority):** the index is part
  of the Y adapter, annealed before ligation and installed **at the ligation step**. Fig. 1
  labels it "Barcoded Y-adapter with Molecular Index (random 8-mer)" and Supplementary
  Table 4 prints it inside A01–A16.
- 🟡 **Main text, Stage II:** "Incorporation of a random 8 bp molecular barcode during the
  amplification process (Supplementary Fig. 1)". Read literally this would place it in the
  nested PCR, and the Online Methods likewise call it "a 8-nt random molecular index"
  while correctly attributing it to the half-functional adapters.

The protocol and the oligo table win: the index is a synthesized part of the ligated
adapter, with the explicitly restricted `NNWNNWNN` alphabet. An earlier version of this note and of
the model repeated the main text's looser phrasing ("nested amplification installs an 8-bp
random molecular barcode") and placed the index inside the insert, i.e. inside Read 1.
That was wrong on both counts.

## Boundary

GUIDE-seq measures breaks that capture the dsODN in living cells. It is distinct from
CIRCLE-seq, where purified genomic circles are cut by Cas9 in vitro.

## Recovered PMC supplement: an older source version

The actual PMC Cloud object `NIHMS647337-supplement-1.pdf` is 214,485,908 bytes.
It was publicly accessible; the former unavailable status arose from the generic
120 MB downloader cap. It was fetched through the standard source `Folder.save`
API with an explicit 215 MB cap for this identified document and validated as a
complete PDF. It is retained separately from corrected publisher MOESM21.

Direct text and visual comparison of **PMC p. 32** and **corrected publisher p. 33**
shows the GSP1 label swap described by the publisher's correction notice. The
older PMC table labels the GTTTA… GSP1 body `+` and ATACC… `−`; corrected MOESM21
labels ATACC… `+` and GTTTA… `−`. The two GSP2 rows are unchanged. The corrected
publisher file also distinguishes discovery and validation cycling, whereas the
older PMC protocol prints only the validation program. The model uses corrected
labels and does not substitute the older protocol. This version difference is
also disclosed on the public schematic.

Both copies print the same 151/8/16/151 sequencing layout and custom Index 1/Read 2
cartridge positions 13/14. The physical Index 2 mechanism comes from the separately
versioned manufacturer documentation above.

## Model questions and linked-source scope

- A01 and P701 are source-printed representatives. The completed library retains
  variable index windows rather than assigning them to an unspecified experimental
  sample. This is an unselected model value, not an assertion of publication absence.
- The current declaration does not assemble the strand-resolved PCR 1 template or
  sequential products of the three-primer PCR 2 reaction. All ordered library
  primers are source-positive; varying genomic flank lengths do not make their
  relative binding geometry inherently impossible. The structured PCR limitations
  describe model work still to do.
- The two discovery orientations are separate reactions. The graph pools their
  products after nested PCR, following the protocol's later library-normalization
  and pooling section; it does not prove an earlier experimental pooling time.
- The exact physical grafted-primer sequence and authors’ historical instrument
  recipe remain unassigned. The current computational read-start anchor must not
  be reinterpreted as an ordered grafted oligo or a 22-nt primer inferred by subtraction.
- The contemporary manufacturer bulletin links the actual MiSeq Custom Primers
  Guide15041638 and other platform-specific guides. Their source-specific scope is
  tracked in the review; a current guide is not automatically the 2015 revision.

The audit accepts no missing-information claim from these model boundaries.
