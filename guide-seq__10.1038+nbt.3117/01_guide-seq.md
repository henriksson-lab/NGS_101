# GUIDE-seq — research notes

## Assay target

The sequenced material is human genomic DNA flanking nuclease-induced double-strand breaks
in living cells. This protocol does not target, amplify or reconstruct a pathogen genome:
every oligo it uses is a host-genome or library oligo (a non-coding 34-bp tag duplex, a
half-functional Y adapter and nested Illumina primers). The note is here so that the
protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen
work.

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
3. **Shearing.** 400 ng genomic DNA, Covaris S2, average length 500 bp.
4. **End repair and A-tailing** — one tube: end-repair mix plus non-hot-start Taq, with
   12 °C / 37 °C / 72 °C steps. The 72 °C step is the dA addition.
5. **Adapter ligation.** Annealed Y adapter plus T4 DNA ligase, 16 °C then 22 °C. The
   adapter's 3′ dT pairs with the fragment's 3′ dA, so both ends of every fragment receive
   the adapter — but only one tail exists, so nothing is amplifiable without the tag.
6. **PCR 1** — GSP1 + P5_1, touchdown programme. For discovery, separate master mixes and
   separate reactions for the `+` (sense) and `−` (antisense) GSP sets.
7. **PCR 2** — nested GSP2 + P5_2 + P7_#, same programme. GSP2 installs the P7 arm;
   P7_# completes it and adds the 8-nt P7 barcode. At least four different P7_# are used
   per run for image registration.
8. **Sequencing.** qPCR quantification, normalization to 1.2 × 10¹⁰ molecules, MiSeq
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
- 🟡 **Index 2 primer.** The protocol spikes custom primers into cartridge positions 13
  and 14 only, so nothing is added for the i5 read. MiSeq is a forward-strand instrument,
  where the grafted flow-cell P5 oligo primes the index-2 read (Illumina "Indexed
  Sequencing Overview Guide" #15057455). The model declares `illumina.P5` for that role.
  The 16-cycle length and the identifiers it covers are 🟢; the identity of the priming
  oligo is standard-instrument behaviour rather than a GUIDE-seq statement.

## `NNWNNWNN` is degenerate, not random

🟢 Supplementary Methods: "The adapters also contain 8-mer NNWNNWNN (N = A, C, T, or G;
W = A or T) molecular indexes." Positions 3 and 6 are restricted to A/T, so the index has
4⁶ × 2² = 16 384 forms, not 4⁸ = 65 536. The molecular feature records
`encoding="degenerate"`; `"random"` would be wrong.

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
adapter, and it is degenerate rather than random. An earlier version of this note and of
the model repeated the main text's looser phrasing ("nested amplification installs an 8-bp
random molecular barcode") and placed the index inside the insert, i.e. inside Read 1.
That was wrong on both counts.

## Boundary

GUIDE-seq measures breaks that capture the dsODN in living cells. It is distinct from
CIRCLE-seq, where purified genomic circles are cut by Cas9 in vitro.

## Open questions

- 🟡 The two discovery orientations are run as separate master mixes with the `+` and `−`
  GSP pairs; the paper does not say whether the two products are pooled before or after
  qPCR normalization. The model pools after the nested PCR, which is the only point at
  which pooling is described ("normalize the total set of libraries").
- 🔴 Neither PCR is drawn with exact primer-binding coordinates. PCR 1 cannot be: the
  sheared genomic flank has no fixed length and the half-functional adapter is on *both*
  fragment ends, so P5_1 has no unique site on one drawn duplex. PCR 2 additionally uses
  three primers. Both are recorded as structured PCR limitations rather than invented
  geometry.
- 🔴 Supplementary Table 4 gives no index-length statement for its own barcodes beyond
  what the printed oligos show; the 8-nt widths used here are read off the printed
  sequences (A01 `TAGATCGC`, P701 `TCGCCTTA`), which is also what the 8-cycle index-1 read
  confirms.
