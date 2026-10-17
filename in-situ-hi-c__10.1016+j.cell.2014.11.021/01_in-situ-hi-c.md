# in situ Hi-C

## 1. What it is

In situ Hi-C converts pairs of chromosomal loci held close together in a crosslinked
nucleus into chimeric DNA molecules whose two ends can be identified by paired-end
sequencing. The defining implementation moved restriction digestion and proximity
ligation into intact nuclei and used a four-base cutter for dense genome-wide sampling.

Defining paper: Rao *et al.*, “A 3D map of the human genome at kilobase resolution reveals
principles of chromatin looping,” *Cell* (2014),
doi:[10.1016/j.cell.2014.11.021](https://doi.org/10.1016/j.cell.2014.11.021).

## 2. Sources read

- 🟢 The PMC author manuscript, PMCID PMC5635824: Results, “In situ Hi-C methodology and
  maps,” and Experimental Procedures, “In situ Hi-C Protocol.”
- 🟢 Extended Experimental Procedures I.a.1 through the ENCODE transcription of the Rao
  supplement. ENCODE says it is a transcription and marks its own deviations; only the
  undeviated Rao procedure is used here.
- 🟢 Bounded pass 2026-10-17 (round 11): the published Cell article with its Extended
  Experimental Procedures (non-publisher course mirror carrying the Cell DOI and S-page
  pagination; `Rao-3D-map-human-genome-cell-2014-extended.pdf.txt` lines 6632–6802,
  section I.a.1 steps 1–70, and lines 6888–6890, I.a.10). This is the original behind the
  ENCODE transcription; the gated PMC supplement zip itself remains unread.
- 🔴 The Rao procedure specifies an “Illumina indexed adapter” and Illumina PCR primers
  but does not identify a kit generation, index count, or oligo sequences. The final
  adapter regions therefore remain role-labelled unknowns rather than canonical TruSeq.

## 3. Restriction digestion as a strand operation

🟢 MboI recognizes palindromic `GATC`: the top strand is cut as `5′-^GATC-3′`
and its paired strand as `3′-CTAG^-5′`.

The cuts occur at different coordinates on the two strands. Each product therefore has a
four-nucleotide 5′ `GATC` cohesive end. `lib/restriction.py` stores those two cut coordinates
and derives both product strands from them; the overhang is not copied into two independent
drawings.

🟢 The paper permits MboI or DpnII. The detailed procedure uses 100 U MboI at 37 °C for at
least two hours or overnight, then heat-inactivates it at 62 °C for 20 minutes.

## 4. Biotin fill-in and the contact junction

🟢 Klenow fills the recessed 3′ ends using dCTP, dGTP, dTTP and biotin-14-dATP. For one
orientation, the new top strand is `GATC`; for the other, the new bottom strand is `GATC`
when written 5′→3′. Each filled end therefore receives one biotin-dA.

🟡 Blunt ligation of two filled MboI ends produces the top-strand motif `GATCGATC`. This
is derived by `Digest.fill_in().junction`; the schematic does not carry a separately typed
junction sequence. The two biotins occur on opposite strands around the ligation boundary.

🟢 T4 DNA ligase acts for four hours at room temperature with slow rotation. Because the
nuclei remain intact, ends constrained near one another are favored over random DNA ends
in bulk solution.

## 5. Junction selection and library preparation

🟢 After proteinase K treatment and crosslink reversal, DNA is sheared to 300–500 bp and
size-selected. Streptavidin T1 beads capture molecules containing biotinylated ligation
junctions. On-bead end repair both repairs sheared ends and removes biotin from unligated
ends; Klenow exo-minus then adds 3′ dA.

🟢 An Illumina indexed adapter is ligated, and the library is amplified from the beads with
Illumina primers. The paper reports paired-end sequencing. The reads start at opposite
genomic ends of the sheared molecule, while the selected contact junction is internal.

## 5a. Step order from the published Extended Experimental Procedures

🟢 I.a.1 (lines 6632–6802) gives this order, now declared in the protocol graph:
lysis in 0.2% Igepal buffer and 0.5% SDS permeabilization at 62 °C quenched by Triton X-100
(steps 7–11); 100 U MboI (NEB R0147) overnight or ≥2 h, 62 °C inactivation (12–13);
Klenow fill-in with biotin-14-dATP and dCTP/dGTP/dTTP (14–15); T4 DNA ligase (NEB M0202)
4 h at room temperature (16–17); proteinase K/SDS 55 °C then NaCl 68 °C (18–19); ethanol
precipitation (21–25); Covaris LE220 shearing to 300–500 bp (26); **double-sided AMPure XP
size selection before pull-down** (0.55×, then +30 µl; 28–36; optionally moved to final
amplicons for <2×10⁶ cells); MyOne Streptavidin T1 capture with 55 °C TWB washes (37–41);
on-bead end repair with T4 PNK, T4 DNA polymerase and Klenow, which also removes biotin
from unligated ends (43–44); Klenow exo minus dA-tailing (48–49); NEB Quick ligase with
3 µl “Illumina indexed adapter”, recording the sample index (53–55); 4–12 PCR cycles off
the beads with “Illumina primers and protocol (Illumina, 2007)”, optionally after 98 °C
release (60); bead removal and AMPure XP 0.7× plus a 70 µl re-cleanup (61–70).
Sequencing used MiSeq, HiSeq 2000 and HiSeq 2500 paired-end runs (I.a.10, line 6890);
read and index cycle counts are not stated in these passages.

Not inspected in this capped pass: the cited Illumina (2007) genomic-DNA sample-prep
guide, the gated PMC supplement zip, and dilution/variant sections I.a.2–I.a.9 beyond
headings.

## 6. Scope

This page is the Rao 2014 MboI in situ Hi-C workflow. Dilution Hi-C, single-cell Hi-C,
DNase Hi-C, Micro-C, capture Hi-C, HiChIP and later Hi-C 2.0/3.0 variants are distinct
protocols and are not silently folded into it.
