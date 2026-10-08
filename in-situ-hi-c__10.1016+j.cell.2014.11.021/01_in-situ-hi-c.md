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

## 6. Scope

This page is the Rao 2014 MboI in situ Hi-C workflow. Dilution Hi-C, single-cell Hi-C,
DNase Hi-C, Micro-C, capture Hi-C, HiChIP and later Hi-C 2.0/3.0 variants are distinct
protocols and are not silently folded into it.
