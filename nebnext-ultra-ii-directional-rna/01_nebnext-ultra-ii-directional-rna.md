# NEBNext Ultra II Directional RNA

## 1. What it is

A commercial strand-specific bulk RNA-seq workflow based on dUTP marking and selective removal of the second cDNA strand.

## 2. Sources and evidence

- 🟢 [NEB E7760/E7765 manual v7.0 (April 2025)](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7760_e7765.pdf), read from `/husky/henriksson/NGS101_papers/sources/nebnext-ultra-ii-directional-rna/manuale7760_e7765.pdf.txt`, and the [official RNA workflow](https://www.neb.com/en-us/products/next-generation-sequencing-library-preparation/library-preparation-for-illumina/rna-library-prep-for-illumina).
- 🟢 The kit supplies a Second Strand Synthesis Reaction Buffer **with dUTP Mix**, USER Enzyme and Strand Specificity Reagent (kit contents, manual p. 1). 🟡 The manual itself does not state the strand-selection mechanism: it adds USER to the ligation mixture (e.g. Step 4.6.5) without saying what is cleaved. Excision of the dU-marked second strand is inferred from those components; the official RNA workflow page that may state it was not re-read in the 2026-10-17 pass.
- 🟢 Strand Specificity Reagent (8 µl) is added to the **first-strand** reaction (Step 4.2.1); its composition is not printed. The earlier page claim that it "reinforces" the USER selection had no source passage and was removed.
- 🟢 The E7760/E7765 kit does **not** include the adaptor or amplification primers: its manual requires a separate NEBNext Oligo Kit or customer-supplied adaptor and primers (manual page 2). A specific oligo kit must therefore be named before exact indexed primers can be assigned to this workflow.
- 🟢 The manual prints Read 1 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and Read 2 `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` as **adapter trimming sequences** (page 3). These are sequences encountered in reads, not an oligo ordering sheet; they do not specify the ligation-ready adapter's complete strands or modifications.
- 🟡 Exact random-primer identities and kit-internal sequence boundaries are not disclosed; the page models their roles without assigning bases.

## 3. Molecular path

Fragmented RNA is random-primed, copied to first-strand cDNA, converted to a dUTP-bearing duplex, end-prepared and adapter-ligated. USER destroys the marked second strand before indexed PCR, so reads retain transcript orientation.

## 4. Reaction order and conditions (manual v7.0_4/25, Section 4)

🟢 Lines 3129–3780 of the manual text (Section 4, purified mRNA or rRNA-depleted RNA; Sections 1–3 and 5 add poly(A) isolation, rRNA depletion or FFPE handling before the same core):

1. Fragmentation and priming: RNA + First Strand Synthesis Reaction Buffer + Random Primers, 94 °C; 8–15 min for intact RNA (Table 4.1.3A), 5–8 min for degraded RNA.
2. First strand: + Strand Specificity Reagent + First Strand Synthesis Enzyme Mix; 25 °C 10 min, 42 °C 15 min, 70 °C 15 min.
3. Second strand: + Second Strand Synthesis Reaction Buffer with dUTP Mix + Enzyme Mix; 16 °C 1 h. No second-strand primer is named.
4. 1.8X SPRIselect / NEBNext Sample Purification Beads; elute in 0.1X TE.
5. Ultra II End Prep: 20 °C 30 min, 65 °C 30 min.
6. Ligation: diluted NEBNext Adaptor (from a separate oligo kit), Ligation Enhancer, Ultra II Ligation Master Mix; 20 °C 15 min.
7. USER Enzyme (3 µl) added to the ligation mixture; 37 °C 15 min. Optional size selection is done after USER digestion (Appendix 6.2).
8. 0.9X bead cleanup.
9. PCR with Ultra II Q5 Master Mix and the oligo kit's i7/Index and i5/Universal primers: 98 °C 30 s; 6–13 cycles of 98 °C 10 s, 65 °C 75 s (input-dependent, Table 4.8.3B); 65 °C 5 min.
10. 0.9X bead cleanup; library peak about 300 bp; peaks at about 80 bp (primers) or 128 bp (adaptor dimer) prompt a repeat cleanup.

The graph draws these steps and cleanups. Every section title specifies a **non-indexed** adaptor, but this manual does not print the adaptor or say what USER does to it; the page therefore does not draw an adaptor-opening step. The dual-index four-read final layout is one oligo-kit configuration (single-index kits use a Universal PCR Primer).

## 5. Open questions

- Which NEBNext oligo kit (and therefore adaptor structure and index layout) is assumed; the kit manual would need to be read before the USER action on the adaptor is drawn.
- Random primer length/composition, second-strand priming and Strand Specificity Reagent composition are not printed in the manual.
- The strand-selection mechanism is inferred, not stated in the manual (see section 2).

