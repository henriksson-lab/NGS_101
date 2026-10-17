# 10x Chromium GEM-X Single Cell 3' Gene Expression v4 research notes

The source audit and chemistry comparison are preserved in the combined
[v2/v3/v4 note](../10x-chromium-single-cell-3-ge-v2/01_10x-3prime-v2-v4.md#v4-gem-x-cg000731-appendix).

The primary chemistry source is 10x Genomics user guide **CG000731 Rev A**.


## Current source and model boundary (16 October 2026)

The declaration follows the manual singleplex GEX workflow, rather than the separate
Kit C automation branch. Both the full source-printed bead oligo and final library
are now checked against the actual defining guide when its original is available.
The source-printed TT index pairs do not by themselves identify complete physical
PCR oligos. Sample-index segments remain placeholders; no plate well is assigned.

The current guide lists cDNA Primers **PN-2000089**. The predecessor guide prints
forward `CTACACGACGCTCTTCCGATCT` and reverse `AAGCAGTGGTATCAACGCAGAG` under that
same part number. These catalog-linked primers are declared locally and consumed
by cDNA PCR; assigning them to the current kit depends on formulation continuity.
The complete current-kit intermediate ends are not established, so exact PCR
geometry remains explicitly unresolved rather than validated from guessed ends.

The graph no longer inherits a shared assertion that **Ligation Mix PN-2001109**
has a known 3′-dT adapter end. Its physical oligos and ligation-ready geometry are
not established in this declaration. A final sequenced junction does not establish
the reagent strands or ends used before ligation. Previous speculative adapter
panels have been removed; the workflow renders this geometry limitation directly.
The remaining precursor drawings are explicitly conditional predecessor-based
illustrations, including any displayed RT-primer length and terminal-transferase tail.

The separate **CG000791 Rev A** automation supplement (October 2024), pp. 5–9 and
11, uses manually generated CG000731/CG000733 cDNA in automated Kit C
**PN-1000774**, with Ligation Mix **PN-2001317**. This different reagent number
cannot identify manual Ligation Mix PN-2001109. Its CG000802 link concerns
aliquot volumes for that automated branch; it is not substituted into the manual
GEX declaration. Independent linked-source review remains in progress, and no
universal nondisclosure claim is accepted here.

## Bounded source pass (17 October 2026, round breadth-first-06)

Read in **CG000731 Rev A** (md5 077492df…, the defining revision; `CG000731_RevA.pdf.txt`):
kit contents lines 224–320 (TSO B PN-2001027, Ligation Mix PN-2001109, cDNA Primers
PN-2000089, Amp Mix PN-2000047, Library Amp Mix PN-2000531); Introduction pp. 12–17
including the rendered figures on pp. 14 and 16; §1.5 GEM-RT (48 °C 45 min, 85 °C 5 min);
§2.1–2.3 (Recovery Agent, Dynabeads MyOne SILANE, cDNA PCR 98/63/72 °C with 13/12/11
total cycles for <500/501–6,000/>6,000 targeted cells, SPRIselect 0.6X); §3.1–3.6
(Fragmentation Mix of Fragmentation Buffer and Enzyme on 10 µl = 25% of cDNA, 32 °C 5 min
then 65 °C 30 min end repair/A-tailing in the same tube; double-sided SPRIselect
0.6X/0.8X; Ligation Mix + DNA Ligase 20 °C 15 min; 0.8X cleanup; sample-index PCR
98/54/72 °C, 5–16 total cycles by cDNA input, with Library Amp Mix or Amp Mix and
Dual Index TT Set A; double-sided 0.6X/0.8X); Step 4 run parameters (28/10/10/90
recommended cycles); appendix p. 80 (bead primer and both strands of the final library).
Rev B's appendix (`CG000731_RevB.pdf.txt` lines 5060–5067) prints the same sequences.

- 🟢 The bead primer and complete final library already matched the appendix; the
  intro text names a 30-nt poly(dT) and notes that a separate Capture Sequence 1 bead
  primer is unused for GEX.
- 🟢 The p. 14 figure draws a TSO ending rGrGrG pairing with a first-strand CCC but
  prints no TSO B bases; the TSO remains an inferred predecessor-style illustration.
- Model corrections: the graph now shows the Recovery Agent/oil removal, the
  Dynabeads cleanup, the 0.6X cDNA cleanup, both double-sided size selections and the
  0.8X post-ligation cleanup with discarded fractions; fragmentation and end
  repair/A-tailing are one Fragmentation Mix; source temperatures and cycle ranges are
  in the step labels. The sequenced product is the size-selected library.
- Still open: TSO B and Ligation Mix structures, physical TT index oligos and the
  current-kit identity of PN-2000089 (catalog-linked only). Not inspected this pass:
  CG000396, CG000479, CG00053, CG000748, CG000791, Rev B body text.
