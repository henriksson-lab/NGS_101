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
