# Omni-ATAC — research notes

Corces et al., *Nature Methods* (2017),
[doi:10.1038/nmeth.4396](https://doi.org/10.1038/nmeth.4396).

## Sources read

- 🟢 The defining article and Online Methods, available as the PMC author manuscript
  (PMCID: PMC5623106).
- 🟢 Supplementary Protocol 1, *Omni-ATAC protocol*, linked by the article and mirrored
  by the Corces laboratory.
- 🟢 Buenrostro et al., *Current Protocols in Molecular Biology* (2015),
  [doi:10.1002/0471142727.mb2129s109](https://doi.org/10.1002/0471142727.mb2129s109),
  for the inherited ATAC amplification workflow. Its methods refer the complete primer
  list to Buenrostro et al., *Nature Methods* (2013),
  [doi:10.1038/nmeth.2688](https://doi.org/10.1038/nmeth.2688); that paper's
  Supplementary Methods (MOESM200, opening lines) print Ad1_noMX and Ad2 members.

## What Omni changes

1. 🟢 Fifty thousand cells are resuspended in 50 µl ATAC-RSB containing 0.1% NP-40,
   0.1% Tween-20 and 0.01% digitonin and held on ice for 3 min. Multiple detergents make
   permeabilization less cell-type-dependent.
2. 🟢 One millilitre ATAC-RSB with 0.1% Tween-20 washes the permeabilized nuclei before
   they are pelleted. This removes released mitochondria rather than carrying them into
   transposition.
3. 🟢 The 50-µl transposition mix is 25 µl 2× TD buffer, 2.5 µl TDE1, 16.5 µl PBS,
   0.5 µl 1% digitonin, 0.5 µl 10% Tween-20 and 5 µl water. It is incubated for 30 min
   at 37 °C with mixing at 1,000 r.p.m. PBS is an intentional Omni-ATAC component.
4. 🟢 Transposed DNA is column-purified. The original protocol performs five initial
   amplification cycles, measures an aliquot by qPCR, and uses the fluorescence curve to
   select the minimum number of additional cycles.

## Shared ATAC library geometry

🟢 Omni-ATAC changes nuclei preparation and the transposition buffer, not the adapter
architecture. TDE1 transfers the canonical 19-base Nextera mosaic end. Of the possible
S5/S5, S7/S7 and S5/S7 end combinations, only the heterologous S5/S7 fragment can acquire
both flow-cell arms in PCR.

🟢 The defining protocol delegates amplification to the original ATAC procedure. The
published example pair is:

- `Ad1_noMX`: `AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTG`
- `Ad2.1_TAAGGCGA`: `CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGT`

The eight bases physically present in Ad2.1 are `TCGCCTTA`; the sample-sheet Index 1
value named by the primer is its reverse complement, `TAAGGCGA`. There is no i5 index in
this historical construct. Later dual-index ATAC implementations should not be projected
back onto the defining Omni-ATAC page.

## Sequencing boundary

🟢 The final library is paired-end Nextera geometry: Read 1 starts after the S5-side
mosaic end, Index 1 reports the eight-base Ad2 index, and Read 2 starts after the S7-side
mosaic end. The defining method uses multiple Illumina runs and does not make one genomic
read length part of the chemistry, so the page does not invent a fixed R1/R2 cycle count.
