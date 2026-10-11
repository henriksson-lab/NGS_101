# HDST — research notes

Vickovic et al., *Nature Methods* (2019),
[doi:10.1038/s41592-019-0548-y](https://doi.org/10.1038/s41592-019-0548-y).

## Sources read

- 🟢 Defining article and Online Methods (PMCID: PMC6765407).
- 🟢 Jemt et al. automated Spatial Transcriptomics library preparation,
  [doi:10.1038/srep37137](https://doi.org/10.1038/srep37137).

## Bead oligo and decoding

🟢 The published precursor is
`/AmC6/UUUUUGACTCGTAATACGACTCACTATAGGGACACGACGCTCTTCCGATCT-Spatial_barcode_Pool1-Bridge1`.
It contains an amine coupling modification, five uracils for later USER release, a T7
promoter, an Illumina handle, the first 15-base spatial block and its bridge.

🟢 Two more split-and-pool ligations add 14/15-base barcode pools and bridges, producing
65 × 211 × 211 combinations on 2-µm silica beads. A six-base UMI and poly(dT) capture
region follow the spatial address. The individual barcode and decoder sequences were an
Illumina product and are explicitly described as proprietary.

🟢 Beads are randomly loaded into more than 1.4 million hexagonal 2-µm wells. Fourteen
successive FAM/Cy3/dark decoder-hybridization images identify all three barcode blocks and
therefore map every bead to an x–y well before tissue is applied.

## Capture and library construction

1. 🟢 Place a 10- or 16-µm frozen section on the decoded array, fix, H&E stain and image.
2. 🟢 Permeabilize and synthesize cDNA overnight at 42 °C with SuperScript III. The
   bead-bound poly(dT) primer covalently links spatial address, UMI and transcript cDNA.
3. 🟢 Proteinase-K digest the section and release barcode–cDNA molecules with USER for
   2 h at 37 °C at the five-uracil tract.
4. 🟢 Process released material with the cited Spatial Transcriptomics CEL-seq-like
   workflow, which uses T7 linear amplification before final Illumina library completion.

## Reads

🟢 Mouse olfactory-bulb libraries used R1 125 and R2 175; breast-cancer libraries used
150 + 150. Read 1 positions 1–15, 31–45 and 61–76 were retained as the three-part spatial
barcode, and positions 77–82 are the six-base UMI. Read 2 is transcript sequence.

🟡 The schematic preserves those authoritative coordinates. It marks the two unavailable
bridge sequences and the compatible final Illumina shell inferred rather than silently
inventing the proprietary bases.
