# dsciATAC-seq — source note

The complete source audit and oligo transcription are retained in the companion
[dscATAC-seq / dsciATAC-seq research note](../dscatac-seq__10.1038+s41587-019-0147-6/01_dscatac-seq.md).

**Defining source:** Lareau CA *et al.* “Droplet-based combinatorial indexing for
massive-scale single-cell chromatin accessibility.” *Nature Biotechnology* 37, 916–924
(2019). doi:[10.1038/s41587-019-0147-6](https://doi.org/10.1038/s41587-019-0147-6).

dsciATAC-seq is the paper's combinatorial-indexing variant: a six-base barcode is placed
between s5 and the mosaic end during plate-based Tn5 assembly before droplet barcoding.

Graph (2026-10-17 bounded pass, Methods "Assembly of indexed Tn5 transposome complexes"
and "Species mixing controls", PMC10299900 lines 173–179): each indexed Read 1 adaptor
and the Read 2 adaptor are annealed to the 5′-phosphorylated, 3′-ddC ME complement
(1:1:2, 100 µM in 50 mM NaCl, 85 °C → 20 °C at −1 °C/min), diluted 1:1 in glycerol,
mixed 1:1 with 14.8 µM Tn5 for 30 min at room temperature; cells are tagmented with one
barcode per well (24-, 48- or 96-plex), pooled and washed in tagmentation buffer without
Tn5, then processed by the shared droplet workflow. Cells are loaded above Poisson, so a
droplet may hold several cells (distinguished by Tn5 barcode; BAP `--tn5-aware`) and
several beads. The in-droplet bead release, inferred gap fill, droplet PCR, emulsion
break, AMPure cleanups and second PCR follow the dscATAC-seq note.
