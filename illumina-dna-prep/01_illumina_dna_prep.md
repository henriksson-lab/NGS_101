# Illumina DNA Prep

## Sources read

🟢 Illumina, *Illumina DNA Prep Reference Guide*, document 1000000025416 v12, pp. 1–5 and 23–32. It defines bead-linked tagmentation, post-tagmentation cleanup, input-dependent limited-cycle PCR, bead cleanup and pooling.

🟢 Illumina, *Illumina Adapter Sequences*, document 1000000002694 v22, supplies the canonical Nextera adapter and sequencing-primer sequences through `lib/nextera.py`.

## Molecular workflow

1. Bead-linked transposome complexes bind genomic DNA and simultaneously fragment and transfer partial adapters during a 15-minute tagmentation.
2. The bead-bound product is washed after tagmentation.
3. Limited-cycle PCR completes both flow-cell arms and introduces unique dual indexes.
4. Illumina Purification Beads remove short products before pooling.

## Input-dependent behavior

🟢 The guide accepts 1–500 ng or more. It recommends 100–500 ng for human and other complex genomes, while small genomes can use 1 ng.

🟢 At 100–500 ng, the finite bead-linked transposome population provides normalized yield suitable for equal-volume pooling. Below 100 ng, input and libraries require quantification and the PCR cycle count changes with input.

🔴 Illumina does not disclose the reagent-level bead/transposome composition. The schematic therefore treats BLT as reaction geometry while using the public Nextera endpoint and computed standard primer sites.
