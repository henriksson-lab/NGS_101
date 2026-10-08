# Earth Microbiome Project ITS amplicon sequencing

## 1. What it is

The Earth Microbiome Project fungal ITS workflow amplifies ITS1 with primers that already contain complete Illumina arms and a twelve-base Golay sample barcode. Unlike the two-stage Illumina 16S workflow, this is a fusion-primer library and uses locus-specific custom sequencing primers extending beyond the PCR-primer boundary.

## 2. Authoritative source

🟢 Earth Microbiome Project, *ITS Illumina Amplicon Protocol* (EMP.ITSkabir), and the linked November 2016 `ITS1f_ITS2` oligo workbook. The workbook publishes the fusion constructs and all three custom sequencing primers.

## 3. Published fusion primers

🟢 Forward ITS1f construct, 5′→3′:

`AATGATACGGCGACCACCGAGATCTACAC GG CTTGGTCATTTAGAGGAAGTAA`

This is the P5-side Illumina arm, a `GG` linker and ITS1f.

🟢 Reverse ITS2 construct, 5′→3′:

`CAAGCAGAAGACGGCATACGAGAT NNNNNNNNNNNN CG GCTGCGTTCTTCATCGATGC`

This is the P7-side arm, a twelve-base Golay barcode, a `CG` linker and ITS2. Spaces above only expose the published components.

## 4. Read geometry

🟢 Read 1: `TTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCC`.

🟢 Read 2: `CGTTCTTCATCGATGCVAGARCCAAGAGATC` (`V` and `R` are degenerate positions).

🟢 Index 1: `TCTCGCATCGATGAAGAACGCAGCCG`.

These produce paired-end reads plus one custom index read over the twelve-base Golay barcode; there is no i5 read. For an exact computed landing, the schematic uses the valid Read 2 pool member `V=G, R=A` and labels that choice rather than treating it as inferred chemistry.

## 5. Scope

This page documents the fusion-primer EMP implementation, not every possible ITS region or primer pair. Primer choice remains taxon- and study-dependent.
