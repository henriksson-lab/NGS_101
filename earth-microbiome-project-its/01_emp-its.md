# Earth Microbiome Project ITS amplicon sequencing

## 1. What it is

The Earth Microbiome Project fungal ITS workflow amplifies ITS1 with primers that already contain complete Illumina arms and a ten-base Golay sample barcode. Unlike the two-stage Illumina 16S workflow, this is a fusion-primer library and uses locus-specific custom sequencing primers extending beyond the PCR-primer boundary.

## 2. Authoritative source

🟢 Earth Microbiome Project, *ITS Illumina Amplicon Protocol* (EMP.ITSkabir). The maintained protocol publishes both ordered primer constructs and explains that the forward and reverse sequencing primers contain respectively 19 and 15 additional 3′ bases beyond the amplification primers.

## 3. Published fusion primers

🟢 Forward ITS1f construct, 5′→3′:

`AATGATACGGCGACCACCGAGATCTACAC GG CTTGGTCATTTAGAGGAAGTAA`

This is the P5-side Illumina arm, a `GG` linker and ITS1f.

🟢 Reverse ITS2 construct, 5′→3′:

`CAAGCAGAAGACGGCATACGAGAT NNNNNNNNNN CG GCTGCGTTCTTCATCGATGC`

This is the P7-side arm, a ten-base Golay barcode, a `CG` linker and ITS2. Spaces above only expose the published components.

## 4. Read geometry

The custom sequencing primers extend into the locus-specific region so that their melting temperatures remain suitable. The public protocol specifies the extensions but does not print all four complete sequencing-primer oligos on the protocol page; the schematic therefore marks those landing regions and does not substitute standard TruSeq primers.

## 5. Scope

This page documents the fusion-primer EMP implementation, not every possible ITS region or primer pair. Primer choice remains taxon- and study-dependent.
