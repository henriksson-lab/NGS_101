# Illumina 16S metagenomic sequencing

## 1. What it is

Illumina's demonstrated 16S workflow is a reusable two-stage amplicon design: PCR1 amplifies the bacterial V3–V4 locus with locus-specific primers carrying Nextera overhangs, and PCR2 uses those overhangs to add P5/P7, i5/i7 indexes and the standard sequencing-primer sites.

## 2. Authoritative source

🟢 Illumina, *16S Metagenomic Sequencing Library Preparation*, document 15044223 B. The guide gives both full PCR1 oligos, identifies the Klindworth V3–V4 primers, specifies a limited-cycle index PCR using Nextera XT dual indexes, and diagrams the final paired-end library.

## 3. Oligos

🟢 Forward PCR1 primer, 5′→3′:

`TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG`

🟢 Reverse PCR1 primer, 5′→3′:

`GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC`

The 5′ portions are universal Nextera overhangs; the 3′ portions bind the 16S locus. Degenerate IUPAC bases are copied verbatim from the guide.

## 4. Construct logic

PCR1 copies the V3–V4 interval and places a different universal overhang on each end. After bead cleanup, PCR2 primers anneal to those overhangs and add the two flow-cell arms plus combinatorial i5/i7 sample indexes. The final library uses the Nextera sequencing-primer set.

## 5. Scope

🟢 Illumina explicitly states that the same tailed-primer design can target other loci. This page nevertheless remains the named V3–V4 protocol; the fungal ITS page documents a materially different EMP one-stage construct and custom primers.
