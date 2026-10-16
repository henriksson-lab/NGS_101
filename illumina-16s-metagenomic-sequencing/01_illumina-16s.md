# Illumina 16S metagenomic sequencing

## 1. What it is

Illumina's demonstrated 16S workflow is a reusable two-stage amplicon design: PCR1 amplifies the bacterial V3–V4 locus with locus-specific primers carrying Nextera overhangs, and PCR2 uses those overhangs to add P5/P7, i5/i7 indexes and the standard sequencing-primer sites.

## 2. Authoritative source

🟢 Illumina, *16S Metagenomic Sequencing Library Preparation*, document 15044223 B. The guide gives both full PCR1 oligos, identifies the Klindworth V3–V4 primers, specifies a limited-cycle index PCR using Nextera XT dual indexes, and diagrams the final paired-end library.

🟢 The same guide's “Dual Indexing Principle” table (p. 23) gives all 12
N701–N712 i7 sequences and all 8 S501–S508 i5 sequences for the named
FC-131-1002 kit; its FC-131-1001 subset uses N701–N706 and S501–S504.
Thus the index families are published in the defining source. The guide
does not choose one pair for a particular sample, so the library schematic
retains variable index segments.

🟢 The linked Klindworth et al. 2013 original (doi:10.1093/nar/gks808,
PMC3592464 JATS line 84) independently prints the selected bacterial
S-D-Bact-0341-b-S-17 `CCTACGGGNGGCWGCAG` and
S-D-Bact-0785-a-A-21 `GACTACHVGGGTATCTAATCC` primer pair. These are exactly
the 3′ locus-specific portions of the two Illumina PCR1 oligos below. The source
does not assign a Nextera XT index pair to any Illumina-guide sample.

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
