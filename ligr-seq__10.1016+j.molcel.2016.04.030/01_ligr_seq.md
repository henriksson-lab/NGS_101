# LIGR-seq — research notes

Defining source: Sharma *et al.*, [10.1016/j.molcel.2016.04.030](https://doi.org/10.1016/j.molcel.2016.04.030).

## Evidence used

- 🟢 LIGR-seq means LIGation of interacting RNA followed by high-throughput sequencing.
- 🟢 AMT photocrosslinking retains base-paired RNA duplexes; RNase R is used to deplete uncrosslinked linear RNA.
- 🟢 Proximity ligation joins the two RNA partners before crosslink reversal and cDNA library construction.
- 🟢 The complete Supplemental Experimental Procedures specify CircLigase ligation at 60 °C, RNase R enrichment, photoreversal, and library construction with the Clontech SMARTer Stranded RNA-Seq kit.
- 🟢 The paper describes a first-strand primer that adds a universal sequence, terminal-transferase activity followed by a universal second-strand primer, and indexed PCR. Its analysis explicitly trims a 5′ random barcode from single-end strand-specific reads.
- 🟢 Contemporary Clontech documentation identifies the PCR endpoint as Illumina/TruSeq-compatible and publishes the index barcode values. Standard Illumina sequencing-primer sites can therefore be derived centrally rather than treated as unknown.
- 🔴 The paper does not print the SMARTer N6 primer, SMART oligo, or complete PCR-primer sequences. The kit manual explicitly labels the first two as “Clontech proprietary sequences”, so their internal library geometry cannot be reconstructed from the primary or vendor sources.

## Protocol boundary

The page retains LIGR-seq's RNase R selection rather than presenting a generic RNA-duplex family schematic.

This is a genuine commercial-reagent disclosure boundary, not a sequence overlooked in the paper supplement. The schematic therefore shows the standard Illumina shell and the paper-reported inline random barcode, while leaving the proprietary priming junctions unresolved.
