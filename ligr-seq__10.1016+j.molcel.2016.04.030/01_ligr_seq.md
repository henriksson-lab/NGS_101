# LIGR-seq — research notes

Defining source: Sharma *et al.*, [10.1016/j.molcel.2016.04.030](https://doi.org/10.1016/j.molcel.2016.04.030).

## Evidence used

- 🟢 LIGR-seq means LIGation of interacting RNA followed by high-throughput sequencing.
- 🟢 AMT photocrosslinking retains base-paired RNA duplexes; RNase R is used to deplete uncrosslinked linear RNA.
- 🟢 Proximity ligation joins the two RNA partners before crosslink reversal and cDNA library construction.
- 🟢 The complete Supplemental Experimental Procedures specify CircLigase ligation at 60 °C, RNase R enrichment, photoreversal, and library construction with the Clontech SMARTer Stranded RNA-Seq kit.
- 🟢 The paper describes a first-strand primer that adds a universal sequence, terminal-transferase activity followed by a universal second-strand primer, and indexed PCR. Its analysis explicitly trims a 5′ random barcode from single-end strand-specific reads.
- 🟢 The available 2013 Clontech *Product Components List* for catalog numbers 634836–634839 says the kit includes SMARTer Stranded RNA-Seq components and an Illumina Indexing Primer Set for Illumina-compatible libraries. Its second page points to a separate user manual.
- 🟢 An online transcript of the 2014 Clontech *SMARTer Stranded RNA-Seq Kit User Manual* (revision 050714, pp. 6–8) marks its SMARTer Stranded Oligonucleotide and N6 primer as proprietary sequences. It prints the twelve six-base Illumina Indexing Primer Set barcode values, including `ATCACG` for index 1. The original vendor PDF is unavailable in this source bundle. A 2019 Takara revision confirms the same kit family and lists a distinct HT index set; neither revision identifies which kit generation or index the LIGR-seq study used.
- 🔴 Neither the original LIGR-seq supplement nor the 2014 manual prints the full SMARTer and N6 oligo bases or identifies which index primer the LIGR-seq libraries used. The exact selected PCR primer and internal junctions remain unresolved.

## Protocol boundary

The page retains LIGR-seq's RNase R selection rather than presenting a generic RNA-duplex family schematic.

The schematic shows an inferred standard Illumina shell and the paper-reported inline random barcode. Exact kit priming junctions and the selected index remain unresolved after the 2014 manual check.
