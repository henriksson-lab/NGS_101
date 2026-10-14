# HRF-Seq

**Defining source:** Kielpinski and Vinther, *Nucleic Acids Research* (2014),
doi:10.1093/nar/gku167. The open article, methods, Figure 1 workflow and
Supplementary Table 1 were read.

## Evidence used

- 🟢 Hydroxyl radicals cleave solvent-exposed RNA backbones; treated and control samples
  distinguish accessibility-dependent termination.
- 🟢 Random primers with a 5′ Illumina tail create cDNAs terminating at radical-induced
  breaks.
- 🟢 A cDNA 3′ adapter carries a seven-nucleotide random barcode used for unique-molecule
  counting and bias correction.
- 🟢 Paired-end sequencing recovers both the probing and random-priming positions.
- 🟢 Supplementary Table 1 prints the complete random and structure-cassette RT
  primers, phosphorylated/sealed seven-base-barcode ligation adapter, PCR forward primer,
  and four complete indexed reverse primers.

## Construct logic

HRF-Seq measures backbone solvent accessibility and tertiary packing. It is therefore
more directly an RNA-accessibility assay than a base-pairing assay alone.

The representative schematic uses the `GGAACT` reverse-primer member. Its Index 1 read
is `AGTTCC`, the reverse complement carried by the completed library. The other three
source-published members differ only at that six-base sample-index window.
