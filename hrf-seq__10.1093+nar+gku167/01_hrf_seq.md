# HRF-Seq

**Defining source:** Kielpinski and Vinther, *Nucleic Acids Research* (2014),
doi:10.1093/nar/gku167. The open article, methods and Figure 1 workflow were read.

## Evidence used

- 🟢 Hydroxyl radicals cleave solvent-exposed RNA backbones; treated and control samples
  distinguish accessibility-dependent termination.
- 🟢 Random primers with a 5′ Illumina tail create cDNAs terminating at radical-induced
  breaks.
- 🟢 A cDNA 3′ adapter carries a seven-nucleotide random barcode used for unique-molecule
  counting and bias correction.
- 🟢 Paired-end sequencing recovers both the probing and random-priming positions.
- 🟡 Exact oligo sequences remain to be transcribed from the source tables.

## Construct logic

HRF-Seq measures backbone solvent accessibility and tertiary packing. It is therefore
more directly an RNA-accessibility assay than a base-pairing assay alone.
