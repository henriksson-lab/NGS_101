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

## Linked-method coverage

The article's library-preparation Methods explicitly inherit [Kielpinski et al.
2013, *Detection of reverse transcriptase termination sites using cDNA ligation
and massive parallel sequencing*](https://doi.org/10.1007/978-1-62703-514-9_13)
(reference 20), including cleanup and terminal-addition processing. This is a
distinct source from the SHAPE-seq predecessor already inspected. Its full
protocol was unavailable in the followup pass, so the source audit remains
incomplete. The original oligo-absence claim is resolved by Supplementary Table
S1; it is not evidence of missing published oligos.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC4005689 JATS38–40, Sequencing library preparation.

Named PrimeScript RT, CircLigase-mediated adapter ligation and Phusion4+20 PCR cycles. Confirmed200–600bp library selection and paired100 reads without changing the existing exact primer declarations.

Deferred: SupplementaryTableS1 complete modified oligos, cDNA-adapter geometry, RT recovery and final gel-selection graph stages; unusual published PCR dwell wording not silently corrected. This pass does not establish a complete source review or publication-wide absence.
