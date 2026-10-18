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

## Native ownership migration, 2026-10-18

Sources read: Published PMC4005689 JATS Methods24–48 (substrates, Fentonprobing, librarypreparation,paired100,readpreprocessing); published supplementaryFile009.pdf TableS1 lines61–110 (random/cassetteRT,5′PHO/3′NHC3 barcodeadapter,forward andfour indexedPCR oligos).

Removed rna_end_readout steps/finisher and invented exact duplex PCR template. Native selectedribosome/randomRT graph has actualsingle-strandRTproduct, source5′phosphorylated7ntbarcodeadapter with nonnucleotide3′aminoblock, covalent3′cDNA–5′adapterligation and cleanup beforePCR. PCR unavailable on actualblockedSSprecursor; final reference retains sourceoligos andfeatureidentity and corrects random-primer coordinate to reverse-complement orientation. Source4+20Phusion cycles,200×dilution,200–600bpfinalselection andPE100 retained. Indexsequence6nt is not guessed as an acquired indexread. RNaseP/cassetteRT remains a separately identified sourcevariant.

Model boundaries: Nativeforward-first/reverse-secondPCR molecularproof; RNaseP/cassetteRT graph; exactindexacquisition; cleanup component partitions and linked predecessor manual not read thispass. Mandatory bundle oligo scraper was run and relevant original passages checked. This capped pass accepts no publication-absence claim.

## Independent review corrections, 2026-10-18

Independent reader inspected Published PMC4005689 Methods33–49. Actual blocked single-strand ligation precursor/PCR geometry is sound, but RT template RNA recovery remains unresolved through delegated RNAClean purification.

Added explicit molecular RNA-template output and unresolved downstream recovery after RT. PCR Scene is explicitly the tracked cDNA target component; it does not assert an RNA-free tube or infer template destruction in RNAClean. This correction used that bounded review and the existing source scope; no additional original was read by the implementation reviewer. No publication-absence claim is accepted.
