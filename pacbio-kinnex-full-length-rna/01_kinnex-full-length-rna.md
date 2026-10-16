# PacBio Kinnex full-length RNA

## 1. What it is

Kinnex full-length RNA converts full-length cDNA molecules into ordered multi-insert arrays, then encloses each array in a SMRTbell for PacBio HiFi sequencing and computational segmentation.

## 2. Source read

🟢 Pacific Biosciences, *Preparing Kinnex libraries using the Kinnex full-length RNA
kit*, document 103-238-700 Rev09, August 2026.

🟢 PacBio, [SMRT Link v13.1 user guide](https://www.pacb.com/wp-content/uploads/SMRT-Link-User-Guide-v13.1.pdf), pp. 100–102: the full-length RNA default is MAS-Seq Adapter v3 (MAS8) for segmentation and Iso-Seq v2 Barcoded cDNA Primers for primer trimming. The guide prints the twelve 5′ barcoded cDNA primer sequences and the 3′ primer in its example. 🟢 The official [Kinnex data-release README](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/README.txt) names the source dataset and analysis workflow.

🟢 The README links PacBio's [Kinnex full-length RNA application note, 102-326-591 Rev03](https://www.pacb.com/wp-content/uploads/Application-note-Kinnex-full-length-RNA-kit-for-isoform-sequencing.pdf). Its p. 3 workflow text identifies Iso-Seq express 2.0 cDNA input, 5′ cDNA-level barcoding, four-plex Kinnex terminal-adapter indexing, and nuclease removal of partial arrays. This supports array topology but does not specify full primer or adapter arms.

## 3. Construct-changing operations

1. 🟢 Poly(A)-primed, template-switch cDNA synthesis is followed by Iso-Seq barcode PCR.
2. 🟢 Eight parallel Kinnex PCRs add orientation-specific segmentation sequences.
3. 🟢 Equal volumes of the eight products are pooled.
4. 🟢 Kinnex enzyme, ligase and a barcoded Kinnex terminal adapter assemble the cDNA
   segments into a linear array at 45 °C for 60 min; DNA repair follows for 30 min.
5. 🟢 Nuclease treatment removes incomplete molecules with exposed ends. The recovered
   product is the SMRTbell library used for sequencing-primer annealing and polymerase
   binding.

## 4. Evidence boundary

🟢 PacBio’s [technical overview](https://www.pacb.com/wp-content/uploads/Technical-overview-Kinnex-kits-for-single-cell-RNA-full-length-RNA-and-16S-rRNA-sequencing.pdf), p. 27, and its [Kinnex index FASTA](https://www.pacb.com/wp-content/uploads/Kinnex_barcoded_adapter_indexes.fasta_.zip) both give the terminal-adapter barcode cores: bcM0001 `ACAGTC`, bcM0002 `ATGACG`, bcM0003 `CACGTG`, bcM0004 `CATCGC`. The [multiplexing resource page](https://www.pacb.com/multiplexing/) links that archive to the Kinnex full-length RNA kit 103-072-000 and distinguishes these sample indexes from the A–H concatenation identities.

🟢 The official [MAS-Seq Adapter v3 (MAS8) FASTA](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-MAS_adapters/MAS-Seq_Adapter_v3/mas8_primers.fasta), records A–H and Q, gives the short **analysis reference motifs** below. Its adjacent `MAS-Seq_Adapter_v3.barcodeset.xml` names `mas8_primers.fasta` as the BarcodeSet resource. The SMRT Link guide p. 102 explains that a segmentation FASTA is an ordered list of adapter recognition sequences. The guide's custom-FASTA *example* prints 16-base prefixes for A–E; the versioned MAS8 dataset FASTA gives 17-base A–H records and a 14-base Q record. These are different reference contexts, not evidence that the physical ligation junction or full PCR oligo has one length or the other.

| Motif | 5′→3′ analysis reference |
|---|---|
| A | `AGCTTACTTGTGAAGAT` |
| B | `ACTTGTAAGCTGTCTAT` |
| C | `ACTCTGTCAGGTCCGAT` |
| D | `ACCTCCTCCTCCAGAAT` |
| E | `AACCGGACACACTTAGT` |
| F | `AGAGTCCAATTCGCAGT` |
| G | `AATCAAGGCTTAACGGT` |
| H | `ATGTTGAATCCTAGCGT` |
| Q | `AGTAGCTGTGTGCA` |

🟢 The SMRT Link guide p. 101 and the linked [Iso-Seq v2 primer FASTA](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta) print the twelve barcoded 5′ cDNA amplification primers and opposite 3′ primer. For example, bc01 is `CTACACGACGCTCTTCCGATCTACTACACGCAATGAAGTCGCAGGGTTGGG`; the 3′ primer is `AAGCAGTGGTATCAACGCAGAGTAC`. The selected bc01–bc12 depends on the sample.

🔴 Complete Kinnex orientation-PCR primer constructs, covalent segmentation-junction geometry, terminal-adapter arms and end geometry, and sequencing-primer binding bases remain unverified. The short analysis references above cannot be substituted for physical oligo specifications. Array order, nuclease selection and SMRTbell topology can still be represented without guessed bases or an overhang assignment. Two independent reads now cover the index archive, MAS8 reference, SMRT Link guide, multiplexing page and dataset README; the linked application note has one located read. The remaining absence claim stays blocked pending its second read and any further linked component documentation.

## 5. Read layout

The PacBio sequencing primer binds a terminal-adapter site whose bases remain unverified. Polymerase circles
the closed array repeatedly to produce HiFi sequence; ordered segmentation junctions split
that consensus molecule back into its eight full-length cDNA segments.
