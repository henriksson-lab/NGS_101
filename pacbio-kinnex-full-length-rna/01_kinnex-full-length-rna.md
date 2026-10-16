# PacBio Kinnex full-length RNA

## 1. What it is

Kinnex full-length RNA converts full-length cDNA molecules into ordered multi-insert arrays, then encloses each array in a SMRTbell for PacBio HiFi sequencing and computational segmentation.

## 2. Source read

🟢 Pacific Biosciences, *Preparing Kinnex libraries using the Kinnex full-length RNA
kit*, document 103-238-700 Rev09, August 2026.

🟢 PacBio, [SMRT Link v13.1 user guide](https://www.pacb.com/wp-content/uploads/SMRT-Link-User-Guide-v13.1.pdf), pp. 100–102: the full-length RNA default is MAS-Seq Adapter v3 (MAS8) for segmentation and Iso-Seq v2 Barcoded cDNA Primers for primer trimming. The guide prints twelve 5′ primer-trimming references and a 3′ trimming reference in its example. Page 100 explicitly states that the 3′ entry excludes poly(T) and is written in reverse-complement; it is not a physical PCR or RT oligo ordering declaration. 🟢 The official [Kinnex data-release README](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/README.txt) names the source dataset and analysis workflow.

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

🟢 The SMRT Link guide p. 101 and the linked [Iso-Seq v2 primer FASTA](https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta) print twelve barcoded 5′ trimming references and a 3′ trimming reference. For example, the bc01 5′ reference is `CTACACGACGCTCTTCCGATCTACTACACGCAATGAAGTCGCAGGGTTGGG`; the 3′ reference is `AAGCAGTGGTATCAACGCAGAGTAC`. SMRT Link p. 100 requires the 3′ reference to exclude poly(T) and be reverse-complemented. The selected bc01–bc12 depends on the sample. These positive read-trimming sequences do not establish the complete supplied cDNA amplification or RT oligos.

The declaration does not assign complete physical Kinnex orientation-PCR oligos,
segmentation-junction geometry, terminal-adapter arms/end geometry or a loaded
sequencing-primer binding sequence. The short analysis references above are retained
as positive source information without substituting them for ordering specifications.
The closed-array model explicitly leaves its insert–adapter junction unresolved; empty
overhang strings do not mean a source-known blunt end. Placeholder widths are schematic,
not assignments of physical sequence length. No publication-absence claim is accepted.

### Actual MAS-Seq predecessor and commercial version boundary

The application note's p. 2 explicitly attributes concatenation to Al’Khafaji et al.,
[High-throughput RNA isoform sequencing using programmed cDNA concatenation](https://doi.org/10.1038/s41587-023-01815-7).
The original article (PMC12236355, body and Methods “Multiplexed array assembly of
cDNA libraries”) describes a **15-member research assembly** with dU-containing PCR
oligos, USER digestion and HiFi Taq DNA ligase. Its original Supplementary Information
Table 1 (PDF pp. 25–26) positively prints the sixteen boundary identities A–P and the
complete 5′-10x research PCR oligos, including internal U notation. These are substantial
physical-oligo disclosures in the research predecessor. They do not establish that the
commercial eight-member Kinnex primer mixes or Kinnex enzyme have the same formulation.
No research USER, biotin-selection or ligase recipe is transferred into this declaration.

The research paper's Code/Data Availability links name Longbow and the MAS-Seq paper
repository. Their original READMEs describe computational annotation, segmentation,
WDL pipelines, result data and demonstrations with research 10/15-member arrays.
Current pinned READMEs are scope evidence, not proof of equivalence to the paper's
historical software versions. The reporting summary names Longbow 0.6.4; the paper-data
README separately describes preprint 0.2.2 and final 0.6.3, with 0.6.6 for UMI correction.
These statements do not identify commercial Kinnex physical components, and reading the
READMEs does not prove that an entire repository lacks some sequence.

The source checklist preserves both original reading passes and these source/version
boundaries. Unassigned commercial component details are model questions; no hypothetical
inaccessible component-disclosure document is introduced to explain them.

### Historical analysis originals named by the reporting summary

The exact historical Longbow **v0.6.4** README and configured models were inspected,
including `model.py` and `model_utils.py`. The official tag resolves to commit
`9844452445cc803f814613e203eaf971814babbc`, also named for the WDL workflow in the
reporting summary. Its configured MAS15/MAS10 identities are computational
recognition sequences; the workflow reads BAMs, annotates ordered motifs and extracts
cDNA reads. This checks the named historical source rather than assuming that a current
README describes that version.

The named long-read-pipelines commit
`82bfe460778e56ac642942e62638d594f22bd6bf` includes the MAS15/5′-10x analysis workflow
and a BAM UMI-tagging script. It also positively contains an earlier **eight-member
research construct**, `cDNA_array_8x.gb`: its LOCUS dates it to 24 April 2020 and its
FEATURES annotate 10x primer sites, modified bases and ligation regions. The companion
alignment FASTA includes `ohA`–`ohH`, poly(T) and TSO references. For example, its
`ohA_10x_Adapter` starts `ACGTACAGTCTACAC...`, whereas the commercial MAS8 reference A
is `AGCTTACTTGTGAAGAT`. Both are positive source designs; an eight-member count alone
cannot equate their reagents, input handles or junctions. The annotated research GenBank
record is not silently promoted to a commercial Kinnex ordering specification.

The reporting summary's exact paper-data commit
`0c166c63ce4d8ee37c835c89c9693bb65becb51c` was also inspected: its README identifies
research SIRV10/15 arrays and links analysis outputs, while the selected R script consumes
count matrices for normalization and clustering. These source/version checks classify the
actual named resources; they do not assert that entire repositories lack other sequences.
Commercial component identity and physical geometry remain unassigned in this declaration.
Direct external readout-reference and documentation children still require original-source
scope checks; the historical file reads do not complete the linked-source audit.

## 5. Read layout

The PacBio sequencing primer binds a terminal-adapter site whose bases remain unverified. Polymerase circles
the closed array repeatedly to produce HiFi sequence; ordered segmentation junctions split
that consensus molecule back into its eight full-length cDNA segments.
