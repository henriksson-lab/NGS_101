# Adapter and sequencing-primer inference audit

This is the work queue created after the MethylC-seq page was found to describe a named
NEXTflex single-index kit as an inferred generic dual-index TruSeq library. For every item,
the job is to trace the defining paper, supplement, detailed protocol and named historical
kit before retaining an `INFERRED` caveat. Exact product identity, index count and width,
adapter installation stage, PCR primers, sequencing primers and final molecular order must
be represented by the protocol model. Do not replace one unsupported generic construct with
another.

## Batch 1 — exact oligo-table RNA and targeted protocols

- [x] PRO-seq
- [x] eCLIP
- [x] Ribo-seq
- [x] smMIP
- [x] Earth Microbiome Project ITS

The repository already says these protocols have exact oligos, oligo tables or published
custom-primer extensions. Replace generic final-library shells where the primary sources
determine the actual construct.

## Batch 2 — bisulfite and chromatin protocols

- [x] scBS-seq
- [x] scRRBS
- [x] scDNase-seq
- [x] scMNase-seq
- [x] CUT&RUN
- [x] Strand-seq
- [x] ChIP-seq

Several existing research notes already contain complete adapter/PCR/read-primer sequences.
Promote source-supported sequence into the model; for named kits, resolve the exact historical
kit/version rather than assuming a current generic TruSeq endpoint.

### Batch 2 audit notes (2026-10-08)

- **scBS-seq:** Quail 2012 Supplementary Table 1, explicitly cited for iPCRTag,
  establishes the historical PE adapter, eight-base indexed primer and dedicated index
  primer. This is single-index, not modern dual-index TruSeq.
- **scRRBS:** the paper specifies standard premethylated indexed Illumina adapters;
  Illumina's adapter guide supplies the matching six-base single-index TruSeq set.
- **scDNase-seq:** the cited Hu 2013 method identifies Multiplexing Sample Prep
  Oligonucleotide Kit ref. 1005709. Its obsolete six-base single-index oligos are public.
- **scMNase-seq:** the paper names universal adapters and indexed primers but does not
  identify their bases or index set. The upstream-supported endpoint remains explicitly
  inferred; no dual-index claim is made.
- **CUT&RUN:** the KAPA library-preparation kit supplied reaction chemistry but not
  adapters. The paper never names the separately supplied adapter/index oligos, so the
  model now uses unavailable role tokens and does not invent a generic TruSeq endpoint.
- **Strand-seq:** Online Methods print the custom indexed PCR/index-read chemistry and
  name the obsolete Illumina PE oligos. The endpoint is a six-base single-index historical
  PE library.
- **ChIP-seq:** the contemporary protocol identifies the original Illumina Genomic DNA
  Sample Prep adapter and primers 1.1/2.1. The 2007 endpoint is unindexed single-read
  Solexa, not paired-end indexed TruSeq.

## Batch 3 — Hi-C family and named historical library kits

- [x] sci-Hi-C
- [x] sn-m3C-seq
- [x] in situ Hi-C
- [x] Micro-C
- [x] Nagano single-cell Hi-C
- [x] Stevens single-cell Hi-C
- [x] snHi-C
- [x] Dip-C

Preserve each protocol's internal contact/barcode topology. Audit the actual final library kit
and sequencing primers independently; do not collapse custom combinatorial indexes into a
generic dual-index library.

## Batch 4 — supplements and incomplete primary-source retrieval

- [ ] CEL-seq
- [x] MARS-seq
- [ ] Drop-seq
- [ ] inDrop v1
- [ ] SPLiT-seq
- [ ] microwell-seq
- [x] sci-ATAC-seq3
- [ ] scDamID
- [ ] LIANTI
- [ ] STRT-seq
- [ ] Small-seq
- [ ] VASA-seq
- [ ] txci-ATAC-seq
- [x] ch-ATAC-seq

Retrieve and read the primary supplement or detailed protocol. Secondary reconstructions may
guide searching but do not turn into green evidence merely because their sequences look
plausible.

### Batch 4 audit notes (2026-10-08)

Completed:

- **MARS-seq:** the original Jaitin 2014 supplementary PDF was already present as the
  byte-for-byte `jaitin-sm.pdf` mirror. Tables S7–S8 establish RT1, N3 + pool-barcode
  geometry, P5/P7 library primers and a paired-end/no-index-read endpoint. The schematic,
  notes and finder properties now follow that source.
- **sci-ATAC-seq3:** the original Table S7 workbook was already present as
  `scg_aba7612_domcke_table-s7.xlsx`. It establishes N5/N7 oligos, splints, PCR oligos and
  all four sequencing primers; the former secondary-source inference styling is removed.
- **CH-ATAC-seq:** the defining paper's original Supplementary Table 1 workbook was already
  present as `scg_CH-ATAC-seq_SupplementaryTable1.xlsx`. It establishes all three barcode
  sets, Tn5/ME oligos, HY oligos, MGI primers and modifications; the former inferred
  reconstruction is now source-backed.

Still blocked or incomplete after source audit:

- **CEL-seq:** the 2012 Cell Reports article/supplement remains inaccessible; CEL-seq2
  establishes its own successor chemistry but cannot retroactively establish the original
  library arms.
- **Drop-seq:** the original Table S6 supplement is gated and absent from the PMC Cloud
  bucket. The later Seq-Well table is not primary evidence for Drop-seq batch A/B oligos.
- **inDrop v1:** the Cell supplement remains gated and absent from the PMC Cloud bucket;
  the available Nature Protocols tables describe v2/v3, not v1.
- **SPLiT-seq:** the original oligo supplement remains unavailable; microSPLiT's later
  table is not evidence for the defining SPLiT-seq constructs.
- **microwell-seq:** the paper prints PCR primers but not the Tn5 adapter; the later
  same-lab reconstruction does not establish the defining adapter.
- **scDamID:** the defining oligo table remains unavailable; the exact adaptor fold is
  still secondary-source only.
- **LIANTI:** the author Supplementary DOCX is behind the PMC gate and absent from the
  Cloud bucket; the main paper establishes topology but not the exact transposon or
  NEBNext endpoint sequences.
- **STRT-seq:** the Genome Research barcode table and Nature Protocols details remain
  gated. Later STRT-C1/2i supplements cannot silently define the original 2011 protocol.
- **Small-seq:** Supplementary Table 1 containing SRX remains unavailable, so its exact
  index oligo and sequencing-primer identity cannot be promoted.
- **VASA-seq:** the defining supplement prints VASA-specific indexes but cites SORT-seq,
  inDrop v3 and TruSeq Small RNA for several capture/library oligos. Those cited primary
  protocols must be resolved before replacing the current role placeholders.
- **txci-ATAC-seq:** the paper's custom Tn5 and i7 oligos are established, but the 10x
  bead/P5 side is only reconstructed by the secondary page; the cited protocols.io/10x
  source is still unavailable.

## Cross-cutting shared-helper audit

Audit every consumer of `batch_ngs.truseq_library()`. The helper currently defaults to a
generic dual-index endpoint and places full P5 and full Read 1 sequences consecutively even
though their `ACAC` boundary overlaps in a physical TruSeq-family oligo. Confirm index count,
index width, primer family and overlap for each protocol before using the helper:

- [ ] 10x Multiome
- [ ] 10x Flex
- [ ] Visium
- [ ] CITE-seq
- [ ] ECCITE-seq
- [x] CUT&RUN
- [ ] Duplex Sequencing
- [x] eCLIP
- [x] PRO-seq
- [x] Ribo-seq
- [ ] SHAPE-MaP
- [ ] SHARE-seq
- [ ] Slide-seqV2

Prefer a constructor that encodes the resolved physical adapter family and validates its own
declared sequencing primers. Do not use a test to compensate for an ambiguous constructor.

## Caveats that are probably legitimate

Do not remove uncertainty merely to reduce the number of dotted regions. Commercially
undisclosed sequence in current 10x internals, BD Rhapsody, PIP-seq, Oxford Nanopore, PacBio,
DNBSEQ, Flex, SureSelect and similar kits should remain a role/length placeholder unless an
authoritative vendor document actually publishes the bases. Variable target arms, sample
indexes and barcode identities also remain placeholders even when their surrounding fixed
sequence is known.

## Completion standard for each checkbox

1. Record primary sources read and any unavailable source in the protocol notes/manifest.
2. Model the actual adapter and final library, including shared/overlapping sequence once.
3. Declare and computationally locate every sequencing primer used by the run.
4. State explicitly when the library is single-index, dual-index, inline-indexed or requires
   custom sequencing primers.
5. Remove `INFERRED` only from regions directly established by authoritative evidence.
6. Update catalogue properties and the checkbox here.
7. Build the affected page and run its focused checks; leave generated HTML uncommitted.
