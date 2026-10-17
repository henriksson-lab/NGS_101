# slide-DNA-seq — source notes

Defining paper: Zhao et al., *Nature* (2022), [doi:10.1038/s41586-021-04217-4](https://doi.org/10.1038/s41586-021-04217-4).

- 🟢 PMC manuscript, Methods and Supplementary Table 1 oligos were read. Supplementary Table 1 is the `MOESM3` workbook (sheet 1, 15 rows): bead oligo, bridge oligo, mosaic end, bulk and slide-DNA-seq PCR forward primers, PCR reverse primer, the four custom sequencing primers, the two Tn5 adapters and the modification key.
- 🟢 A pre-sequenced 10-µm bead array supplies spatial coordinates. HCl removes histones; custom Tn5 adapters tag tissue DNA.
- 🟢 Bead barcodes are photocleaved and ligated to proximal genomic fragments before PCR.
- 🟢 The bead oligo contains a photocleavable linker (`PCT`) and two 7-nt
  variable blocks around a fixed linker. Supplementary Table 1 also prints
  the bridge oligo, 5′-phosphorylated/3′-dideoxy blocked mosaic bottom,
  custom Tn5 adapters, indexed PCR primers, and custom Read 1, Read 2,
  index 1 and index 2 sequencing primers. The page now lists these
  additional exact oligos and the published modification notation.
- 🟢 Supplementary Methods §1.5, lines 100–117: 0.5% Triton-X permeabilization,
  0.1 N HCl 5 min, tagmentation at 37 °C for 1 h (repeated four times in the
  increased-tagmentation variant), bridge oligo 10 µM hybridized at 37 °C for 2 h
  in 4× SSC, 365 nm photocleavage for 3 min, T4 ligase 12 U/µL for 1 h at room
  temperature, overnight Proteinase K at 55 °C, NucleoSpin cleanup, Buenrostro
  2015 PCR with 65 °C annealing and a 0.9× AMPure cleanup.
- 🟡 The bridge oligo is a **splint**, derived by pairing the printed strands:
  `CGTTAGCAAAGCTC` is complementary to the 5′-phosphorylated Tn5 adapter 1 handle
  and `CGCCCATGTGTC` to the bead oligo's 3′ handle `GCTCGGACACATGGGCG`, so T4
  ligase seals the bead 3′ OH to the adapter 5′ phosphate. The ligation geometry
  is now declared splinted instead of unavailable; the source states the
  hybridization and ligation but does not draw the junction.
- 🟢🟡 **The index-primer binding sites are now resolved** (previously an open model
  limitation). The final strand is assembled from the printed oligos rather than
  retyped, which exposed two missing pieces in the earlier model: the Tn5 adapter 1
  handle begins `GAG`CTTTGCTAACGGTCGAG while the Read 1 *primer* begins
  `GCTTTG…`, and the PCR reverse primer contributes twelve bases
  (`GTTGGCACCAGG`) beyond the Tn5 adapter 2 arm. With both restored:
  index 1 (`GCTGGTGCAACATCCGTAAGCCTGGTGCCAAC`) is the exact reverse complement of
  the PCR reverse primer's 32-base arm and reads the 8-base i7 sample index; index 2
  (`CTCGACCGTTAGCAAAGCTCCGCCCATGTGTC`) matches the bead-handle/Tn5-adapter-1
  junction exactly and reads **backwards across the spatial barcode**, which is why
  its read is 44 cycles rather than 8: 5 bases of the bead handle + block B (7) +
  linker (18) + block A (7) + 7 bases of the Read 1 arm. Index 2 is therefore a
  spatial-barcode read, not a sample-index read, and the page now says so.
- 🟢 Supplementary Methods §1.5, lines 113–117, specifies the slide-DNA-seq
  NextSeq/NovaSeq recipe: Read 1 57 bp, Read 2 57 bp, Index 1 8 bp, Index 2
  44 bp. Section 1.6, lines 129–131, gives a separate 40/40/8/8 recipe for
  bulk whole-genome controls, which also use the **bulk** PCR forward primer
  (`…GACACATGGGCGGAGCTTTGCTAACGGTCGAG`, priming across the same junction the
  slide-DNA-seq index 2 primer uses). These must not be conflated.
- 🔴 Open: the i5/i7 sample-index bases themselves are not published, the PCR
  cycle number is left to the Buenrostro 2015 qPCR-guided procedure, and the
  bead oligo's `TTTT-PCT` 5′ spacer is drawn only as a modification note.
- 🟢 The reporting summary links the bead-indexing code to the authors'
  `broadchenf/Slideseq` repository and the downstream analysis code to
  `buenrostrolab/slide_dna_seq_analysis`. Its sample/library statistics table
  contains no oligo or reaction recipe. The linked frozen-tissue nuclei
  extraction protocol is for the companion snRNA-seq assay, not the
  slide-DNA-seq DNA library. Zhao cites the February 2020 v1
  (`10.17504/protocols.io.bck6iuze`); the initially fetched
  `10.17504/protocols.io.bi62khge` is a July 2020 v2.
- 🟢 The paper-linked [Zenodo v1.0 code archive](https://doi.org/10.5281/zenodo.5553305)
  contains `scripts/preprocess/extractBCfromR2.py` and
  `scripts/run/slidedna_alignment_template.sh`. The scripts use three
  demultiplexed FASTQ streams named R1/R2/R3, extract two 7-nt spatial
  barcode blocks from the middle stream, attach that barcode to the two
  genomic-read headers, and align those genomic reads as a pair. The
  scripts support barcode handling; they do not resolve the physical
  index-primer binding sites missing from this page model.
- 🟢 The cited Buenrostro 2015 ATAC-seq PCR Methods (PMC4374986,
  "PCR Amplification," lines 111–167) describe 72 °C end extension,
  five initial cycles, a qPCR side reaction to choose additional cycles,
  and 63 °C annealing. Zhao Supplementary Methods §1.5 explicitly changes
  the annealing temperature to 65 °C for slide-DNA-seq. This is a
  source-declared adaptation, not a discrepancy between the papers.
- 🟢 Rodriques 2019 Slide-seq supplementary Methods print RNA-capture
  bead oligos with split-pool barcode blocks and a poly(T) tail. Zhao
  Supplementary Methods §1.3 says a different oligo was coupled to
  slide-DNA-seq beads for DNA ligation; Zhao Supplementary Table 1 is
  therefore the authority for that DNA-capture bead sequence.
