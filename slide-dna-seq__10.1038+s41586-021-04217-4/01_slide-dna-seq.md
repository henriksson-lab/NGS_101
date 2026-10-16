# slide-DNA-seq — source notes

Defining paper: Zhao et al., *Nature* (2022), [doi:10.1038/s41586-021-04217-4](https://doi.org/10.1038/s41586-021-04217-4).

- 🟢 PMC manuscript, Methods and Supplementary Table 1 oligos were read.
- 🟢 A pre-sequenced 10-µm bead array supplies spatial coordinates. HCl removes histones; custom Tn5 adapters tag tissue DNA.
- 🟢 Bead barcodes are photocleaved and ligated to proximal genomic fragments before PCR.
- 🟢 The bead oligo contains a photocleavable linker (`PCT`) and two 7-nt
  variable blocks around a fixed linker. Supplementary Table 1 also prints
  the bridge oligo, 5′-phosphorylated/3′-dideoxy blocked mosaic bottom,
  custom Tn5 adapters, indexed PCR primers, and custom Read 1, Read 2,
  index 1 and index 2 sequencing primers. The page now lists these
  additional exact oligos and the published modification notation.
- 🟡 The final page places the two printed genomic read sites exactly. Both
  index-primer sequences are now displayed, but their binding sites have not
  been reconstructed on the representative final strand; this is a model
  limitation rather than missing source bases.
- 🟢 Supplementary Methods §1.5, lines 113–117, specifies the slide-DNA-seq
  NextSeq/NovaSeq recipe: Read 1 57 bp, Read 2 57 bp, Index 1 8 bp, Index 2
  44 bp. Section 1.6, lines 129–131, gives a separate 40/40/8/8 recipe for
  bulk whole-genome controls. These must not be conflated.
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
