# Methyl-HiC — source notes

Defining paper: Li et al., *Nature Methods* (2019), [doi:10.1038/s41592-019-0502-z](https://doi.org/10.1038/s41592-019-0502-z).

## Sources read

- 🟢 PMC full text, main description and Methods.
- 🟢 Original Supplementary Protocol, including the reagent list, bulk on-bead
  adapter-ligation procedure and complete scMethyl-HiC oligo table.
- 🟢 The oligo scraper was run on the full Li source bundle; its P5L/P7L hits map to the separately headed single-cell oligo table, not the bulk adapter.
- 🟢 Two independent readers checked the [Rao et al. in situ Hi-C predecessor](https://doi.org/10.1016/j.cell.2014.11.021) original Extended Experimental Procedures, section I.a.1 step 54 (local extracted lines 6771–6775). It names an indexed Illumina adapter and says to record the sample-index combination, but prints no specific member. The linked ENCODE transcription agrees at its lines 227–261. Rao step 60 refers to an Illumina 2007 protocol whose version-matched text is still unexamined.
- 🟡 The [Illumina TruSeq DNA guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseqdna/TruSeq_DNA_SamplePrep_Guide_15026486_C.pdf) lists multiple indexed adapter families. The Li bulk method does not choose one. The Zymo D5005 conversion guide is linked for the named bisulfite reagent, not for adapter identity. A second independent source and link check is pending; the detailed bulk source checklist is retained in `catalogue/source_reviews/` for maintainers.

## Bulk Methyl-HiC

- 🟢 Fixed nuclei are digested overnight with DpnII. Ends are filled with biotinylated nucleotides and ligated in situ.
- 🟢 After crosslink reversal the DNA is sonicated to about 400 bp and pulled down with streptavidin beads.
- 🟢 Adapters are ligated on beads, then the molecules undergo EZ DNA Methylation-Gold bisulfite conversion and KAPA HiFi HotStart Uracil+ amplification.
- 🟢 Libraries were sequenced paired-end on HiSeq 2500 or 4000 with balanced libraries or PhiX because of bisulfite base imbalance.
- 🟢 The bulk supplementary protocol explicitly lists `TruSeq PCR Primers`, ligates an
  indexed `Illumina Adapter` on beads, and then performs uracil-tolerant PCR. Standard
  TruSeq Read 1, Index 1 and Read 2 entry sites can be located on an inferred compatible
  final shell. The article's explicit cytosine-methylated TruSeq adapter belongs to its
  separate whole-genome bisulfite sequencing control, not the bulk Methyl-HiC method.
- 🟡 The exact ligated adapter strands and selected index member are unresolved. The page
  marks the canonical TruSeq-compatible bases as inferred rather than pretending they were
  transcribed from the checked Li article, supplement or Rao predecessor; the ligation-ready
  end remains unresolved while the independent source check is pending.

## Related single-cell branch

- 🟢 scMethyl-HiC sorts nuclei after proximity ligation, bisulfite-converts each nucleus,
  uses indexed P5 random priming, Exonuclease I/SAP cleanup, and Adaptase addition of P7
  before dual-index PCR. The supplement prints all four P5-tailed random-primer families,
  eight custom P5 indexing primers and twelve custom P7 indexing primers.
- 🟢 The separately documented [scMethyl-HiC protocol](../scmethyl-hic__10.1038+s41592-019-0502-z/scmethyl-hic.html) models this branch, including its published random-primer and indexing-primer families. This page models only the bulk library.
