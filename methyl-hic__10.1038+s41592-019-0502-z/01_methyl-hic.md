# Methyl-HiC — source notes

Defining paper: Li et al., *Nature Methods* (2019), [doi:10.1038/s41592-019-0502-z](https://doi.org/10.1038/s41592-019-0502-z).

## Sources read

- 🟢 PMC full text, main description and Methods.
- 🟢 Original Supplementary Protocol, including the reagent list, bulk on-bead
  adapter-ligation procedure and complete scMethyl-HiC oligo table.
- 🟢 The oligo scraper was run on the full Li source bundle; its P5L/P7L hits map to the separately headed single-cell oligo table, not the bulk adapter.
- 🟢 Two independent readers checked the [Rao et al. in situ Hi-C predecessor](https://doi.org/10.1016/j.cell.2014.11.021) original Extended Experimental Procedures, section I.a.1 step 54 (local extracted lines 6771–6775). It names an indexed Illumina adapter and says to record the sample-index combination, but prints no specific member. The linked ENCODE transcription agrees at its lines 227–261. Rao step 60 cites an Illumina 2007 protocol for PCR. Two readers checked a [period Illumina ChIP guide](https://support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/cluster_station/11257047_ChIP_Sample_Prep.pdf): its adapter mix is unindexed, so it cannot identify Rao's indexed adapter or Li's selected bulk adapter.
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

## Bounded source pass — 2026-10-17

Read published PMC6765429 JATS Methods79–105: bulk Methyl-HiC, separate
single-cell and WGBS branches, sequencing and mapping. Inspected published
MOESM2 reporting summary identity/sections1–59 and the original published
Supplementary Protocol `supplement-5.pdf.txt` reagent list1–57, shared
fill/ligation191–209 and bulk211–257,304–347,374–415,437–485. Manifest
identifies all three originals. Ran the full bundle oligo scraper; its hits
were assigned to the separately headed single-cell oligo table, not used
as bulk adapter evidence.

🟢 The graph now includes proteinase/reverse crosslinking, ethanol recovery,
post-sonication SPRI cleanup, on-bead end repair and dA tailing, Quick Ligase,
same-adapter-ligated lambda control before conversion, and final double AMPure
selection. The optional capture-supernatant background control is retained
as a stored branch rather than declared discarded. Its subsequent library
preparation is outside this graph.

🟡 Source handling differences remain public: shared fill-in step1 names
MboI although the main Methods and reagent list specify DpnII. Both bulk
and single-cell pages retain DpnII and display the unresolved label difference.
The bulk protocol sends converted supernatant to a purification column and
then says PCR directly on beads; the main article describes purified converted
DNA. The graph follows that purified-DNA route through a warned unresolved
recovery boundary, rather than inventing the bead/column mechanism.

🟡 Converted intermediates no longer display the unchanged contact state.
The final insert uses inferred C/T classes for cytosines after conversion
and PCR; methylation values remain unassigned. The canonical adapter shell
is still inferred. Source-specific primer/index identity is not supplied by
the single-cell oligo table or separate WGBS control.

Deferred: detailed source tables/figures and uninspected protocol parts;
Rao/vendor originals; explicit crosslink/precipitation/repair and conversion
recovery molecules; whether the protocol's 0.4 nM biotin stock wording is
intended; exact bulk adapters/indices and acquired run lengths. The shearing
recipe's “For PE2x100” is a preparation note, not independently sufficient
for a study-wide acquired cycle declaration. Counterpart single-cell page
received only the shared-source discrepancy notice, not its pending full
bounded pass. No absence claim is accepted or complete validation asserted.
