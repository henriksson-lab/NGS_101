# ChEC-seq

## Current native graph — batch08 (2026-10-18)

The current graph describes the defining ChEC preparation. DNA, RNA-derived material, protein-derived material and GlycoBlue are independent populations. The defining Methods require phenol/chloroform extraction, ethanol recovery, RNase treatment, negative 2.5:1 AMPure selection retaining the small-DNA supernatant, and a second extraction/precipitation before library preparation. RNase hydrolysis does not itself physically discard RNA-derived material.

The defining paper cites the 2011 and 2013 library methods and positively names KAPA amplification. Their detailed historical assembly is outside this capped pass. That preparation is an explicit boundary followed by PCR on its actual native precursor, with original templates/carrier/RNA-derived recovery unresolved. The later Orsi Basic Protocol3 oligos remain a separately labelled source reference, not the defining study's reagents. Defining paired-end25 acquisition is retained; the later reference's six-base physical index does not assign historical index cycles.

Source scope: Defining published PMC4618392 JATS Methods ChEC and Library preparation/sequencing129–143; published SI opening FigS1 digestion caption/table inventory; Orsi2015 PMC4410783 BasicProtocol3 sourceoligos812–832 and steps4–39/text950–1205. Mandatory scraper adapter/PCR hits checked directly in the Orsi original.

This is a capped migration pass, not full source review or accepted missing information.

## Earlier source audit and model description (historical)

The following records preserve earlier source readings and the former model. Their Orsi graph/endpoint descriptions are superseded by the native graph above.

Zentner et al., *Nature Communications* (2015),
[doi:10.1038/ncomms9733](https://doi.org/10.1038/ncomms9733).

## Sources read

- 🟢 The defining article and supplementary information (PMCID: PMC4618392).
- 🟢 Zentner et al., *Genome-wide Mapping of Protein-DNA Interactions with ChEC-seq
  in Saccharomyces cerevisiae* (PMCID: PMC5608238), for the detailed experimental
  sequence and controls.
- 🟢 Orsi et al., *Mapping regulatory factors by immunoprecipitation from native
  chromatin* (PMCID: PMC4410783), Basic Protocol 3, for the exact paired-end adapter
  and enrichment-primer oligos in a separately documented Henikoff-laboratory
  low-input library variant. The defining ChEC paper does not assign a concrete
  sample-index member to each published library.
- 🟢 GEO series GSE67453, including GSM1647311, for the extraction, size-selection and
  25-base paired-end run metadata of the defining data set.
- 🟢 Its deposited **ChEC-seq Protocol**, dated March 2015, confirms the unbound
  fraction selection and the two library-method PMIDs with the KAPA exception.
  The entire small RTF was independently read twice after registered recovery.

## Published corrections and source boundaries

- 🟢 The [2015 corrigendum](https://doi.org/10.1038/ncomms10264), sole correction
  paragraph, states that columns E, F, G, H and K of the original Supplementary
  Data 1 **Reb1** sheet were sorted independently. The authors replaced the table
  and state that the incorrect version was not used for the reported analyses.
- 🟢 The [2017 corrigendum](https://doi.org/10.1038/ncomms15723), correction
  paragraphs, identifies mismatched motif-centering of the Fig. 7 random controls.
  The authors withdraw the conclusion that this analysis establishes
  shape-dependent, sequence-independent binding at low-scoring sites. They state
  that this does not invalidate ChEC-seq for identifying binding sites. The public
  schematic renders both corrections; neither supplies a different library adapter.
- 🟢 Defining Methods **Library preparation and sequencing** (JATS line 135)
  cites Henikoff 2011 (reference 46, PMC3215028) and Zentner/Henikoff 2013
  (reference 58, PMC3889552). Orsi 2015 is a separate same-laboratory protocol,
  linked by the later detailed ChEC procedure, and must not silently substitute
  for those two directly cited implementations. Their original comparison remains
  open in the source checklist.
- 🟢 Actual GSE67453 family metadata contains 56 ChEC samples with the same
  extraction/library description and HiSeq 2500, plus a separate ChIP comparator
  using HiSeq 2000. The ChEC metadata confirms 25-base paired-end mapping and the
  two directly cited library methods. It does not by itself assign an i7 member
  or an acquired technical-index cycle recipe.

## Molecular logic

1. 🟢 A yeast chromatin protein is genetically fused to 3×FLAG-MNase. Nuclear calcium is
   normally too low for strong MNase activity, so the fusion can bind chromatin before
   the timed cleavage reaction.
2. 🟢 Cells are washed and permeabilized for 5 min at 30 °C in 0.1% digitonin. Calcium is
   then added to 2 mM, activating tethered MNase in situ. Aliquots are withdrawn across
   short time points and stopped with EDTA, EGTA, high salt and SDS.
3. 🟢 Rapid cuts are enriched around sites occupied by the fused factor. Two nearby cuts
   release a small DNA fragment containing the protected binding region. Later cleavage
   increasingly includes transient encounters and bulk accessibility.
4. 🟢 Nucleic acid is phenol/chloroform-extracted, ethanol-precipitated and RNase-treated.
   A 2.5:1 AMPure-bead:sample ratio binds large genomic DNA; the unbound fraction carrying
   small ChEC fragments is recovered for library preparation.
5. 🟢 The separately illustrated Orsi library variant combines T4 DNA polymerase,
   Klenow and T4 PNK in one end-polish incubation, cleans up, then adds 3′ dA with
   Klenow exo-minus and cleans up again. Its protected Y adapters are ligated,
   bead-bound ligated DNA is retained, KAPA PCR follows, and bead-bound amplified
   DNA is retained. These positive library cleanups retain the opposite fraction
   from the initial ChEC negative selection. The original 2015 experiment names
   the 2011/2013 library methods; the displayed Orsi member is a scoped
   illustration while that full historical comparison remains open.

## Exact library oligos

The separate Orsi low-input procedure publishes the ordering forms, including modifications:

- Universal Adapter:
  `5′-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′`
- Indexed Adapter:
  `5′-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTT*G-3′`
- P5 enrichment primer: `5′-AATGATACGGCGACCACCGAG-3′`
- P7 enrichment primer: `5′-CAAGCAGAAGACGGCATACGAG-3′`

🟢 `*` is the terminal phosphorothioate linkage, `/5Phos/` is the indexed adapter's
5′ phosphate, and `NNNNNN` is the six-base i7 sample barcode in this illustrated
variant, which has no i5 index. The model preserves the complete 58-base Universal
Adapter, including its shared P5/Read 1 overlap once. The former shared endpoint
duplicated four bases at that boundary; this was a repository error.

🟡 The illustration uses six Index 1 cycles to display that variant's barcode.
This is not a recovered original 2015 acquisition mask or selected barcode member.
The initial physical Y-arm PCR template is not yet fully constructed; its missing
model relation is explicit and does not imply that the source omits the primers.

## Read interpretation and controls

🟢 The defining data set used 25-base paired-end reads on HiSeq 2500. Both mate starts
mark repaired MNase boundaries around the selected short fragment. The paired interval,
not merely local coverage, therefore represents the two-cut ChEC product.

🟢 A soluble nuclear MNase strain under the matched promoter measures MNase sequence bias
and accessible-chromatin background. No-calcium and early-time samples distinguish basal
damage and rapid tether-proximal cleavage. ChEC-seq is antibody-free and unfixed, but it
requires an expressed fusion whose residence and linker geometry may affect cleavage.

## Bounded curation round 2 (2026-10-16)

Inspected: Published original PMC4618392 JATS119–135.

Added proteinase K deproteinization before extraction. Recorded2.5:1 negative AMPure selection retains short-DNA supernatant, followed by repeated extraction/precipitation; existing model agrees on25-cycle paired reads.

Deferred: Cited original library procedures and oligos; second extraction/precipitation graph state and actual sample-index acquisition; detailed kinetic control validation. This capped pass does not establish a complete source review or publication-wide absence.
