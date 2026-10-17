# DM-tRNA-seq / DM-TGIRT-seq

**Defining source:** Zheng *et al.*, *Nature Methods* (2015),
[doi:10.1038/nmeth.3478](https://doi.org/10.1038/nmeth.3478).

## Source-led chemistry

- 🟢 The defining author-manuscript Methods, “TGIRT template-switching reactions
  and sequencing” (JATS text line 90), print the RNA starter
  `AGAUCGGAAGAGCACACGUCUAGUUCUACAGUCCGACGAUC/3SpC3/`. Its 3′ spacer blocks
  extension. A complementary **32P-labeled DNA primer** has a single 3′ `T`
  overhang; “32P” denotes the isotope, not a DNA-primer length.
- 🟢 That method names **GsI-IIC MalE rigid fusion RT**, followed by NaOH treatment,
  neutralization, gel purification of cDNA, CircLigase II circularization and
  Phusion-HF amplification. The RNA starter must therefore be distinguished from
  the extended DNA strand and the purified single-stranded DNA product.
- 🟡 The page reconstructs the minimal DNA starter as the reverse complement of
  the printed RNA plus its specified 3′ `T`. The complementary core is marked
  inferred: full-length pairing without additional DNA extensions is a modelling
  assumption. The label position is unassigned.
- 🟢 Both PCR primers are printed in the defining Methods. Their opposed binding
  suffixes account for the entire RNA starter; the remaining PCR-primer bases are
  5′ tails. A full PCR primer need not occur as a contiguous sequence in the starter.
- 🟡 Circular-to-linear PCR product assembly, selected barcode, barcode width and
  loaded HiSeq sequencing-primer reagent remain unassigned **model data**. The page
  does not assert that all linked sources omit these facts.

## Source versions and linked material

- The cited Mohr *et al.* 2013 original (doi:10.1261/rna.039743.113), body and
  relevant Methods, distinguish GsI-IIC-MRF and TeI4c-MRF enzymes. Its SOLiD
  template-switching example is a separate library architecture.
- The cited Shen *et al.* 2015 original supplement, “TGIRT-seq of RQC-associated
  RNAs” (CDC Stacks original DOCX text lines 72–75), positively describes a
  complementary 42-nt DNA starter, equimolar single-base overhangs, and
  **TeI4c-MRF with a C-terminal truncation**. The RNA and PCR sequences agree with
  the DM source, but its enzyme and overhang choices do not replace DM’s GsI-IIC
  and selected `T`.
- The cited Katibah *et al.* 2014 original (doi:10.1073/pnas.1412842111), main
  Methods, points to SI Methods for its Illumina library modifications. The actual
  SI original remains a pending source obligation: publisher HTTP and HTTPS
  supplement links returned 403, while the PMC attachment route returned a
  download gate. Author-copy recovery remains to be checked.
- The later author-lab small-RNA/CircLigase PDF (filename dated 2016-06-14) and
  [InGex TGIRT-III protocol](https://fnkprddata.blob.core.windows.net/domestic/data/datasheet/IGX/TGIRT10.pdf)
  (updated 2017-01-23) print a DNA primer with a degenerate `N` overhang and
  explicitly describe 5′ labeling. These are later variants; neither assigns
  the 2015 study’s label position or selected run reagent.
- The original Epicentre CircLigase II manual, Lit. 298, Rev. A, December 2012,
  pp. 2–4, requires single-stranded DNA with 5′ monophosphate and 3′ hydroxyl.
  The 2015 source explicitly extends the manufacturer incubation. The manual’s
  optional exonuclease cleanup is not imported into the DM graph.
- The current official HiSeq 2500 primer-box page lists alternative primer mixes
  and links compatibility/instrument guides. It is a candidate comparator, not
  evidence that the authors used a HiSeq 2500 or any particular mix.

The source review remains **in progress**. Publisher supplements, the separate
NIHMS results supplement and later vendor method retain their earlier independent
checks. Newly recovered originals and direct linked-source dispositions are recorded
in the maintainer source checklist; no missing-information claim is accepted.

## Boundary

This is the 2015 demethylase/TGIRT protocol. It is not mim-tRNAseq or YAMAT-seq.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC4624326 JATS84–94, demethylation and TGIRT template-switching.

Specified 2h room-temperature wt/D135S AlkB treatment, CircLigaseII60°C5h and 12-cycle PhusionHF PCR. Confirmed printed 41nt blocked RNA starter; reconstructed DNA remains inferred.

Deferred: Source DNA starter identity and radiolabel position, original supplementary material, deacylation provenance and read acquisition; source total-RNA D135S molar wording not silently corrected. This pass does not establish a complete source review or publication-wide absence.

### Published D135S dose inconsistency

The published Methods, “Demethylation reactions”, estimates about 40 pmol tRNA
in the total-RNA input and gives D135S as both fourfold and 200 pmol
(PMC4624326 JATS86; PDF278–280). Four times 40 pmol is 160 pmol.
The intended amount or ratio cannot be selected from this passage; the model
leaves dose unassigned. This notice is also shown on the schematic page.
