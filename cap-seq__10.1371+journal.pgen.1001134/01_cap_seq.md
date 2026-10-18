# CAP-seq — CXXC affinity purification sequencing

## Current native graph scope — 2026-10-18 batch11

Two originals were inspected: published defining PMC2944787 JATS, DNAchromatography
and library Methods218–226, and the directly cited Skene2010 author-manuscript
PMC4338610, Methods109–121. The predecessor points to SupplementalData; those
detailed preparation sections were not inspected. Mandatory scraper:35hits in2files.
No additional originals were opened after the capped investigation.

🟢 Solexa adapter ligation precedes CAP purification. DNaseI treats the recombinant
CXXC reagent preparation, not the sample genomicDNA. CXXC is immobilized on nickel
Sepharose; sampleDNA binds in0.1MNaCl, is washed at600mM for the selected human
route, eluted in1MNaCl, pooled, concentrated and precipitated. Mouse560mM is a
separate source condition. Affinity enriches native unmethylatedCpG clusters; it
does not change a methylatedC into unmethylatedC, nor prove population purity.

🟡 The native graph owns independent originalDNA/opposite-strand/residual classes
and adapter roles, without the old shared endpoint/affinity recipes. Exact input
preparation and postCAP Solexa construction remain scoped predecessor boundaries;
no PCR/size-selection order or canonical adapter/index shell is invented. Single-end
sequencing and mapping-score≥30 are separate physical/computational stages.

Deferred: cited Skene2010 SupplementalData preparation/library sections and any
applicable native adapter/primer/end/run identities; exact class-specific retention.
Historical audits remain history. This capped pass is not a full source review or
accepted missing-information claim.


## Name boundary

This page covers **CXXC Affinity Purification plus sequencing**, the named DNA method introduced by Illingworth and colleagues. It is not the similarly named capped-RNA “Cap-seq”/CIP–TAP workflow; that is a distinct protocol and should receive its own page if added.

## Sources read

- 🟢 Illingworth *et al.*, “Orphan CpG islands identify numerous conserved promoters in the mammalian genome,” *PLOS Genetics* (2010), DOI `10.1371/journal.pgen.1001134`, PMCID `PMC2944787` (defining paper and methods).
- 🟢 Earlier MAP/CAP preparation methods cited by the defining paper were used only to interpret the named chromatography steps; no unprinted oligo was imported from them.

- Round 2026-10-17-breadth-first-07 re-read PMC2944787 JATS Methods lines 218 (CAP) and
  224–226 (library preparation and Solexa sequencing), Results line 79 (600→560 mM wash
  optimisation), Table S1 (s010) header rows, and the predecessor PMC4338610 JATS lines
  113–135 (its Supplemental Data pointer; `NIHMS62106-supplement-01.pdf` is listed in the
  XML but is absent from PMC Cloud and the PMC page returns HTML, so it stays unread).
  The cached bioRxiv `2024.09.10.612220_v3` file in the source bundle is an unrelated
  2025 long-read microbiome method and is not used.

## Molecular interpretation

- 🟢 CAP enriches genomic fragments containing clusters of **unmethylated CpGs** by binding them to recombinant CXXC immobilized on nickel-charged Sepharose.
- 🟢 The sequencing adapters are ligated to input DNA **before** affinity purification. This matters: selection acts on an already sequenceable library rather than on naked DNA followed by a separate library build.
- 🟢 Binding occurs in column buffer containing 0.1 M NaCl. Human libraries are washed at 600 mM NaCl; the optimized mouse workflow uses 560 mM because mouse CpG islands are, on average, less CpG-dense. Bound material is eluted at 1 M NaCl.
- 🟢 Recombinant CXXC (expressed as in reference 7) was pretreated with DNaseI to remove bacterial DNA and bound to nickel-charged Sepharose at 50 mg/ml bead volume (Methods line 218).
- 🟢 Input DNA was "prepared as for MAP-seq" (reference 57). The page shows this as a warned delegated step; fragmentation and fragment ends are not assigned.
- 🟢 Eluted fractions were pooled, concentrated and precipitated as for MAP; CAP was performed once per sample, with at least two technical replicates per DNA sample. Whether any amplification or size selection followed elution is delegated to reference 57 and not drawn.
- 🟢 The enriched libraries were sequenced single-end on the Solexa platform at the Wellcome Trust Sanger Institute; reads were mapped with MAQ and retained at mapping score ≥30 (a computational filter). Read length is not stated in the inspected Methods.
- 🔴 The defining paper's “Library preparation and Illumina Solexa sequencing”
  section delegates the library recipe to its reference 57, Skene *et al.*,
  *Molecular Cell* (2010), PMCID `PMC4338610`. That open author manuscript
  itself defers detailed Solexa protocols to Supplemental Data, which were
  absent from the PMC Cloud archive checked here. Neither accessible article
  prints CAP-seq's ordered adapters or sequencing primer, so the diagram
  labels these ends by role. A contemporary generic Solexa adapter would be
  an inference, not a transcription of this experiment.

## Final readout

A mapped read marks an adapter-ligated genomic fragment retained because its clustered unmethylated CpGs bound CXXC under the chosen salt conditions. CAP-seq reports affinity enrichment, not a base-by-base conversion call.
