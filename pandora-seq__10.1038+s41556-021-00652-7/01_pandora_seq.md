# PANDORA-seq

**Defining source:** Shi *et al.*, *Nature Cell Biology* (2021),
doi:10.1038/s41556-021-00652-7. The PMC full text, accessible supplement,
and the [NEB E7330 instruction manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7300_e7330_e7560_e7580.pdf)
(version 9.0, July 2025) were read.

## Evidence used

- 🟢 Methods establish 15–50-nt size selection before treatment, T4 PNK end repair,
  AlkB demethylation and subsequent NEBNext small-RNA library preparation.
- 🟢 Methods line 319 prints a **Takara** preadenylated, 3′-blocked SR adapter for a
  ligation-validation assay. Methods line 331 separately identifies **NEBNext E7330S**
  for the sequenced libraries; it names the 3′ adapter, RT primer, 5′ adapter and PCR
  primer cocktail in reaction order. The assay adapter alone does not establish the
  library adapter's supplier or chemistry.
- 🟢 The NEB E7330 manual, p. 20, prints the complete 3′ and 5′ adapter, RT primer,
  SR PCR primer and index 1 PCR primer sequences, including rApp, 3′ amine, RNA bases
  and 3′ phosphorothioate bonds. The 3′ adapter has the same bases as the Takara assay
  oligo. [NEB's sequencing FAQ](https://www.neb.com/en/faqs/do-i-have-to-use-custom-sequencing-primers-for-my-libraries-generated-with-nebnext-small-rna-library-prep-kit-for-illumina)
  states that these libraries use standard Illumina sequencing reagents, with no
  custom sequencing primer required.
- 🟡 The available vendor manual is dated 2025. It establishes the published E7330
  kit formulation but does not independently prove that every oligo lot used in the
  2021 experiment was chemically identical. The paper does not identify each sample's
  index-primer number; index 1 on the public page is a kit example.
- 🟢 Sequencing was single-end 50 cycles generally and 75 cycles for sperm-head samples.

## Boundary

The defining operation is the combination of terminal repair and internal demethylation.
Neither the untreated small-RNA library nor either single-enzyme control is PANDORA-seq.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC8236090 JATS309–331, AlkB/PNK treatment, ligation capability assay and library preparation.

Specified T4PNK37°C20min and AlkB37°C30min, plus 11–22 library PCR cycles. Kept the Takara ligation-capability test oligo separate from NEBNextE7330S library adapters; confirmed SE50 versus sperm-head SE75.

Deferred: Original preceding RNA-size-selection/order and Fig1 workflow, manual-version adapter assignments, interposed RT-primer annealing, cleanup states and sample-specific run choice. This pass does not establish a complete source review or publication-wide absence.

Vendor manual now resides in this protocol source bundle; adapter-source guard passes. Public caveat distinguishes later manual from unverified historical kit revision.

## Native graph migration (2026-10-18, batch03)

Inspected: Published PMC8236090 JATS Introduction/Results154–168 (modifications and preselection); Methods RNA isolation/size selection around297–307, AlkB/PNK and separate adapter-validation assay309–319, and Library preparation329–331. Whole-bundle scraper completed: one23MB workbook skipped by ordinary size cap, no negative conclusion drawn; assay/2025-manual hits were not promoted into historical library chemistry.

Native graph selects PNK-first/AlkB-second, with initial15–50nt selection, explicit TRIzol recovery after each pretreatment and source library-stage15–45nt reselection. Source RT-primer annealing occurs between3′ and5′adapter additions.2025 E7330 ordering oligos remain clearly labeled comparator inventory only; no2025 adapter bases/modifications, example index member or fixed6nt index width enters actual2021 molecular states. Takara validation-assay chemistry remains separate. Historical adapter-ligated RNA/primer populations have separate strands; actual ss first-strand target feeds11–22cycle PCR, with RNA-template fate molecular/unresolved rather than invented removal. Source postPCRPAGE and SE50/general versus SE75/sperm alternatives are explicit.

Deferred: Historical2021 E7330S revision-specific oligos/end modifications and native binding/copying/ligation geometry; selected sample index members/widths/acquisition; exact finalPAGE window and sample-specific run assignment remain unassigned. Later vendor comparator was not newly read or accepted as revision-equivalent. This capped pass does not establish complete molecular proof or publication-wide absence.
