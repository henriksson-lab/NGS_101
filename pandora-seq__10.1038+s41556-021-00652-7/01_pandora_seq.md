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
