# AlkAniline-Seq

**Defining source:** Marchand *et al.*, *Angewandte Chemie International Edition*
(2018), doi:10.1002/anie.201810946. The article, graphical abstract, supplementary-method
description and later detailed protocol were read. The author-hosted
[paper PDF](https://www.lafontainelab.com/wp-content/uploads/Marchand2018.pdf)
and the [NEB small-RNA kit manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale7300_e7330_e7560_e7580.pdf)
(version 9.0, July 2025) were also checked. Wiley's 2018 Supporting Information
currently returns HTTP 403 to the source fetcher.

## Evidence used

- 🟢 Mild alkaline hydrolysis fragments RNA and forms abasic sites at m7G, m3C and other
  susceptible modifications.
- 🟢 Antarctic phosphatase removes all existing 5′ and 3′ phosphates before aniline.
- 🟢 Treatment with 1 M aniline, pH 4.5, for 15 min at 60 °C cleaves abasic sites and
  creates a new 5′ phosphate at N+1.
- 🟢 NEBNext Small RNA 5′-adapter ligation positively selects those new phosphates;
  RT, indexed PCR and SR50 Illumina sequencing follow.
- 🟢 The later author protocol, Marchand *et al.*, *Methods in Enzymology* (2021),
  doi:10.1016/bs.mie.2021.06.001, specifies the NEBNext Multiplex Small RNA kit
  family (Sets 1/2 or E7330). NEB's manual, p. 20, publishes the 3′ and 5′
  adapters, the SR RT primer, the SR PCR primer and indexed PCR primers, including
  rApp, terminal amine, RNA bases and phosphorothioate bonds. The public page now
  carries those ordered oligos; its index 1 primer is a kit example, because the
  original experiment's index selections were not specified.
- 🟢 [NEB's sequencing FAQ](https://www.neb.com/en/faqs/do-i-have-to-use-custom-sequencing-primers-for-my-libraries-generated-with-nebnext-small-rna-library-prep-kit-for-illumina)
  says this kit family uses standard Illumina sequencing reagents without a custom
  sequencing primer.
- 🟡 The manual's 2025 version identifies the current kit formulation, while the
  2018 Supporting Information remains unavailable. Exact historical lot equivalence
  is not established by the accessible sources.

## Construct logic

The phosphatase step is the selectivity reset: after it, the intended ligatable
5′ phosphate arises from aniline cleavage. The original article places the
first mapped Read 1 base at N+1. The local model declares that source relation;
its chemistry checker does not independently prove cleavage selectivity or yield.

## Bounded primary-source check (2026-10-16)

Inspected: Original Marchand2018 author PDF text90–145 and figure1 description658–661.

No source-supported correction identified in inspected scope.

Deferred: Original supporting methods; exact protocol timings, adapter assignment, run and full library states were not re-verified; main PDF confirms alkaline treatment, background dephosphorylation, aniline cleavage and N+1 selective-ligation principle. This is a capped check, not a full source audit or confirmation of missing information.

## Native ownership migration — 2026-10-18

Re-read the original Marchand2018 article, text lines99–109 and 141–188, plus
rendered Figure 1a on PDF page2. It establishes alkaline hydrolysis in bicarbonate
pH9.2 at 96 °C, extensive background 5′/3′ dephosphorylation, aniline cleavage,
adapter ligation and dsDNA amplification. Read 1 begins at N+1. The graph now
declares those steps locally, including ligation-competent downstream 5′-P/3′-OH
RNA and distinct upstream/uncleaved outcomes with unresolved physical recovery.
Both adapter boundaries are marked on the ligated RNA candidate.

The exact adapters and primers in NEB manual v9.0 (July2025), including rApp,
RNA bases, terminal amine and phosphorothioates, remain a labelled comparator
inventory. That later formulation does not establish the historical 2018 kit
implementation and supplies no native graph reagent or run assignment. The
2018 Supporting Information remains an unresolved source target. Native RT/PCR
reagents and primer footprints therefore remain unassigned on the actual cDNA
precursor; the output duplex is not substituted as a PCR input.

The earlier notes report SR50, but the original main passages re-read here
establish Read 1 start semantics without its cycle count. That cycle assignment
is now a deferred run question rather than an inherited acquisition default.
Similarly, no sample-index width or member is imposed from the 2025 comparator.
The shared recipe and shared readout constructor no longer choose these source
decisions. No absence conclusion or historical-kit equivalence is established.

Two originals inspected: the 2018 article (including Figure 1) and the 2025 NEB
manual's printed version and oligo table, text961–1073. The whole cached-bundle
primer scraper was run and relevant modified-oligo hits inspected.
