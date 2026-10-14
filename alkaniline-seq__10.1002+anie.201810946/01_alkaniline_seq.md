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

The phosphatase step is the selectivity reset: after it, a ligatable 5′ phosphate should
only arise from aniline cleavage. Consequently the first mapped read base is one nucleotide
downstream of the responsible modification, a relationship enforced by the shared model.
