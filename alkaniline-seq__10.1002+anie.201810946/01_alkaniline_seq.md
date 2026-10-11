# AlkAniline-Seq

**Defining source:** Marchand *et al.*, *Angewandte Chemie International Edition*
(2018), doi:10.1002/anie.201810946. The article, graphical abstract, supplementary-method
description and later detailed protocol were read.

## Evidence used

- 🟢 Mild alkaline hydrolysis fragments RNA and forms abasic sites at m7G, m3C and other
  susceptible modifications.
- 🟢 Antarctic phosphatase removes all existing 5′ and 3′ phosphates before aniline.
- 🟢 Treatment with 1 M aniline, pH 4.5, for 15 min at 60 °C cleaves abasic sites and
  creates a new 5′ phosphate at N+1.
- 🟢 NEBNext Small RNA 5′-adapter ligation positively selects those new phosphates;
  RT, indexed PCR and SR50 Illumina sequencing follow.
- 🟡 The commercial kit's unprinted oligos remain molecular roles; its public 3′ adapter is shown.

## Construct logic

The phosphatase step is the selectivity reset: after it, a ligatable 5′ phosphate should
only arise from aniline cleavage. Consequently the first mapped read base is one nucleotide
downstream of the responsible modification, a relationship enforced by the shared model.
