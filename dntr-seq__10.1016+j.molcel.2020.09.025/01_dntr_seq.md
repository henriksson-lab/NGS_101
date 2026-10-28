# DNTR-seq — published-source reading blocked

Zachariadis et al. (2020), *A Highly Scalable Method for Joint Whole-Genome Sequencing and Gene-Expression Profiling of Single Cells*, [10.1016/j.molcel.2020.09.025](https://doi.org/10.1016/j.molcel.2020.09.025).

## Actual scope of this bounded pass

The standard source fetcher obtained no published article or chemistry supplement.
The published ScienceDirect page's indexed title/summary was available for identity
screening; opening its full article returned HTTP 403. The indexed summary supports
only the identity of a joint single-cell genome/transcriptome method called direct
nuclear tagmentation and RNA sequencing. It does not establish oligos, reaction
conditions or material fates. **No chemistry Methods were read in this pass.**

The source-bundle oligo scraper was run and returned no inspectable chemistry
source; this is an availability result, not evidence that sequences are absent.
The defining article is identified as Molecular Cell 80, 541–553.e5 (2020).
Its title/DOI identity is distinct from G&T-seq and DR-seq, but this pass does not
claim a source-verified molecular comparison to those protocols.

The fetcher also recorded the superseded preprint DOI 10.1101/2020.03.04.976530.
Its cooldown/manual rows are excluded from active chemistry and are not requests
for an alternative authority. No preprint recipe, current vendor kit or uncited
protocol revision fills the published-source gap.

## Deferred originals

- Published [main article and STAR Methods](https://www.sciencedirect.com/science/article/pii/S1097276520306559).
- Published chemistry supplements linked from that article; their identities,
  sections and individual retrieval URLs remain to inspect.
- Any directly cited detailed protocol and historical reagent/oligo sources,
  after the published Methods establishes their identity and applicability.

There is no molecular graph or preview: building one from the indexed abstract
would overstate the available evidence. Status remains **notes**, source review
**blocked**, with no accepted missing-information claim and no ancestry edge.
