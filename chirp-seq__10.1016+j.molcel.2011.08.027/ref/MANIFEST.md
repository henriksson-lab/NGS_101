# Source manifest

Originals live outside the checkout at `/husky/henriksson/NGS101_papers/sources/chirp-seq__10.1016+j.molcel.2011.08.027`. Re-obtain using `python3 tools/get_sources.py --slugs chirp-seq__10.1016+j.molcel.2011.08.027`, then run `tools/scrape_primers.py` on that directory.

Inspected originals:

- [j.molcel.2011.08.027_PMC3249421.xml](https://pmc-oa-opendata.s3.amazonaws.com/PMC3249421.1/PMC3249421.1.xml): Associated accepted manuscript, Methods86–124/overview; selected comparison only, subsequently obtained final PDF controls.
- [source_pdf_87b24b6fe07e40b97170.pdf](https://qulab.ustc.edu.cn/wp-content/uploads/2017/06/Genomic-Maps-of-Long-Noncoding-RNA-Occupancy-Reveal-Principles-of-RNA-Chromatin-Interactions.pdf): Published Molecular Cell44,667–678,18 November2011; rendered originalp1 title/DOI/date checked; Methods text1143–1275, probe design1148–1165, DNAelution1251–1263/library1271–1275.

Deferred/blocked targets:

- https://pmc.ncbi.nlm.nih.gov/articles/instance/3249421/bin/NIHMS340672-supplement-01.pdf: Actual supplement gated by existing PMC host cooldown; probe sequence TableS4 applicability uninspected.
- Johnson et al.2007 ChIP-seq library method cited in final Methods: Direct predecessor not followed within cap; specific adapter/index/primer/read geometry unassembled.

The source-review JSON records actual manifest URLs and availability. This is a capped selected-source reading, not full publication review or confirmed absence.
