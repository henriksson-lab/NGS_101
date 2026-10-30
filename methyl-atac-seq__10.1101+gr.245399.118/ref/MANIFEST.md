# Source manifest

Originals live outside the checkout at `/husky/henriksson/NGS101_papers/sources/methyl-atac-seq__10.1101+gr.245399.118`. Re-obtain with `python3 tools/get_sources.py --slugs methyl-atac-seq__10.1101+gr.245399.118`; then mandatory `tools/scrape_primers.py`. No third-party originals or extracted text are checked in.

Selected published originals/sections:

- [gr.245399.118_PMC6581052.xml](https://pmc-oa-opendata.s3.amazonaws.com/PMC6581052.1/PMC6581052.1.xml): PublishedVORPMC6581052 mainMethods102–116, analysis140; selected preparation source and sequencing; genotype/chip/RNA comparator recipes distinguished.
- [gr.245399.118_supp_supp_gr.245399.118_Supplemental_Table_S1-S3.pdf](https://pmc-oa-opendata.s3.amazonaws.com/PMC6581052.1/supp_gr.245399.118_Supplemental_Table_S1-S3.pdf): PublishedTablesS1–S3; S3 oligos extractedlines96–102 and renderedoriginalp1; S1 sample inventory/S2 validationprimer inventory, not active library inputs.
- [gr.245399.118_supp_supp_gr.245399.118_Supplemental_Fig_S1-S6.pdf](https://pmc-oa-opendata.s3.amazonaws.com/PMC6581052.1/supp_gr.245399.118_Supplemental_Fig_S1-S6.pdf): PublishedFigS1–S6 captions, especiallyS5 titration andS6 before/afterconversionQC; no exhaustive numericalplot review.

Deferred directly named sources:

- MainMethods102 Picelli2014 Tn5production andAdey/Shendure2012 assembly: Namedpredecessor originals not followed within cap; only directly printed assembly parameters selected.
- MainMethods108/112 Corces2017OmniATAC andLu2015T-WGBS/Suzuki2018X-WGBS: Actual citedcomponent predecessors; only local sourceprinted lysis/inactivation/sizeselection/endrepair used.
- Main112ZymoD5030T/PfuTurboCx600410;108AMPureA63880 historical guides: Exacthistoricalkit/end/acquisitiongeometry unreviewed, not absence.

The source-review JSON retains actual failed-manifest URLs and additional uninspected inventory. A download or scraper match alone does not establish a chemical fact. No accepted missing-information claims.
