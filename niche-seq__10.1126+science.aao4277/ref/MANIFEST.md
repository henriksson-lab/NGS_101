# Reference retrieval

Run `python3 tools/get_sources.py --slugs niche-seq__10.1126+science.aao4277` and mandatory `tools/scrape_primers.py` against `$CHEM_DATA/sources/niche-seq__10.1126+science.aao4277/`. Originals and extracted text stay outside the checkout under the configured source root.

- Associated author manuscript: https://pmc-oa-opendata.s3.amazonaws.com/PMC7234837.1/PMC7234837.1.xml — external `science.aao4277_PMC7234837.xml(.txt)`, main47/49/51 and references7/17. Header identifies EMS86422/Funders Author Manuscript; final-body equivalence unverified.
- Actual detailed supplementary material: https://pmc.ncbi.nlm.nih.gov/articles/instance/7234837/bin/EMS86422-supplement-Supplementary_Materials_incl_Figures.pdf — PMC host gate/cooldown; not supplied as a Cloud TDM object. Reference17 points to supplementary Methods.
- Actual study tables: https://pmc.ncbi.nlm.nih.gov/articles/instance/7234837/bin/EMS86422-supplement-Supplementary_Tables.xlsx — gated, unread; chemistry role unreviewed, not assumed an oligo disclosure.

Only one original was read. Detailed Methods and whole library/read states remain open; source acquisition failures do not establish missing information.
