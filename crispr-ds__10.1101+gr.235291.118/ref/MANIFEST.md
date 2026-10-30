# Source retrieval

Run `python3 tools/get_sources.py --slugs crispr-ds__10.1101+gr.235291.118`. Originals/text twins live in the external source store `/husky/henriksson/NGS101_papers/sources/crispr-ds__10.1101+gr.235291.118/`, never in this checkout. The source-store MANIFEST.tsv identifies URLs, versions and rejected responses.

- https://pmc-oa-opendata.s3.amazonaws.com/PMC6169890.1/PMC6169890.1.xml — Published main Methods140–172,selectedFig1C originalPDFp3 /tmp/crisprds22-fig1-actual.png;referenced2014identitywithinmain
- https://pmc-oa-opendata.s3.amazonaws.com/PMC6169890.1/supp_gr.235291.118_Supplemental_Tables.pdf — SupplementalTablesS1crRNA+PAM/tableS2captureprobe inventories;resultstableunreviewed
- https://pmc-oa-opendata.s3.amazonaws.com/PMC6169890.1/supp_gr.235291.118_Supplemental_Materials.docx — Administrativecontentsonly;thirdoriginalcountsincap,notunderlyingfigure/data/codechemistryreview

Other fetched files remain uninspected unless explicitly listed. No challenges were bypassed. Unsuccessful source-index aliases do not establish missing original chemistry.
