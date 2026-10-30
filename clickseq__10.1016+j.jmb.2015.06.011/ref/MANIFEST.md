# Source manifest

Retrieve with `python3 tools/get_sources.py --slugs clickseq__10.1016+j.jmb.2015.06.011`.
All originals, text twins and URL/checksum/failure records remain outside the checkout.

- Defining publication: https://doi.org/10.1016/j.jmb.2015.06.011
- Associated accepted manuscript: https://pmc.ncbi.nlm.nih.gov/articles/PMC4523409/
- Main JATS original obtained from https://pmc-oa-opendata.s3.amazonaws.com/PMC4523409.1/PMC4523409.1.xml
- Actual named supplementary original, declared in catalogue/extra_sources.tsv: https://pmc.ncbi.nlm.nih.gov/articles/instance/4523409/bin/NIHMS705616-supplement.docx

The main is an author manuscript under TDM, not the final published methods. The exact DOCX request reached a preparing-download gate; the main's empty scraper result does not establish absent oligo disclosure. No downloaded material is checked in.
