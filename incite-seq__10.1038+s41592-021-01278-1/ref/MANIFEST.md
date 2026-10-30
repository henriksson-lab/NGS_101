# Reference retrieval

Originals and text twins live outside the checkout in `$CHEM_DATA/sources/incite-seq__10.1038+s41592-021-01278-1/`. Re-run `python3 tools/get_sources.py --slugs incite-seq__10.1038+s41592-021-01278-1` and mandatory `tools/scrape_primers.py` on that bundle.

- s41592-021-01278-1_PMC8532076.xml — https://pmc-oa-opendata.s3.amazonaws.com/PMC8532076.1/PMC8532076.1.xml — AssociatedNIHMS1745895identity; selectedMethods173–209, TotalSeq-Aformat197, directnuclearstaining/capture/physicallibraryfractions.
- s41592-021-01278-1_41592_2021_1278_MOESM1_ESM.pdf — https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41592-021-01278-1/MediaObjects/41592_2021_1278_MOESM1_ESM.pdf — Published reporting-summary identity/study/antibody/flowcytometryscope. Verified reportingsummary, NOTan oligotable or central chemistryrecipe.
- s41592-021-01278-1.html — https://www.nature.com/articles/s41592-021-01278-1 — Published finalNatureidentity/abstract/EDcaptions read; actualfullMethods notexposed: subscriptionpreview.

Actual needed sources:

- https://www.nature.com/articles/s41592-021-01278-1 — Actual finalpublishedmain contains subscriptionpreview only; fullMethods neededforfinal/AMversionequivalence and remainingprepcontext.

Other automatic downloads, including superseded preprints, were not adopted as recipe authority. No originals belong in this checkout.
