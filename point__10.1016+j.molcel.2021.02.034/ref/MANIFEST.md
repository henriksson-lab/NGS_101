# Source manifest

Retrieve with `python3 tools/get_sources.py --slugs point__10.1016+j.molcel.2021.02.034`. Originals/text conversions belong in `$CHEM_DATA/sources/` (default `/husky/henriksson/NGS101_papers/sources/`), never in the checkout.

- article: https://pmc-oa-opendata.s3.amazonaws.com/PMC8122139.1/PMC8122139.1.xml — read.
- chemistry_supplement: https://pmc-oa-opendata.s3.amazonaws.com/PMC8122139.1/mmc1.pdf — read.
- available_original_5: https://pmc-oa-opendata.s3.amazonaws.com/PMC8122139.1/mmc2.pdf — available.

Source identity: VOR; CC BY, author-manuscript false. Exact inspected sections and named deferrals are in `01_point.md` and corresponding `catalogue/source_reviews/` JSON. This capped pass does not establish a full review or missing disclosure.

The auto-listed2020preprint is superseded; its failed download is excluded from needed current sources.
