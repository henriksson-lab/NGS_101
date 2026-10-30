# Reference retrieval

Run `python3 tools/get_sources.py --slugs molseq__10.1038+s41598-026-41186-y` and mandatory `tools/scrape_primers.py` against `$CHEM_DATA/sources/molseq__10.1038+s41598-026-41186-y/`. The default external root is `/husky/henriksson/NGS101_papers/sources/`; no originals or conversions belong in the checkout.

- Published main: https://pmc-oa-opendata.s3.amazonaws.com/PMC13136315.1/PMC13136315.1.xml — external `s41598-026-41186-y_PMC13136315.xml(.txt)` and same-original PDF. Main54–70/100–110/124–168; actual MethodsPDFp9 confirms132nM stock wording.
- Published SI: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41598-026-41186-y/MediaObjects/41598_2026_41186_MOESM1_ESM.pdf — external `s41598-026-41186-y_supp_41598_2026_41186_MOESM1_ESM.pdf(.txt)`. Figs1/2/4, actual Table1p6–8 with whole modified orderings, photocleavable-linker semantics and final enrichment primers.

Two originals read. Automatically discovered superseded preprint sources are excluded; rejected responses are not evidence. There is no failed needed original in this selected pass. Further linked-source and full molecular work remains open without an absence claim.
