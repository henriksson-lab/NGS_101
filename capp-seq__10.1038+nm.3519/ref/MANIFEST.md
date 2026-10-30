# Source retrieval

Run `python3 tools/get_sources.py --slugs capp-seq__10.1038+nm.3519`. Originals and extracted text live externally under `/husky/henriksson/NGS101_papers/sources/capp-seq__10.1038+nm.3519/`, never in this checkout. Read the source-store MANIFEST.tsv for URL/version identity.

- https://pmc-oa-opendata.s3.amazonaws.com/PMC4016134.1/PMC4016134.1.xml — Author manuscript Methods128–146; selected plasmaDNA library/capture branch
- https://pmc-oa-opendata.s3.amazonaws.com/PMC4016134.1/NIHMS538716-supplement-1.pdf — Associated SI A3/A4text1683–1716 and diagnostic primer source1731–1745; historicalqualified
- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnm.3519/MediaObjects/41591_2014_BFnm3519_MOESM16_ESM.pdf — Published SI A3text1707–1735 controls inspected optimization; other SIsections unreviewed

Other fetched files are not automatically evidence read. This is a capped source pass, not a two-reader absence audit.
