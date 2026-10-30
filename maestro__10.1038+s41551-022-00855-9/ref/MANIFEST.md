# Source retrieval

Run `python3 tools/get_sources.py --slugs maestro__10.1038+s41551-022-00855-9`. Originals and extracted text live externally under `/husky/henriksson/NGS101_papers/sources/maestro__10.1038+s41551-022-00855-9/`, never in this checkout. Read the source-store MANIFEST.tsv for URL/version identity.

- https://pmc-oa-opendata.s3.amazonaws.com/PMC9089460.1/PMC9089460.1.xml — Author manuscript Methods163–205; title/version identity verified; full final main not compared
- https://pmc-oa-opendata.s3.amazonaws.com/PMC9089460.1/NIHMS1775951-supplement-Probe_Sequences_File.xlsx — First438SNVpanel header/first rows and exact first probe/coordinates
- https://pmc-oa-opendata.s3.amazonaws.com/PMC9089460.1/NIHMS1775951-supplement-Supplementary_Figures.pdf — Selected probe-extension explanation around extractedline297; remaining figures not audited

Other fetched files are not automatically evidence read. This is a capped source pass, not a two-reader absence audit.
