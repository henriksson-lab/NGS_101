# Source retrieval

Run `python3 tools/get_sources.py --slugs tam-seq__10.1126+scitranslmed.3003726`. Originals and extracted text live externally under `/husky/henriksson/NGS101_papers/sources/tam-seq__10.1126+scitranslmed.3003726/`, never in this checkout. Read the source-store MANIFEST.tsv for URL/version identity.

- https://s3-us-west-1.amazonaws.com/ptab-filings%2FIPR2025-01434%2F1010#final-article — Printed finaltitle/DOI/date; main extracted421–439 preparation paragraphs
- https://s3-us-west-1.amazonaws.com/ptab-filings%2FIPR2025-01434%2F1010#supplement — SIpp2–3;FigS1PDFp20 andTableS2PDFp26 original renders;TableS1handle definitions;firstFLD0001row

Other fetched files are not automatically evidence read. This is a capped source pass, not a two-reader absence audit.
