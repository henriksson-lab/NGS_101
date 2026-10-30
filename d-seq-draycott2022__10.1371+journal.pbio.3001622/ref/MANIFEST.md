# Reference retrieval

Run `python3 tools/get_sources.py --slugs d-seq-draycott2022__10.1371+journal.pbio.3001622` and mandatory `tools/scrape_primers.py` against `$CHEM_DATA/sources/d-seq-draycott2022__10.1371+journal.pbio.3001622/` (default `/husky/henriksson/NGS101_papers/sources/`). Originals/text never go in the checkout.

- Publishedmain XML: https://pmc-oa-opendata.s3.amazonaws.com/PMC9129914.1/PMC9129914.1.xml; PDF https://pmc-oa-opendata.s3.amazonaws.com/PMC9129914.1/PMC9129914.1.pdf — Methods155/159–164,Results108,actualFig1PDFp3.
- OriginalPLOS main: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001622 — same published article, not another chemistry source.
- Actualsupplement URLs/identities in externalMANIFEST.tsv; none opened as an original in this capped pass. LargeEPS/raw-image files were downloaded but not mistaken for inspected sources.

Oneoriginalpositivepass; no failed needed original identified. Remaining readable-but-uninspected sources are ordinary deferred work, not absence evidence.
