# Source retrieval

Originals and extracted text live under `$CHEM_DATA/sources/goldenbraid__10.1371+journal.pone.0021622/`, default `/husky/henriksson/NGS101_papers/sources/goldenbraid__10.1371+journal.pone.0021622/`. No third-party material is checked in.

Re-obtain with `python3 tools/get_sources.py --slugs goldenbraid__10.1371+journal.pone.0021622`. The defining published DOI is 10.1371/journal.pone.0021622; the fetched manifest identifies article PDF/JATS `journal.pone.0021622_PMC3131274.xml` and every additional file. Run `python3 tools/scrape_primers.py $CHEM_DATA/sources/goldenbraid__10.1371+journal.pone.0021622/` when reading the bundle. Availability and inspected sections are recorded in the protocol source-review ledger; downloaded does not mean inspected.
