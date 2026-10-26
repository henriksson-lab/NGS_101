# Source retrieval

Originals and extracted text live under `$CHEM_DATA/sources/in-fusion__10.1093+nar+gkm047/`, default `/husky/henriksson/NGS101_papers/sources/in-fusion__10.1093+nar+gkm047/`. No third-party material is checked in.

Re-obtain with `python3 tools/get_sources.py --slugs in-fusion__10.1093+nar+gkm047`. The defining published DOI is 10.1093/nar/gkm047; the fetched manifest identifies article PDF/JATS `nar_gkm047_PMC1874605.xml` and every additional file. Run `python3 tools/scrape_primers.py $CHEM_DATA/sources/in-fusion__10.1093+nar+gkm047/` when reading the bundle. Availability and inspected sections are recorded in the protocol source-review ledger; downloaded does not mean inspected.
