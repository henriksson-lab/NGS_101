# Source manifest

Retrieve with `python3 tools/get_sources.py --slugs cutac__10.7554+eLife.63274`, then run the oligo scraper on the full bundle before reading.

- [Defining eLife63274 original](https://doi.org/10.7554/eLife.63274), PMC7721439; original dataset workbook and transparent-reporting DOCX.
- Registered [8yrhxv6 in-house pA–Tn5 loading protocol](https://doi.org/10.17504/protocols.io.8yrhxv6), November26,2019, integer29425. The standard fetcher recovers its original PDF and anonymous public API table.
- Registered [bmbfk2jn CUTAC Version2](https://doi.org/10.17504/protocols.io.bmbfk2jn), September18,2020, integer42055; original PDF/API.
- Cited custom index table: [Buenrostro2015](https://doi.org/10.1038/nature14590), MOESM36 XLSX, retrieved with `python3 tools/get_sources.py --doi 10.1038/nature14590`.

Exact independent read scopes and pending linked sources are recorded in `catalogue/source_reviews/`. No third-party originals or extracted twins are committed.
