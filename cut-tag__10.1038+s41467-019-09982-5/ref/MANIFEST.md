# Source manifest

Re-obtain with `python3 tools/get_sources.py --slugs cut-tag__10.1038+s41467-019-09982-5` and run the oligo scraper on the resulting bundle before reading.

- Defining [2019 paper](https://doi.org/10.1038/s41467-019-09982-5), PMC6488672, and original supplements.
- Cited adapter-construction original [Picelli2014](https://doi.org/10.1101/gr.177881.114), PMC4248319: registered PMC Cloud route retrieves main XML/PDF and supplemental information.
- Historical author protocols.io v1 wnufdew (integer 18868) and v2 z6hf9b6 (22441): registered version-specific `/v1` and `/v2` exports and public API preserve identity. The default export is later v3 bcuhiwt6 (33385).
- Cited index-source [Buenrostro2015](https://doi.org/10.1038/nature14590), original MOESM36 workbook: retrieve with `python3 tools/get_sources.py --doi 10.1038/nature14590`.

The machine-readable review in `catalogue/source_reviews/` records original locations, independent reading scopes and remaining linked-source work. No third-party originals, images or extracted twins are committed.
