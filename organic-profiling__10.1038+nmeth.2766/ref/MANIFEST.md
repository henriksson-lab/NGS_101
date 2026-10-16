# Source manifest

Recover registered originals with:

```sh
python3 tools/get_sources.py --slugs organic-profiling__10.1038+nmeth.2766
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/organic-profiling__10.1038+nmeth.2766/
```

- Defining source: doi:10.1038/nmeth.2766, PMC3929178. The author manuscript's extracted body omits Online Methods; published MOESM254 p. 3 is the original detailed workflow figure.
- Original study metadata: GSE45672, registered GEO family metadata URL. Its series line 51 links the small author protocol `GSE45672_ORGANIC_profiling_protocol.rtf.gz`; do not fetch the adjacent RAW archive to read this protocol.
- The RTF's yeast native-protocol precursor is doi:10.1371/journal.pgen.1003317, PMC3585014. Main original and two small DOC supplements are recovered through the registered DOI route; remaining TIFF supplement availability is recorded in the source manifest.
- Historical library precursor: doi:10.1073/pnas.1110731108, PMC3215028. The `pmc:` route retrieves the original main HTML; the Duke course mirror is an additional original-PDF copy. The actual PNAS `STXT` SI anchor and PMC instance PDF routes are registered; failed retrieval is not scientific absence.
- Later exact library illustration: doi:10.1002/0471142727.mb2131s110, PMC4410783, Orsi2015 Basic Protocol3. Its oligos and recipe are a separately attributed later implementation.
- Original author methods entry and publisher preview are registered for source/version and link inventories. The current author page links PubMed articles; it does not establish historical recipe identity.

Third-party originals, extracted twins and inspection images belong under ignored `$CHEM_DATA/sources/` and `_data/tmp/`. Exact pass locations, withdrawn historical claims and outstanding originals are in `catalogue/source_reviews/organic-profiling__10.1038+nmeth.2766.json`.
