# Sources

No third-party source files are committed. Re-obtain originals with:

```sh
python3 tools/get_sources.py --slugs nice-seq__10.1186+s13059-017-1247-6
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/nice-seq__10.1186+s13059-017-1247-6
```

Registered originals include the defining PMC5488340 article and S1–S17 DOCX,
NEB E7370 manual (fetched revision 9.1_1/23, not a historical study oligo assignment),
and original GSE97889 family metadata. The registered current E7370 product page returned
HTTP 403 during this audit. Sources, extracted text and embedded figure derivatives stay
under `$CHEM_DATA`, never in this directory. `catalogue/source_reviews/` records exact
reading locations, independent passes and actual remaining source/link obligations.
