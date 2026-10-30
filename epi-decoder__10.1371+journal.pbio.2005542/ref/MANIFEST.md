# Re-obtaining inspected originals

All originals remain outside checkout under `/husky/henriksson/NGS101_papers/sources/epi-decoder__10.1371+journal.pbio.2005542/`. Run `OMP_THREAD_LIMIT=1 python3 tools/get_sources.py --slugs epi-decoder__10.1371+journal.pbio.2005542`.

- https://pmc-oa-opendata.s3.amazonaws.com/PMC6059479.1/PMC6059479.1.xml
- https://pmc-oa-opendata.s3.amazonaws.com/PMC6059479.1/pbio.2005542.s012.xlsx
- https://journals.plos.org/plosbiology/article/file?id=10.1371%2Fjournal.pbio.2005542.s001&type=supplementary

Main XML/PDF are representations of one published original. S7 XLSX has exact source orders. S1 is a validated TIFF saved as file.html, not an HTML chemistry source. Parent extra_sources declares its actual URL. No source files are checked in.
