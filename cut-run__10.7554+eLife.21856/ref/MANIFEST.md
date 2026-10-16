# CUT&RUN originals

Retrieve registered originals with `python3 tools/get_sources.py --slugs cut-run__10.7554+eLife.21856`. Downloads and manifest remain under `$CHEM_DATA/sources/cut-run__10.7554+eLife.21856/`, never committed here.

The seven GSE84474 protocol downloads are original `.rtf.gz` files. Decompress each with Python standard-library `gzip` into the adjacent `.rtf`, then run `python3 tools/doctext.py FILE.rtf`. RTF conversion requires optional LibreOffice. These decompression/conversion steps were run separately for this audit; fetching the gzip alone does not imply its body was converted or read.

Run `python3 tools/scrape_primers.py SOURCE_DIRECTORY` before reading sources. Read the original-source text twins and source-specific linked material; the checklist in `catalogue/source_reviews/` records located readings and remaining obligations. The KAPA article URL currently serves an unrelated shop page: its rejected response is outside the source bundle and recorded as a manual identity failure, not a valid vendor original.
