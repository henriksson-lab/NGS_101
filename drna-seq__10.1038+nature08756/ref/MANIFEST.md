# Reference retrieval

Run `python3 tools/get_sources.py --slugs drna-seq__10.1038+nature08756` and `python3 tools/scrape_primers.py $CHEM_DATA/sources/drna-seq__10.1038+nature08756/`. Originals/text remain outside the checkout under `$CHEM_DATA/sources/drna-seq__10.1038+nature08756/` (default `/husky/henriksson/NGS101_papers/sources/`).

- Published main: https://www.nature.com/articles/nature08756 — saved HTML is subscription preview; full article still needed.
- Published Supplementary Methods: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature08756/MediaObjects/41586_2010_BFnature08756_MOESM206_ESM.pdf — Methods237–282 plus oligo-role/companion-library boundary read.
- Published Supplementary Tables: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature08756/MediaObjects/41586_2010_BFnature08756_MOESM208_ESM.pdf — actualTable2p6/MLbarcodes read; remaining tables not inspected as source originals.
- Other exact supplement identities/URLs are in externalMANIFEST.tsv and remain unreviewed.

Three-original capped positive pass, not full linked-source or absence audit. No third-party original/text is checked in.
