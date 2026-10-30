# SISSOR source retrieval

Run `python3 tools/get_sources.py --slugs sissor__10.1073+pnas.1707609114`.
Sources are stored at `/husky/henriksson/NGS101_papers/sources/sissor__10.1073+pnas.1707609114/`.
The external `MANIFEST.tsv` records exact source URLs, checksums and original identity.

Selected published originals: PMC5703283 main PDF/JATS and
`pnas.1707609114.sapp.pdf` SI Appendix. Read the `.txt` twins and inspect original
mainPDFp6 for the printed3′phosphorylated-random-hexamer terminology. Do not
substitute inferred phosphorothioate protection or an ordinary3′OH primer.

Run `python3 tools/scrape_primers.py` on the external bundle. Published TableS9
primers are QC assays, not sequencing-library oligos. Automatically fetched
bioRxiv135384v1 and its supplement are superseded and excluded from chemistry
use; cached availability is not source applicability.

Directly cited Peters2012 (mainRef18/SIRef1, Nature487:190–195) and Adey2010
(mainRef22, GenomeBiology11:R119) are actual historical method follow-ups, not
assumed source-complete replacements. Exact vendor/primer revisions and full
molecular end/run geometry remain unreviewed. No third-party originals, extracted
text or rendered source figures are stored in this checkout.
