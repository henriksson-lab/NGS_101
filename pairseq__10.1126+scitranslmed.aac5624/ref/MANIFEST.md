# Reference retrieval

Run `python3 tools/get_sources.py --slugs pairseq__10.1126+scitranslmed.aac5624`, followed by `python3 tools/scrape_primers.py $CHEM_DATA/sources/pairseq__10.1126+scitranslmed.aac5624/`. The default source root is `/husky/henriksson/NGS101_papers/sources/`; no originals belong in the checkout.

- Defining article: https://doi.org/10.1126/scitranslmed.aac5624 — standard routes found no open full text.
- Published erratum: https://www.science.org/doi/full/10.1126/scitranslmed.aad5480 — HTTP403.
- Indexed corrected-main mirror: https://www.therenegadecell.com/wiki/images/c/cb/High-throughput_pairing_of_TCR_alpha_and_beta_sequences.pdf — HTTP404; no source was received or read. This failed mirror is not an additional requested original.

The article and erratum are needed; chemistry supplement identities must come from the actual article. No guessed supplement URL is requested. The mandatory scraper found no source content; no absence claim follows.
