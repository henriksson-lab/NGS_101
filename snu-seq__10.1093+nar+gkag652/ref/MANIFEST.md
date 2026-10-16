# SNU-seq originals

Run `python3 tools/get_sources.py --slugs snu-seq__10.1093+nar+gkag652`, then run the full bundle oligo scraper before reading. Original source material, manifests and text twins remain under `$CHEM_DATA/sources/` and are not committed.

The defining 2026 NAR article and supplement are PMC13309788. The separately registered bioRxiv **v1** PDF prints July 14, 2021; the generic DOI route currently recovers **v2**, printed October 24, 2025. Their cache filenames differ and their chemistry/version claims must stay separate.

The exact dated manufacturer originals are Lexogen **012UG036V0130, June 13, 2018** (Ion Torrent) and **015UG009V0260, August 16, 2021** (Illumina FWD, with a separately numbered critical-recommendations insert). Their positive architectures do not identify an experimental sample's selected revision or barcode/index member.

Original GEO self metadata and the article-linked author repository landing are registered. The Dölken 2008 capture predecessor is also registered by DOI; the already-recovered institutional PDF and its text twin were independently read in the TT-seq source bundle, obtainable through `python3 tools/get_sources.py --slugs tt-seq__10.1126+science.aad9841`. See the SNU review JSON for exact source locators, read locations and still-uninspected source children. Raw reads and result archives are not needed for this bounded audit.
