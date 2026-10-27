# RecE/RecT ET sources

Run `python3 tools/get_sources.py --slugs et-recombineering__10.1038+2417`. Source identities and retrieval results live in `/husky/henriksson/NGS101_papers/sources/et-recombineering__10.1038+2417/MANIFEST.tsv`.

The defining published abstract was inspected at PMID9771703. The full article DOI10.1038/2417 remains manual retrieval; no source bases are adopted. The bundle scraper ran and found no files, which is not evidence that sequences are absent. No third-party file is committed. See root `geturl.txt` for outstanding requests.

A direct-original follow-up identified an [author-uploaded published copy](https://www.researchgate.net/publication/13513898_Zhang_Y_Buchholz_F_Muyrers_JP_Stewart_AF_A_new_logic_for_DNA_engineering_using_recombination_in_E_coli_Nat_Genet_20_123-128), attributed to Adrian Stewart. Its printed title/journal/page fragment agrees with the defining identity. The [PDF](https://www.researchgate.net/profile/Adrian-Stewart-2/publication/13513898_Zhang_Y_Buchholz_F_Muyrers_JP_Stewart_AF_A_new_logic_for_DNA_engineering_using_recombination_in_E_coli_Nat_Genet_20_123-128/links/0deec517ba6a802d3a000000/Zhang-Y-Buchholz-F-Muyrers-JP-Stewart-AF-A-new-logic-for-DNA-engineering-using-recombination-in-E-coli-Nat-Genet-20-123-128.pdf) is registered in `catalogue/extra_sources.tsv`, but retrieval returned HTTP403; its methods and figures remain uninspected. The fragmented web preview is not a substitute for that review.

27 October 2026 resumed-goal availability check: the registered direct publisher PDF `https://www.nature.com/articles/ng1098_123.pdf` returned HTML and was rejected by PDF validation. The author-upload PDF again returned HTTP403. No original methods or figures were inspected; ET remains private.
