# Source retrieval

Run `python3 tools/get_sources.py --slugs lace-seq__10.1038+s41556-021-00696-9`. Originals and text conversions stay in the external source store, with URL identity in MANIFEST.tsv.

- https://www.nature.com/articles/s41556-021-00696-9 — publisher HTML retrieved, but full Methods not accessible in its body.
- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41556-021-00696-9/MediaObjects/41556_2021_696_MOESM2_ESM.xlsx — published workbook retrieved; TableS1 and selected TableS10 inspected.
- https://doi.org/10.21203/rs.3.pex-1499/v1 — exact cited detailed protocol revision1; standard source route records Cloudflare barrier, original unread.

No third-party documents belong in this checkout. No recipe inferred from comparison methods or a newer protocol revision.
