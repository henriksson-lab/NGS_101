# CAP-seq microbial long-read sources

Run `python3 tools/get_sources.py --slugs cap-seq-long-read-microbiome__10.1101+2024.09.10.612220`.
The external bundle is `/husky/henriksson/NGS101_papers/sources/cap-seq-long-read-microbiome__10.1101+2024.09.10.612220/`.
Read its `MANIFEST.tsv` for original URLs, checksums and duplicate-source entries.

Selected main original: bioRxiv10.1101/2024.09.10.612220v3, printed13October2025,
`2024.09.10.612220_v3.full.pdf` and its `.txt` twin. The same exact URL declared
in `extra_sources.tsv` is a duplicate reference, not another original. Verify the
printed revision before using a re-fetched copy. Do not merge ResearchSquare or
older preprint revisions; use a qualified published successor if one is identified.

The exact v3 `.supplementary-material` index is separately declared/fetched. Its
HTML title and links were inspected, but no actual SupplementaryTable1 original
was located. This is an unresolved source identification issue, not proof that the
source omits primers. The code repository is labelled bioinformatics code by the
main; no presumed oligo disclosure document is invented inside it.

Run `python3 tools/scrape_primers.py` on the external bundle. Its one main-text
sequence hit is a computational adapter-search motif, not a whole primer. Render
original mainp22 for mixture quantities and capsule/bead distinctions. Historical
material, bead and ONT references remain directly linked follow-ups with no recipe
inheritance yet. Originals/extracted text and rendered pages stay outside the repo.
