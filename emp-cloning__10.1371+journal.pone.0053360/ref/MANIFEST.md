# EMP cloning source originals

Originals live in `$CHEM_DATA/sources/emp-cloning__10.1371+journal.pone.0053360/`, default root
`/husky/henriksson/NGS101_papers/`. The external `MANIFEST.tsv` retains original URLs,
retrieval outcomes and file identities. No downloaded originals are committed.

- [Defining published article](https://doi.org/10.1371/journal.pone.0053360); [PMC3534072](https://pmc.ncbi.nlm.nih.gov/articles/PMC3534072/).
- Main inspected cached original: `journal.pone.0053360_PMC3534072.pdf` and its `.txt` twin.
- Article-linked figures and supplements keep their original identity in the manifest.

```sh
python3 tools/get_sources.py --slugs emp-cloning__10.1371+journal.pone.0053360
python3 tools/scrape_primers.py "$CHEM_DATA/sources/emp-cloning__10.1371+journal.pone.0053360/"
```

The initial `--doi` bundle was moved intact from its temporary canonical DOI directory
to this registered protocol source directory; no duplicate compatibility copy was made.
Retrieval does not imply complete source review. See the research note and source-review JSON.
