# SLIM source originals

Originals live in `$CHEM_DATA/sources/slim__10.1093+nar+gnh172/`, default root
`/husky/henriksson/NGS101_papers/`. The external `MANIFEST.tsv` retains original URLs,
retrieval outcomes and file identities. No downloaded originals are committed.

- [Defining published article](https://doi.org/10.1093/nar/gnh172); [PMC535700](https://pmc.ncbi.nlm.nih.gov/articles/PMC535700/).
- Main inspected cached original: `nar_gnh172_PMC535700.html` and its `.txt` twin.
- Article-linked figures and supplements keep their original identity in the manifest.

```sh
python3 tools/get_sources.py --slugs slim__10.1093+nar+gnh172
python3 tools/scrape_primers.py "$CHEM_DATA/sources/slim__10.1093+nar+gnh172/"
```

The initial `--doi` bundle was moved intact from its temporary canonical DOI directory
to this registered protocol source directory; no duplicate compatibility copy was made.
Retrieval does not imply complete source review. See the research note and source-review JSON.
