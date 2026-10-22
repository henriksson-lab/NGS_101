# in4mer reference originals

Defining source: [Esmaeili Anvar et al. 2024](https://doi.org/10.1038/s41467-024-47795-3),
[PMC11055879](https://pmc.ncbi.nlm.nih.gov/articles/PMC11055879/).

The originals and their plain-text twins live under
`$CHEM_DATA/sources/in4mer__10.1038+s41467-024-47795-3/`.
The external `MANIFEST.tsv` retains URL identities, hashes and duplicate records.
The publisher ESM route and the explicit `pmccloud:PMC11055879` extra-source row
retrieve the article and published supplementary files. No third-party files belong
in this checkout.

Supplementary Data 10 is publisher file
`41467_2024_47795_MOESM13_ESM.xlsx`, sheet `Supplementary Data 10`.
Large result workbooks are not substitutes for this oligo workbook.

After catalogue registration:

```sh
python3 tools/get_sources.py --slugs in4mer__10.1038+s41467-024-47795-3
python3 tools/scrape_primers.py "$CHEM_DATA/sources/in4mer__10.1038+s41467-024-47795-3/"
```

See `01_in4mer.md` and the pooled-cloning round ledger for inspected sections,
branch distinctions and outstanding original-vector checks.
