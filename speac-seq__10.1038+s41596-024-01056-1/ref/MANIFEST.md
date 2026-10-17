# Reference manifest

- Defining source: https://doi.org/10.1038/s41596-024-01056-1 (Faust Akl et al.,
  *Nature Protocols*). Readable as the PMC author manuscript PMC11805652; the
  verbatim oligos are in its "Oligonucleotides" list.
- Supporting original study: https://doi.org/10.1126/science.abq4822 (PMC10066924).

Retrieve with

```sh
python3 tools/get_sources.py --slugs speac-seq__10.1038+s41596-024-01056-1   # Science original
python3 tools/get_sources.py --doi 10.1038/s41596-024-01056-1                # defining protocol
```

The second call writes to `$CHEM_DATA/sources/doi_10.1038+s41596-024-01056-1/`;
its files were copied into this protocol's source directory (the catalogue routes
this slug by the Science DOI). The PMC supplements of both papers are behind a PMC
download gate and remain `(manual)` rows.
