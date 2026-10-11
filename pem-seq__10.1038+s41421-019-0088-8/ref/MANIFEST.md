# PEM-seq source manifest

No third-party source files are committed.

Retrieve the defining article, Supplementary Table S11, and the authors' detailed protocol with:

```sh
python3 tools/get_sources.py "PEM-seq"
```

The source registry includes the two PMC archives and the protocol's Table 2 image; `MANIFEST.tsv` records their origins and checksums under `$CHEM_DATA`.
