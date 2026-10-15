# PEM-seq source manifest

No third-party source files are committed.

Retrieve the defining article, Supplementary Table S11, and the authors' detailed protocol with:

```sh
python3 tools/get_sources.py "PEM-seq"
```

`MANIFEST.tsv` under `$CHEM_DATA` records the origins and checksums of the two PMC archives
(PMC6434046, the defining paper and its supplementary PDF; PMC9019705, the STAR Protocols
procedure).

## Reading the STAR Protocols oligo table

Tables 1 and 2 of PMC9019705 are **vector objects, not text**. `pdftotext` and the JATS XML
return only their captions, and the registered `fx6.gif` is a 100×160-px thumbnail that
cannot be read. Render the pages instead:

```sh
pdftoppm -r 400 -f 10 -l 11 -png \
  "$CHEM_DATA/sources/pem-seq__10.1038+s41421-019-0088-8/PMC9019705_PMC9019705.pdf" table
```

Page 10 is Table 1 (bait DSBs; PAM in red), page 11 is Table 2 (DNA sequences). The colour
coding in Table 2 carries information that no text extraction preserves: the I5/I7 index
bases are red, and the bridge-adapter random molecular barcode is orange and underlined.
The 4-nt inline I5 index boundary in `tools/chem.py` was confirmed this way and
cross-checked inside Supplementary Table S11 — see `01_pem_seq.md`.
