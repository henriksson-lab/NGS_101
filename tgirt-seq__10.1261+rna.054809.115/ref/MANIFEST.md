# Source manifest

Third-party originals remain under `$CHEM_DATA/sources/`, never in the repository.

```sh
python3 tools/get_sources.py --slugs tgirt-seq__10.1261+rna.054809.115
python3 tools/doctext.py /husky/henriksson/NGS101_papers/sources/tgirt-seq__10.1261+rna.054809.115/
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/tgirt-seq__10.1261+rna.054809.115/
```

The registered `pmccloud:PMC4691826` route now recovers the defining article
XML/PDF, original Supplemental Figures S1–S9 EPS files and data XLSX. The older
empty-bundle or PMC-gate reports do not describe this recovered set.
`MANIFEST.tsv` records each original URL and checksum.

`tools/doctext.py` extracts EPS text using optional Ghostscript followed by the
PDF text reader. Figure S1’s complete oligos are in
`PMC4691826_supp_supp_054809.115_FigS1.eps.txt`. Source-dependent checks skip
when that twin or its optional extraction tools are absent. Inspect the original
figure visually as well; text extraction is not itself a source-reading pass.

`catalogue/extra_sources.tsv` also registers the actual publisher-linked
`SuppLegends_Tables.docx` (publisher route currently 403), official NEB M0319
product/protocol URLs and the later PMC8678547 comparator. The NEB product and
adenylated-linker FAQ have readable browser scientific bodies. Actual protocol
opens currently yield title/footer shells, and an installed-browser attempt
reached a challenge page; their bodies remain pending. Record route-specific
responses separately from original-content reads.

The cached later detailed protocol is also retrievable with:

```sh
python3 tools/get_sources.py --doi 10.21769/BioProtoc.4239
```

It lives in `$CHEM_DATA/sources/doi_10.21769+BioProtoc.4239/`. This is a later
version, not the defining 2016 architecture. Cited predecessor original locations,
independent passes and named pending links are recorded in
`catalogue/source_reviews/tgirt-seq__10.1261+rna.054809.115.json`.
