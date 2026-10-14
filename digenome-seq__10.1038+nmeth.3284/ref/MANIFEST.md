# Digenome-seq source manifest

No third-party source files are committed.

Retrieve the fetchable sources with:

```sh
python3 tools/get_sources.py "Digenome-seq"
```

`catalogue/extra_sources.tsv` registers four of them:

| Source | Why |
|---|---|
| Springer supplementary PDF for the defining paper | Supplementary Notes, 5′-end plots and Supplementary Tables 1–6. It contains no library-preparation text. |
| `pmccloud:PMC4772022` — Kim et al., *Genome Research* (2016) | The open-access companion; the only readable source that states the fragmentation, input and sequencing depth. |
| Illumina TruSeq DNA PCR-Free Reference Guide 1000000039279 v00 | Vendor authority for the TruSeq reaction order and the index-adapter components offered. |
| Illumina Adapter Sequences 1000000002694 v22 | Canonical TruSeq adapter, sequencing-primer and index bases; fixes the single-index length at six. |

Two sources remain out of reach and are **not** registered:

- **The defining article itself**, Kim et al., *Nature Methods* (2015). Subscription-only
  and not in PMC; `MANIFEST.tsv` keeps it as a `(manual)` row with its DOI. Only its
  supplement is open.
- **Kim, Kang and Kim, *Nature Protocols* (2021), doi:10.1038/s41596-020-00453-6.**
  Paywalled, no PMCID, no legitimate fetch route. Claims in `01_digenome_seq.md` that
  depend on it are marked 🟡. If a copy is obtained by hand, drop the PDF into
  `$CHEM_DATA/sources/digenome-seq__10.1038+nmeth.3284/` and re-run
  `python3 tools/get_sources.py --slugs digenome-seq__10.1038+nmeth.3284 --revalidate`.

The note therefore identifies which claims come from the defining paper, which are
corroborated by the authors' open follow-up, and which remain inferred.
