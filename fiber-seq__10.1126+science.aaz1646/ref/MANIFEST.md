# Source manifest

- Stergachis et al. 2020, DOI `10.1126/science.aaz1646`; obtain with
  `python3 tools/get_sources.py --doi 10.1126/science.aaz1646` (publisher full text may require manual access).
- Author-maintained bench protocol: `https://fiberseq.github.io/protocol/bench-protocol.html`.
- PacBio SMRTbell prep kit 3.0 procedure; fetched through the existing commercial protocol entry.
- GEO GSE146942 sample records (SOFT text, all GSM):
  `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146942&targ=gsm&form=text&view=brief`,
  saved by hand as `geo_GSE146942_gsm_brief.txt` in the source bundle with a `MANIFEST.tsv` row.
- Science Supplementary Materials (`aaz1646_stergachis_sm.pdf`):
  `https://www.science.org/doi/suppl/10.1126/science.aaz1646/suppl_file/aaz1646_stergachis_sm.pdf`;
  publisher returns HTTP 403 to the fetcher, so download by hand on an institutional network
  and place it in the source bundle (done 2026-10-17; MANIFEST.tsv row records md5).
