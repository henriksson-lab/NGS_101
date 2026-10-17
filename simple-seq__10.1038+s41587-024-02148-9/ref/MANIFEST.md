# Source manifest

Retrieve sources with `python3 tools/get_sources.py "SIMPLE-seq"`.

The authors' bench protocol is not yet covered by a fetcher; download it by hand from
the code repository linked in the article's Code availability statement:
<https://github.com/cxzhu/SIMPLE-seq/blob/main/Protocol/Protocol_SIMPLE-seq_Feb2024.pdf>
(commit `13e187f`, 8 pages, md5 `51d76194653a99a270e04adf323e5985`), saved as
`github_cxzhu_SIMPLE-seq_Protocol_SIMPLE-seq_Feb2024.pdf`, then run `tools/doctext.py`.
The main article is subscription-only (no PMC copy). The publisher PDF was hand-downloaded
on the institutional network from <https://doi.org/10.1038/s41587-024-02148-9> and saved as
`s41587-024-02148-9.pdf` (md5 `892edd12cf8e048df548749fbedd3b0c`, 11,631,024 bytes,
printed DOI checked) with a `MANIFEST.tsv` row; then run `tools/doctext.py`.

GEO record text (authors' short protocol and read lengths) was hand-fetched as SOFT
text: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<ACC>&targ=self&form=text&view=full`
for GSM5929312, GSM5929313 and GSM7291623 (and `view=brief` for GSE197740), saved as
`geo_acc_<ACC>_self_<view>.txt` with `MANIFEST.tsv` rows.
