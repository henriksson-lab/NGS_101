# Reference manifest

- Defining source: https://doi.org/10.1038/s41587-021-01065-5

Retrieve with `python3 tools/get_sources.py --doi 10.1038/s41587-021-01065-5`.

That call writes to `sources/doi_10.1038+s41587-021-01065-5/`; the published files
(PMC10549956 JATS XML, Springer MOESM1–3) were moved into
`sources/phage-atac__10.1038+s41587-021-01065-5/` beside the superseded bioRxiv v1
files, with their MANIFEST.tsv rows. The PMC "Integrated supplementary information" PDF
is download-gated (manual row); MOESM1 is the publisher copy read instead.
