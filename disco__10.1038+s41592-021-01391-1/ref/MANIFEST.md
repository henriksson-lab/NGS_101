# Reference manifest

- Defining source: https://doi.org/10.1038/s41592-021-01391-1

Retrieve with `python3 tools/get_sources.py --doi 10.1038/s41592-021-01391-1`.
- Author step-by-step protocol (library oligos, reactions, read lengths):
  https://www.epfl.ch/labs/deplanckelab/wp-content/uploads/2021/09/DisCo_protocol.pdf
  (fetched by `python3 tools/get_sources.py --slugs disco__10.1038+s41592-021-01391-1`
  via its `catalogue/extra_sources.tsv` row).
- Superseded bioRxiv preprint (lead only, not used by the model):
  `python3 tools/get_sources.py --doi 10.1101/2020.05.19.103812`.
