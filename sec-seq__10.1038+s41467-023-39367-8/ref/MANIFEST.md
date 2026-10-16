# Sources

Retrieve the article and supplements with `python3 tools/get_sources.py --doi 10.1038/s41467-023-39367-8`.

Official CG000317 Rev E is repeatably fetched into the SEC-seq source bundle with
`python3 tools/get_sources.py --no-papers --slugs sec-seq__10.1101+2022.08.25.505190`.
The source catalogue retains the preprint slug; the public protocol directory
uses the published defining DOI. The guide URL is recorded in
`catalogue/extra_sources.tsv`. Appendix pp. 71–72 supplies both dual-index
endpoint families used by the main human ASC model.
