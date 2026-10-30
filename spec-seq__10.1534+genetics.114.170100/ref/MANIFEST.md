# Source retrieval

Run `python3 tools/get_sources.py --slugs spec-seq__10.1534+genetics.114.170100`. Originals and texttwins remain in external source store;canonical manifest URL/version controls identity.

- Unavailable: https://academic.oup.com/genetics/article/198/3/1329/6065695 — Defining published Spec-seq original identified by primary publisher index; why: Generic DOI/PMC routes failed;verify returned printed title/DOI and Methods before use;no challenge bypass.; HTTP 403
- Unavailable: https://doi.org/10.1534/genetics.114.170100 — journal article; no open full text found (publisher page only)
- Unavailable: https://pmc.ncbi.nlm.nih.gov/articles/PMC4224169/ — PMC article page: full text; host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after host cool-down after reCAPTCHA challenge page (5.7 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (5.7 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (6.0 h left; not requested) (5.9 h left; not requested)
- Unavailable: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4224169/fullTextXML — full-text JATS XML (Europe PMC); HTTP 500

See [selected scope/deferred targets](../01_chem.md).
