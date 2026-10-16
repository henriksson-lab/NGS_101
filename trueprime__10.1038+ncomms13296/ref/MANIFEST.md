# Reference manifest

- Picher et al., *Nature Communications* (2016), [doi:10.1038/ncomms13296](https://doi.org/10.1038/ncomms13296), PMCID PMC5141293.
- Re-obtain with `python3 tools/get_sources.py --doi 10.1038/ncomms13296` or `pmccloud:PMC5141293`.
- The later [4basebio TruePrime WGA handbook](https://www.4basebio.com/s/3x0xxx_TruePrime-WGA-Kit-1.pdf) corroborates purified-DNA handling but does not establish the 2016 kit revision. The `catalogue/extra_sources.tsv` route retrieves it with `get_sources.py`.
- The paper-linked [ENA study SRP085855](https://www.ebi.ac.uk/ena/browser/api/xml/SRP085855), run table and four experiment XMLs are retrieved by the `catalogue/extra_sources.tsv` routes; they do not name the downstream NEBNext kit.

No third-party source files are stored here.
