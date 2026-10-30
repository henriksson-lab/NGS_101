# Source retrieval

Run `python3 tools/get_sources.py --slugs poly-a-seq-yu2020__10.1371+journal.pone.0234696`. Originals and extracted twins live under `/husky/henriksson/NGS101_papers/sources/poly-a-seq-yu2020__10.1371+journal.pone.0234696/`, outside the checkout.

- journal.pone.0234696_PMC7297374.xml — Published final PMC7297374/PLOS article; inspected Methods77–99, Results126–130/299: initial and later platform preparations, spike timing, nativeRNA versus copiedDNA/readout distinctions.
- journal.pone.0234696_supp_pone.0234696.s001.xlsx — PublishedS1Table; inspected Header and firstsample rows only: own-study mapping statistics, not an oligo table.
- journal.pone.0234696_supp_pone.0234696.s002.xlsx — PublishedS2Table; inspected Header and firstsample rows only: downloaded comparison-data mapping statistics, not a library recipe.

Source checklist and cap are recorded in `catalogue/source_reviews/poly-a-seq-yu2020__10.1371+journal.pone.0234696.json`. Filenames alone are not provenance.
