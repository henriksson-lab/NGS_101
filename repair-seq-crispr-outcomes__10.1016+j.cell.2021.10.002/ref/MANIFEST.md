# Reference retrieval

Originals and text twins live outside the checkout in `$CHEM_DATA/sources/repair-seq-crispr-outcomes__10.1016+j.cell.2021.10.002/`. Re-run `python3 tools/get_sources.py --slugs repair-seq-crispr-outcomes__10.1016+j.cell.2021.10.002` and mandatory `tools/scrape_primers.py` on that bundle.

- j.cell.2021.10.002_PMC9074467.xml — https://pmc-oa-opendata.s3.amazonaws.com/PMC9074467.1/PMC9074467.1.xml — AssociatedNIHMS1750972identity; selectedMethods264–270/292, UMIs106 and supplementalTableS2role469–477.

Actual needed sources:

- https://pmc.ncbi.nlm.nih.gov/articles/instance/9074467/bin/NIHMS1750972-supplement-9.xlsx — ActualTableS2oligonucleotide/TaqManworkbook(mainlabels9/TableS2) neededforwholeUMIadapter/primerorderings; standardfetch gated/TDMnotCloudobject.

Other automatic downloads, including superseded preprints, were not adopted as recipe authority. No originals belong in this checkout.
