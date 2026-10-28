# Source manifest

Re-obtain registered originals with `python3 tools/get_sources.py --slugs hydro-trnaseq__10.1016+j.celrep.2017.07.029`. Originals and text conversions live in `$CHEM_DATA/sources/`, default `/husky/henriksson/NGS101_papers/sources/`. No third-party sources are kept here.

- article: https://pmc-oa-opendata.s3.amazonaws.com/PMC5564215.1/PMC5564215.1.xml — read. accepted author manuscript; not publisher VOR
- library_predecessor: https://pmc-oa-opendata.s3.amazonaws.com/PMC3508525.1/PMC3508525.1.xml — read. accepted author manuscript; not publisher VOR
- blocked_supplement_0: https://pmc.ncbi.nlm.nih.gov/articles/instance/5564215/bin/NIHMS894202-supplement-1.xlsx — unavailable.
- blocked_supplement_1: https://pmc.ncbi.nlm.nih.gov/articles/instance/5564215/bin/NIHMS894202-supplement-2.xlsx — unavailable.
- blocked_supplement_2: https://pmc.ncbi.nlm.nih.gov/articles/instance/5564215/bin/NIHMS894202-supplement-3.xlsx — unavailable.
- blocked_supplement_3: https://pmc.ncbi.nlm.nih.gov/articles/instance/5564215/bin/NIHMS894202-supplement-4.xlsx — unavailable.
- blocked_supplement_4: https://pmc.ncbi.nlm.nih.gov/articles/instance/5564215/bin/NIHMS894202-supplement-5.xlsx — unavailable.
- blocked_supplement_5: https://pmc.ncbi.nlm.nih.gov/articles/instance/5564215/bin/NIHMS894202-supplement-6.xlsx — unavailable.
- blocked_supplement_6: https://pmc.ncbi.nlm.nih.gov/articles/instance/5564215/bin/NIHMS894202-supplement-7.pdf — unavailable.

The explicitly inherited library original was fetched with `python3 tools/get_sources.py --doi 10.1016/j.ymeth.2012.07.030`; bundle `doi_10.1016+j.ymeth.2012.07.030`.

Exact inspected sections, source versions and deferred targets are in the main protocol note and corresponding `catalogue/source_reviews/` JSON. This is one capped pass, not an accepted source-absence audit.
