# Source retrieval

Run `python3 tools/get_sources.py --slugs xyzeq__10.1126+sciadv.abg4755`. Originals and text twins belong only in the external source store.

- https://pmc-oa-opendata.s3.amazonaws.com/PMC8059935.1/PMC8059935.1.xml — published main; selected Methods and Results read.
- https://pmc-oa-opendata.s3.amazonaws.com/PMC8059935.1/PMC8059935.1.pdf — same published main; rendered Fig1A/B p3 inspected.
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8059935/supplementaryFiles — retrieved ZIP contains main figure images; not a substitute for detailed SI.
- http://advances.sciencemag.org/cgi/content/full/7/17/eabg4755/DC1 — exact main-printed supplementary index; standard request returned403. Actual SI original filenames remain unassigned.

No third-party originals are checked into this directory. Source download outcomes do not establish that information is absent.
