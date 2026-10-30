# BAsE-Seq original-source retrieval

Run `python3 tools/get_sources.py --slugs base-seq__10.1186+s13059-014-0517-9`.
The standard source store is `/husky/henriksson/NGS101_papers/sources/base-seq__10.1186+s13059-014-0517-9/`.

The external `MANIFEST.tsv` records exact URLs, checksums and download identities.
Published original: DOI10.1186/s13059-014-0517-9, PMC4269956; Cloud JATS/PDF.
Published supplements: `13059_2014_517_MOESM1_ESM.pdf` and
`13059_2014_517_MOESM2_ESM.pdf`. The second is the detailed protocol printed
last updated2July2014; verify that revision when retrieving it again.

Run `python3 tools/scrape_primers.py` on that external directory. Read source-line
contexts and inspect original tablep3 for underlining, the red AT insertion and
physical modifications. The scraper's truncated internal-biotin parse is not the
modified primer sequence. No upstream drawn page supplies an oligo or reaction.

The historical NexteraXT Sample Preparation Guide is cited by detailed Step8;
its exact source revision/applicability remains unresolved. A current guide is not
a silent replacement. No requested or unavailable oligo document is invented.
No third-party originals or extracted text are stored in this checkout.
