# Source manifest

Retrieve the defining original and registered linked originals with:

```sh
python3 tools/get_sources.py --slugs dm-trna-seq__10.1038+nmeth.3478
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/dm-trna-seq__10.1038+nmeth.3478/
```

Downloads remain under `$CHEM_DATA/sources/dm-trna-seq__10.1038+nmeth.3478/`;
`MANIFEST.tsv` records actual URLs, checksums and failed routes. Repeatable extra
routes are in `catalogue/extra_sources.tsv`: cited Mohr/Shen/Katibah originals,
InGex TGIRT10, the author-lab later small-RNA/CircLigase PDF, distributor-mirrored
original Epicentre CircLigase II manual, current official HiSeq 2500 primer page,
and the CDC Stacks original Shen supplement.

Katibah's readable main article was obtained separately by
`python3 tools/get_sources.py --doi 10.1073/pnas.1412842111`, under
`$CHEM_DATA/sources/doi_10.1073+pnas.1412842111/`. Its actual SI link is registered;
HTTP and HTTPS publisher attempts returned 403 and PMC attachments returned a
download gate. These are route-specific failures, not a completed absence or
all-route unavailability conclusion. Author-copy recovery remains pending.

Source-specific locations, independent readers and open links are recorded in
`catalogue/source_reviews/dm-trna-seq__10.1038+nmeth.3478.json`. Neither source
files nor their extracted text are committed.
