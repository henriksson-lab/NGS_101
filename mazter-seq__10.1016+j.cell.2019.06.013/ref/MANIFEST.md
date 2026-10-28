# MAZTER-seq source retrieval

Originals and extracted text live only under `$CHEM_DATA/sources/mazter-seq__10.1016+j.cell.2019.06.013/` (default external source store).

```sh
python3 tools/get_sources.py --slugs mazter-seq__10.1016+j.cell.2019.06.013
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/mazter-seq__10.1016+j.cell.2019.06.013/
```

`catalogue/extra_sources.tsv` declares the institutional publisher-PDF original and the publisher article/supplement index. The former fetched and was read at the sections recorded in the note/review; the latter returned403. Table S7 is explicitly referenced, but its exact attachment URL has not been resolved. Failed superseded-preprint retrievals are excluded from the current blocker queue. No downloaded originals or text twins are stored here.
