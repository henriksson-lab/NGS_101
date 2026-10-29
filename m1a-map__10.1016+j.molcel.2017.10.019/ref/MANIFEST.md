# m1A-MAP source retrieval

Originals and extracted text remain outside the checkout in `$CHEM_DATA/sources/m1a-map__10.1016+j.molcel.2017.10.019/`.

```sh
python3 tools/get_sources.py --slugs m1a-map__10.1016+j.molcel.2017.10.019
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/m1a-map__10.1016+j.molcel.2017.10.019/
```

Canonical manifest identifies associated accepted-manuscript PMC5722686 XML, not a verified whole-VOR equivalent. Methods S1 is explicitly identified as `NIHMS915125-supplement-2.pdf` and its exact PMC request reached a download gate. Other actual attachments remain in the manifest/checklist; superseded preprint failures are excluded. Reading sections and pending linked originals are recorded in the note/source review. No downloaded original or text twin is checked in.
