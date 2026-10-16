# Sources

Retrieve registered originals without committing them:

```sh
python3 tools/get_sources.py --slugs gro-seq__10.1126+science.1162228
python3 tools/scrape_primers.py /husky/henriksson/NGS101_papers/sources/gro-seq__10.1126+science.1162228
```

- Defining article: https://doi.org/10.1126/science.1162228 (PMC2833333).
- Original Supporting Online Material: registered public mirror
  https://www.cecerelab.com/_files/ugd/594009_9144e1ea58ad4573871ec024827290f3.pdf
  (cover explicitly identifies Core/Waterfall/Lis, Science 2008, DOI and DC1).
  The PMC attachment route was gated; this does not negate the recovered original.
- Later detailed GRO-seq protocol PMC5880289 is a separate implementation, not the
  source of historical 2008 adapter orders or poly(A)/circularization operations.

Downloaded originals and extracted text stay under `$CHEM_DATA/sources/`. The review JSON
records exact independent readings, former availability history and remaining source scope.
