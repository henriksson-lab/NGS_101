# Source manifest

Defining published original and downloads: `$CHEM_DATA/sources/arm-seq__10.1038+nmeth.3508/MANIFEST.tsv`, default external root `/husky/henriksson/NGS101_papers`.

Re-obtain with `python3 tools/get_sources.py --slugs arm-seq__10.1038+nmeth.3508`. No third-party originals or extracted text are kept here.

- `nmeth.3508_PMC4553111.pdf`: https://pmc-oa-opendata.s3.amazonaws.com/PMC4553111.1/PMC4553111.1.pdf — inspected Abstract28–35; Introduction85–108; Discussion272–297; Methods315–402; processing403–426; differential458–465.
- `nmeth.3508_41592_2015_BFnmeth3508_MOESM407_ESM.pdf`: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3508/MediaObjects/41592_2015_BFnmeth3508_MOESM407_ESM.pdf — inspected Figure1caption1–13; Figures2–5captions15–65.
- `nmeth.3508_41592_2015_BFnmeth3508_MOESM408_ESM.xlsx`: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnmeth.3508/MediaObjects/41592_2015_BFnmeth3508_MOESM408_ESM.xlsx — inspected First sheet heading and opening rows only.

Exact positive findings and deferred sources are in `01_arm_seq.md` and the corresponding `catalogue/source_reviews/` JSON. This is one capped pass, not two independent absence readers.
