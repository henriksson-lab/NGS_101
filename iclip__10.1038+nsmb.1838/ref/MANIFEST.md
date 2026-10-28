# Source manifest

Defining published original and downloads: `$CHEM_DATA/sources/iclip__10.1038+nsmb.1838/MANIFEST.tsv`, default external root `/husky/henriksson/NGS101_papers`.

Re-obtain with `python3 tools/get_sources.py --slugs iclip__10.1038+nsmb.1838`. No third-party originals or extracted text are kept here.

- `nsmb.1838_PMC3000544.pdf`: https://pmc-oa-opendata.s3.amazonaws.com/PMC3000544.1/PMC3000544.1.pdf — inspected Introduction71–83; Results87–110; Discussion288–296.
- `nsmb.1838_41594_2010_BFnsmb1838_MOESM21_ESM.pdf`: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnsmb.1838/MediaObjects/41594_2010_BFnsmb1838_MOESM21_ESM.pdf — inspected Supplementary Methods2048–2195; original PDF15–16 visually inspected.
- `nsmb.1838_supp_NIHMS31547-supplement-1.pdf`: https://pmc-oa-opendata.s3.amazonaws.com/PMC3000544.1/NIHMS31547-supplement-1.pdf — inspected Corresponding Methods2024–2155.

Exact positive findings and deferred sources are in `01_iclip.md` and the corresponding `catalogue/source_reviews/` JSON. This is one capped pass, not two independent absence readers.
