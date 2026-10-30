# Source manifest

Originals/text twins live under `/husky/henriksson/NGS101_papers/sources/`. Run `python3 tools/get_sources.py --slugs miclip-nsun2-m5c__10.1016+j.celrep.2013.06.029` for the published article and mmc9 consolidated Supplemental Information. Exact URLs/hashes are in its external MANIFEST.tsv.

The directly cited König2010 published supplement is cached in external `iclip__10.1038+nsmb.1838/`: `nsmb.1838_41594_2010_BFnsmb1838_MOESM21_ESM.pdf`. Re-obtain with `python3 tools/get_sources.py --slugs iclip__10.1038+nsmb.1838`. Use the Springer published supplement, not the associated-manuscript counterpart. Original: https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnsmb.1838/MediaObjects/41594_2010_BFnsmb1838_MOESM21_ESM.pdf

No originals are stored in this checkout. Source-dependent checks skip when unavailable; selected scopes/limits are in the source review.
