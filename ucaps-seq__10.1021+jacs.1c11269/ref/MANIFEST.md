# Ucaps-seq source request

Run `python3 tools/get_sources.py --slugs ucaps-seq__10.1021+jacs.1c11269`.

Sources belong in `$CHEM_DATA/sources/ucaps-seq__10.1021+jacs.1c11269/` (default `/husky/henriksson/NGS101_papers/sources/`), never the checkout. The manifest records the defining DOI route and the declared independent published-PDF mirror. Neither returned a validated original in this pass; mirrorHTTP202 was rejected rather than used as a PDF. No downloaded originals or extracted text are committed.

Needed: published defining article/fullMethods and its actually linked chemistry Supporting Information. Exact SI download URL remains unresolved; no guessed URL is a source request. Identity/indexed excerpts are not a complete source reading. Run mandatory scraper again once originals/text conversions arrive, then inspect actual tables/figures. No absence claim or source-equivalence assumption.
