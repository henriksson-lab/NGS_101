# Source manifest

Run `python3 tools/get_sources.py --doi 10.1038/nmeth.3960`. PMCID PMC5045719 and its PMC
Cloud package contain the author manuscript, figures, Supplementary Figures and
Supplementary Tables. Do not commit downloaded source files.

## Wanted: a text source for the oligo set

Every DSBCapture oligo is published only as raster artwork inside Fig. 1a and
Supplementary Fig. 1c. Nothing fetched above contains a single base of it: `grep -oiE
"[ACGT]{15,}"` over all extracted `.txt` twins returns no hit, and `strings` over the raw
`supplementary_figures.doc` finds no oligo, because the figure is an embedded image. The
bases this repository models therefore have no machine-checkable witness, which
`lib/dsb_capture.py` records per graph material in `SOURCE_OLIGO_EXCEPTIONS` and
`01_dsbcapture.md` explains under "Oligo provenance".

**The missing document is the authors' step-by-step protocol.** The Methods cite it twice
with an unfilled placeholder — "A detailed protocol can be found on the protocol exchange
(Insert URL)" and "the step by step method on the protocol exchange (Insert URL)" — so the
paper as published does not identify the entry. No source on disk names it, and the only
URL the article does give, `http://breakome.utmb.edu`, is cited for the BLESS method
DSBCapture was adapted from, not for DSBCapture. A guessed Protocol Exchange URL is
therefore not added to `catalogue/extra_sources.tsv`.

If the entry is located (Lensing et al. 2016, *Nature Methods* 13:855; Protocol Exchange
entries carry `10.1038/protex.*` DOIs), add it there as
`dsbcapture__10.1038+nmeth.3960 <tab> doi:10.1038/protex.… <tab> Authors' step-by-step
DSBCapture protocol <tab> Only possible text source for the modified P5/P7 duplexes, the
control modified P5 adapter and the PCR primers.` Recovering it would convert the oligo
set, the index length and the control adapter from figure-image evidence to text, let the
`SOURCE_OLIGO_EXCEPTIONS` entries be deleted, and let
`tools/selftest_source_oligos.py` check the transcription for real.
