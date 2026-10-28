# Restart handoff

Snapshot: 26 October 2026, repository `/data/henriksson/github/chem`.
Read `AGENTS.md` and `DESIGN.md`, then inspect the current worktree before continuing.
The next substantial requested work is the eight-method batch in
`plans/cloning-examples-04.json`. The earlier ET source blocker remains unresolved.

## Immediate working state

At this snapshot, HEAD is `0ee7a45` (source-download helper). The preceding commit
`d62d4a9` contains the preceding protocol batch and computed index completeness filter.
The user makes the commits; no commit or push was made for the latest CI fix.

`tools/selftest_annealing.py` has a staged one-line fix: add `source-downloads` to
the temporary site's HTML fixtures. Preserve that staged change. This handoff is
also staged when saved. There are no ongoing tool processes to resume.

## Latest CI failure and fix

The Pages run for `0ee7a45` failed in the shared annealing self-test:
https://github.com/henriksson-lab/NGS_101/actions/runs/40406737576

The new mandatory standalone page was missing from the annealing test's temporary
site. `build_site.assemble({})` raised `FileNotFoundError` for
`source-downloads.html`. The fixture now includes the page. Production build
requirements were retained; the fix does not silently omit the helper.

Verified locally after the fix:

- `python3 tools/selftest_annealing.py`: 47 checks pass.
- `python3 tools/selftest.py`: 1,505 checks pass.
- Staged whitespace check passes.

The patch has not been pushed or verified by a new remote CI run. Useful local
logs, if still present: `/tmp/ci-annealing-fixed.log` and `/tmp/ci-shared-fixed.log`.
These temporary logs are convenience pointers, not permanent evidence.

## Completed index filter

The public index has a primary **Model completeness** filter, computed from
versioned graph-renderer reports rather than manual per-protocol ratings.
Implementation: `lib/model_completeness.py`, `lib/chemdraw.py`, `build_index.py`;
checks: `tools/selftest_model_completeness.py` and the shared suite.

The levels are no schematic, unassessed model, unresolved steps/material fates,
molecular/readout limits, inferred regions, and no declared model gaps. This is
model scope, separate from source review, chemistry-check coverage or experimental
validation. Malformed, missing and old reports remain unassessed. Variable inserts
and barcodes alone do not create a gap.

The index's search, combined filters, URL persistence/reload and Clear were tested
in Chrome. The feature's full site build passed: 448 diagram pages, 440 finished
public schematics and every relative link. Current counts should be recomputed if
protocols or model declarations change.

## Unlisted download page

`build_source_downloads.py` generates `source-downloads.html` offline from the
outstanding-only ledger `geturl.txt`. `build_site.py` builds and explicitly publishes
it. Direct deployed path: `/NGS_101/source-downloads.html`.

The page has no inbound site-navigation/index links or protocol-search entry and
has a `noindex,nofollow` directive. It remains publicly accessible by direct URL.
It publishes document links and request context, not source originals or the raw
ledger. Current input has 16 requests and 37 links. Remove fulfilled requests from
`geturl.txt`; do not maintain a second request list in the builder.

Ten helper checks and the site publication/link/asset audit passed. The source
helper addition is committed in `0ee7a45`; its missed annealing fixture is the staged
CI fix described above. Other site-assembly fixtures in `tools/selftest.py` already
include the new page.

## Previous protocol batch and source blocker

Seven public illustrative examples are implemented: BioBrick/RFC10, RF cloning,
circular TAR, original cloning CATCH, historical In-Fusion Dry-Down (2007), original
GoldenBraid (2011), and native-DNA CATCH nanopore (2018). The latter keeps its
historical MAP007 kit chemistry explicitly unassigned and separate from the
human-target amplified branch and current ONT kits.

RecE/RecT ET (`et-recombineering__10.1038+2417`) is private status `notes` with an
explicitly hypothetical neutral replacement graph. Do not promote it on the
strength of an abstract or matching intended DNA alone. The defining paper is
Zhang et al. 1998, DOI `10.1038/2417`. Publisher retrieval is inaccessible; an
identified Adrian Stewart author-uploaded PDF also returned HTTP403. Its fragmented
web preview did not verify the original substrate/topology branch. The source
cache contains only its manifest, not a retrieved original.

The earlier tracked goal is **blocked** on this source. It was not completed.
Its blocker recurred over three goal turns before that status was set. Avoid
repeating the same blocked fetches or bypassing publisher challenges.

Author-copy and publisher links are in `geturl.txt`, request 15. Request 16 is the
archived RFC10 text needed for formal revision comparison. The BioBrick defining
article and inspected OpenWetWare standard are already available. Other outstanding
requests remain in that ledger. Batch evidence:
`catalogue/curation_rounds/2026-10-24-cloning-examples-03.tsv` and the corresponding
`catalogue/source_reviews/*.json`. These bounded passes are not completed two-reader
absence reviews; no absence claims were accepted.

## User-authorized next batch

The user requested adding **all eight**, including the two lower-priority suggestions:

| Method | Defining DOI | Main modelling addition |
|---|---|---|
| Twin-Primer Assembly (TPA) | 10.1093/nar/gkx132 | PCR products and reannealing alternatives |
| NE-LIC | 10.1371/journal.pone.0057943 | Nicking and retained/displaced strands |
| MASTER ligation | 10.1093/nar/gkt122 | Methylation-dependent recognition and cleavage |
| MetClo | 10.1093/nar/gky596 | Methylation protection during hierarchical assembly |
| Original SSRTA | 10.1038/srep00141 | Source-specific integrase sites and recombination products |
| Cre/lox FLEx | 10.1038/nbt811 | Orientation-dependent inversion/excision and byproducts |
| Hot Fusion | 10.1371/journal.pone.0115318 | Source-specific assembly without added ligase |
| Golden Mutagenesis | 10.1038/s41598-019-47376-1 | Declared variant population and ambiguity |

The complete objective and requirements are saved in
`plans/cloning-examples-04.json`; implementation has not started. A new tracked goal
could not be created because the earlier blocked goal remains unfinished and the
goal tools do not support replacing it. Do not falsely complete ET to clear that
tracker. A fresh conversation may allow a new goal; otherwise the repository plan
still records the authorized scope.

Use neutral invented examples with structured `illustrative_example` metadata,
protocol-owned graphs and constructed-DNA endpoints. Preserve source-specific
implementations, molecular geometry, copying, methylation, rearrangement products
and unresolved host boundaries. Golden Mutagenesis needs an explicit variant
population without a coverage/performance claim. Exclude invented sequences from
annotation and sample-sheet references. Register catalogue/properties/citations,
source routes/reviews and outstanding requests; verify and stage completed work.

The user has authorized parallel workers in this protocol-curation work and asked
for a parallel worker on the completeness filter. Move to the next method when
source-blocked. Do not turn bounded source review into a fetcher/architecture
migration or treat unavailable material as confirmed absent information.

## Storage and build reminders

Originals and extracted text live under
`/husky/henriksson/NGS101_papers/sources/<slug>/`, never in this checkout. Use the
registered source tools, preserve document/version identity, run primer scraping
and verify actual source sections. Respect source-review caps and the two-reader
missing-information gate in `AGENTS.md`.

Keep Python/builds offline and standard-library-only. Generated HTML and `_site/`
are ignored output; downloaded originals must not be staged. Site delivery is via
CI after the user commits/pushes. Stage completed changes; do not commit unless
requested. Existing user preferences favor concise public pages and SVG figures.

## Resumed goal update on 27 October 2026

The user resumed the earlier eight-method goal. Its tracker was reactivated; the
ET blocker remained during the fresh blocked-goal audit. No new
original was found in its source cache or `up2/`. The author PDF again returned403;
the direct publisher PDF route returned HTML and failed PDF validation. The route,
source review, reference manifest and outstanding request were updated. No ET
chemistry was promoted. Preserve the previously staged CI fixture fix.

The resumed audit has now confirmed the same missing ET original over three
consecutive turns. The source cache still has only `MANIFEST.tsv`, `up2/` has no
new source files, and the last fetch completed with zero valid originals. All
other entries in the earlier eight-method goal retain their public example and
endpoint declarations. No source-qualified ET implementation can be completed
without the full published methods and figures. The blocked-status threshold
is met; do not mark this older goal complete or repeatedly retry blocked hosts.

## Restriction analysis tab added on 27 October 2026

The new `restriction.html` tab is implemented and ready for deployment after commit/push.
`build_restriction.py` builds it and the site-local catalogue; `build_site.py` includes
its page, worker and data. Shared navigation includes Restriction sites. Source and
UI code are `assets/restriction_analysis.js`, `restriction_worker.js`,
`restriction_ui.js` and `restriction.css`.

It supports unique/zero/arbitrary count ranges, separate cut/recognition counts,
linear/circular topology, both orientations, IUPAC ambiguity, type IIS/two-cut
coordinates, end-type/length filters, custom enzymes, context/flank views, reference
fragments and TSV/SVG export. CpG/Dam/Dcm sequence-context and sensitivity filters
are independent; sensitivity does not automatically alter cut counts. Unknown
sensitivity remains unknown. Exact-count defaults exclude uncertain extra cuts.

The checked-in reference snapshot is `catalogue/restriction_enzymes.json`: 1,088
Biopython 1.88 records, including 754 defined cut geometries, derived from REBASE 404
(2024), with pinned source SHA/revision, copyright and full permission notice.
`catalogue/restriction_methylation.tsv` has 36 exact-name factual annotations with
primary source URLs. No variant or isoschizomer inherits another annotation.
The refresh tool parses literal data without executing upstream Python or adding
Biopython as a dependency. Builds are offline. The generated
`assets/restriction_catalogue.json` is ignored output.

Validation: 45 Node matcher checks, 13 Python catalogue/build checks, 1,520 shared
checks, 47 annealing checks, existing browser sequence/samplesheet Node suites,
full 448-diagram site build/link audit and 25 actual Chrome smoke assertions pass.
`tools/benchmark_restriction_browser.py --site _site` repeats browser checks, including
catalogue-load retry, deliberately failed worker import, cancellation/retry,
GitHub Pages paths, CpG modes and exports. Worker errors are handled explicitly so
cancellation during import does not produce an uncaught browser error.
A local 100,000-base/1,088-enzyme matcher benchmark completed in about 3.3 seconds.

All prior staged work, including the earlier CI fixture fix and ET source evidence,
is preserved. The annealing and publication fixtures also include the new tab and
its explicit worker/data assets. The older ET goal remains blocked; the next
protocol batch in `plans/cloning-examples-04.json` remains unimplemented.

### Interactive restriction display — 28 Oct

The restriction tab now opens with a clickable whole-molecule overview and a
zoomable paired-sequence viewport. Recognition spans and top/bottom cut boundaries
come from the matcher; the aligned column/polarity template comes from
`Scene.duplex`, and the existing browser IUPAC helper derives the complement.
Select a span, enzyme or site; pan or choose a 40/80/160 bp window; export the
current duplex as standalone SVG. The enzyme table is retained under a collapsed
“Enzyme table” disclosure. Circular-origin, outside-recognition and incomplete
linear-end geometry retain their existing coordinate meaning. Dense displays and
retained-position limits are labelled.

New asset: `assets/restriction_view.js`, loaded explicitly by the restriction page
and copied through the usual site asset walk. Verified: 45 matcher checks, 13
catalogue/build checks, 1,520 shared checks, and 34 actual Chrome checks against the
assembled site, including complement alignment, both EcoRI cut boundaries, span
selection, pan/zoom, reverse type IIS, circular sites and standalone SVG export.
All assembled-site links/publication checks passed. Browser helper accepts an
optional `--screenshot PATH` for visual inspection. Changes staged; no deployment.

### Restriction viewer revised to full-sequence scrolling — 28 Oct

Supersedes the earlier overview/window design above. The overview, pan/zoom/slider
controls and tall per-enzyme tracks are removed. The query wraps into continuous
paired rows in a native scroll area; virtualized rows preserve the entire sequence
and its absolute coordinates without creating 200,000 base elements for a 100 kb
query. Motifs sharing recognition and cut geometry combine into compact labels.
Hover labels or bases (including overlapping motifs) to see associated cuts and
coordinates; click a label to pin, or clear the selection. Tooltip sits above the
hover target so it does not cover the cut. Every cut path has vertical strand marks
and a horizontal middle segment. The matcher now retains raw (unwrapped) cut
coordinates alongside normalized coordinates to preserve origin-crossing geometry.
The collapsed enzyme table, filters and TSV exports remain; SVG exports visible
sequence rows and active cut geometry.

Verified 45 matcher checks, 13 catalogue/build checks, 1,520 shared checks, and 37
Chrome checks, including actual pointer hover, horizontal connectors, circular
origin splits, reverse type IIS sites, and native scrolling to a motif after base
10,000 with bounded DOM size. Site assembly/link/publication checks pass. No deploy.
