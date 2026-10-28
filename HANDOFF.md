# Restart handoff

Snapshot: 28 October 2026, repository `/data/henriksson/github/chem`.
Read `AGENTS.md` and `DESIGN.md`, then inspect git status and the active goal.
The user handles commits; do not commit or push.

## Paused work

The active goal expands the NGS collection toward at most 1,000 actual protocols,
or stops earlier at defensible candidate exhaustion. It is paused at the user's
explicit request; do not restart curation without explicit resumption. Do not
restore the old cloning/ET goal or treat a capped source pass as a full review.
`plans/ngs-expansion.json` records baseline, counts, admission policy and rounds.
Use parallel disjoint protocol workers as explicitly authorized by the user.
Keep shared catalogue, citation, round-ledger and progress edits with the coordinator.

Baseline commit: `fe8eac0146c54e557fbc2f91c8bfb37856108c48`.
The baseline had 451 total entries, including 414 NGS entries (413 published,
408 finished published). Five bounded rounds added46 source-read NGS methods
(460 collected including WIP,459 published). One more registered name,DNTR-seq,
is identity-only and excluded from that count. Catalogue has498 entries total,
461 registered NGS names including WIP. New43 source-read additions stay notes-only.
All have source-reading notes and local graphs; GLOE-seq native-end single-index,
bulk ACT-seq and the NHF1/CS-B XR-seq branch have finished schematics (411 finished published NGS). Other additions remain notes-only.
Status, computed model completeness and a completed source/absence audit are
separate. No new accepted-missing claims.

The three `2026-10-28-expansion-discovery-*.tsv` ledgers hold screened candidates;
they do not count as collected. Select the next bounded batch from them after
reviewing the five current round ledgers and their deferred targets. Prioritize
breadth, 10 minutes/three originals per protocol, published originals first,
mandatory primer scraper and source verification. Move past blocked sources.
Do not recursively chase citations or fill chemistry from related protocols.

## Pause checkpoint

The user requested pausing the goal, stabilizing the checkout and staging files.
Verified work through round05 is preserved. Round06 selected13candidates and
started source retrieval, but no protocol directories/models/review records were
created. Its provisional catalogue/properties/citation additions were restored
to the verified staged round05 versions so no empty notes entries remain.
`plans/ngs-expansion.json` saves the exact selected registration records and
preliminary reading pointers; the round06 ledger marks them paused/unregistered.
Do not count them as collected. On explicit resumption, re-register and refresh
catalogue/citations before continuing bounded source checks.

All root/worker fetch processes are terminal. External caches remain intact.
Xdrop retains an interrupted `MANIFEST.tsv.part` beside its ordinary manifest;
validate source identities/manifests before cache reuse on resumption. No temporary
source file or downloaded original is staged. NOMe source retrieval completed,
but no NOMe chemistry was read or scraped. FLEP/CLASH and DNaseHi-C received only
the preliminary inspections listed in the paused ledger; no conclusions saved.

## Relationships feature

User also requested method ancestry, multiple parents, inferred links, sequence
similarity, a full graph tab and local views. Implemented in:

- `catalogue/protocol_lineage.json`: qualified curated relationships.
- `lib/protocol_lineage.py`, `build_lineage.py`: validation, similarity and views.
- `assets/protocol_lineage.js` / `.css`: pan/zoom, focus/search, evidence, SVG export.
- Site assembly inserts a collapsed local graph on eligible protocol pages.

Ancestry never supplies molecular constructors. Exact known sequence-component
similarity is an optional independent layer, with inverse-frequency weighted
Jaccard and strand deduplication; it never becomes ancestry automatically.
Only finished public non-illustrative protocols appear. Notes-only edges remain
validated but private. Curation coverage is partial, not evidence of independent
origins. Follow `MAINTAINERS.md` when adding relationships during source reading.

Private iCLIP→irCLIP inference is explicitly accepted-manuscript-based. Its
2 h versus publisher-figure1 h click timing is a source-version comparison;
ordinary revision differences have not been excluded. Do not call it a confirmed
publication error. Other source/revision limits live in each note/review JSON.

## Validation and staging

Validation logs are transient `/tmp/ngs-expansion-*` and
`/tmp/protocol-lineage-browser.json`; the durable feature report is
`plans/protocol-lineage-implementation.json`. The full viewer passed 21 semantic
and 17 real Chrome checks, including prefix hosting, evidence filters, multi-parent
focus, local/full links, SVG export, failed data loading and mobile layout.
Shared suite, catalogue, source-tool checks and full assembled-site link checks
are recorded in the expansion plan as completed. Inspect final git status for
staging; generated HTML, `_site/` and source originals must remain untracked.

Source originals live at `/husky/henriksson/NGS101_papers/sources/`.
No third-party materials belong in the checkout. Recent source-specific model
gaps include carbon-spacer polymerase stopping (irCLIP), AtRNL-compatible
ribonucleotide circularization (Ribose-seq), reversible G-adduct occupancy
(KAS-seq), and Cas-bound exonuclease protection (CaBagE). Local typed declarations
are source-scoped and do not claim unsupported reaction/yield validation.

## Round03 integration

Added HydraPsiSeq, hydro-tRNAseq, XR-seq, Nick-seq, QuRIE-seq, DAb-seq and
PhIP-Seq; only XR is newly public. PhIP's detailed author supplement remains
provisional against the final article; its exact VOR endpoint is now queued in
extra_sources.tsv, not declared inaccessible. The other new notes preserve
kit/end/read-layout boundaries and source discrepancies.

The local irCLIP C18 backbone now supports explicit closure and conditional
first-extension traversal with22 passing checks. It never flattens spacers into
DNA or assumes rolling-circle copying. Full two-primer PCR/readout stays deferred.
Both iCLIP and irCLIP private builders use the supported literal OUT assignment.

Source guard whole-strand fallback fixes mixed known-N/unknown-width barcode
windows, with6 new positive/negative regression checks;1592 shared checks pass.
Current catalogue84/source-tools227 checks pass. All464 builders were verified,
then assembly/link checks completed after fixing the new private irCLIP OUT
parser integration;443 public pages and411 documented NGS graph contracts pass.
No global shared reaction claim validates the conditional spacer mechanism.

Select further candidates from the three discovery ledgers; e.g. POINT/nano-COP,
HydEn-seq/emRiboSeq and LIBRA-seq/scID-seq remain eligible for another bounded
round, after checking current registration. Do not count them before source notes.

## Round04 integration

Twelve additional source-read models/notes, plus DNTR identity-only source-blocked
notes with no graph. Registered NGS total447;collection count446. All12 new
models remain notes-only. The original emRiboSeq uses Ion Torrent (new valid
platform facet) with instrument emPCR,not an invented library-PCR series.

CPT's local ContiguityState preserves source protein-held fragment association
through pooling/dilution,then SDS dissociation in a physical well. Corrected
Table4 supplies20 transferred oligos and4 other stock/run oligos. Thirty-one local
checks pass. Four code components use distinct IDs with a common group;custom
physical index-cycle mapping and historical PCR primers remain unresolved.

irCLIP multipart IDs were also corrected:W3/N3/N5 have unique IDs sharing one
UMI group. Actual C18 nodes stay present and Construct compatibility is checked;
24 local checks pass. Full PCR/readout and final-VOR applicability still deferred.

PhIP's final publisher HTML is subscription-preview only,and its exact PDF URL
returns HTML. Its source review is now blocked;the existing source-read author
model is still counted as provisional notes. DNTR has no chemistry originals
and is explicitly uncollected in plans/ngs-expansion.json. Geturl adds9 groups
with22 actual blocked documents;none of16old groups was verifiably fulfilled.

Validation:1614shared/84catalogue checks;472current diagram builders,443public
pages,411NGS contracts;all relative links resolve. Scoped round04sourceguard
22checked/0failed,1source-notation exception,21unrepresented inputs. Source
reviews across new33registered names:30in_progress,3blocked,noaccepted absence.

The source queues still contain defensible candidates. Select another bounded
round after checking registration;ID-seq (10.1038/s41467-018-04761-0) was added
as a citation-identity-only discovery lead,not admitted/curated. Do not treat
empty initial candidate queues as proof of global NGS-method exhaustion.

## Round05 integration

Fourteen further source-read methods, all notes-only: m6ACE, BID, Nm, CeU,
RRBS2008, MiGS, FAIRE, Boyle DNase, Han TCR/phenotype, epicPCR, targeted DART,
TAP, MAZTER and m6A-REF. Collection460 incl WIP; registered461 because DNTR
remains identity-only/uncollected. No new public promotion or accepted absence.

Ten new private previews render model-derived workflows without completion
markers. Four RNA models expose scoped graphs/notes without a preview. CeU
has only source-inspected probe label/click events; its full Methods remain
blocked. Actual Han corrigendum and Boyle chemistry supplement blockers stay
explicit. MAZTER TableS7 is source-referred but attachment identity unresolved.

Independent implementation review fixed two omitted source heat operations
(MAZTER RNA+adapter preligation and m6A-REF first pre-FTO heat), plus MiGS
conjugate washing and possible blocking-DNA persistence through protease/cleanup.
Source-version differences, exact unassembled library geometries and chemical
versus computational fates remain explicit. The fourteen graphs have8established
narrow pairing relations and166unresolved claims;0contradictory does not mean
complete chemistry validation. Source guard15exact checked,0failed,0exemptions,
34inputs without exact states. All local suites pass;1647shared/84catalogue/
227source-tool checks pass. Final full-site result/staging is recorded in the plan.

TAP's published [X]8 notation now normalizes narrowly to exactly8variable positions
in source_oligo_check, with seven semantic regressions. No protocol was changed
to satisfy matching. Existing whole-strand gap and fixed-flank protections remain.

Lineage registry22 relations, including5new private scoped component adaptations;
public8links unchanged. RNA discovery ledger adds two unread leads from actual
references: DNA method10.1093/nar/gky1104 and RBS-seq10.1073/pnas.1817334116.
Do not count or register them without boundary/source review. Remaining verified
candidates in all three ledgers are still available; global exhaustion is unproven.
