# Design decisions

This file records project-level choices that should remain visible and reviewable. It is
not a protocol reference and does not establish scientific facts. Dates record when a
decision was adopted; revise entries when the project changes direction rather than
letting implementation conventions become accidental policy.

## 2026-10-09 — Citation counts are dated discovery metadata

**Decision.** The protocol index may show citation counts for defining papers and filter
on them, but normal builds never query a citation service. A manually invoked script
resolves the DOI set against OpenAlex and replaces a checked-in snapshot containing the
count, OpenAlex work ID and retrieval date. The index labels the values as a dated
snapshot.

A missing paper or unresolved DOI is not zero citations. Such protocols remain visible by
default when a citation threshold is applied, with a separate control to exclude entries
without citation data. Commercial protocols without a defining paper are labelled as
commercial rather than assigned a synthetic count. Citation numbers rank discoverability,
not scientific quality or evidentiary authority.

**Why.** Citation count is useful as the catalogue grows into niche protocols, but it is
mutable external metadata and varies by provider. A deliberate snapshot gives useful
filtering without making builds network-dependent or pretending the number is live.
DOI-exact lookup avoids title matching.

**Review when.** Refresh the snapshot when maintainers want a newer view, or reconsider
the provider if its DOI coverage or count definition stops being useful. Do not schedule
automatic refreshes without revisiting this decision.

## 2026-10-09 — Identifier semantics belong to the molecular model

**Decision.** Barcode and UMI meaning is stored as format-neutral structured metadata on
the `Segment` that carries the bases. A molecular feature has a stable logical identity,
a controlled role, an encoding/source class, and optional multipart group, part and
whitelist reference. These semantics are independent of the displayed bases, placeholder
letters, colour tag and any external interchange format.

Construction and strand-transform helpers preserve feature identity. A final library's
identifier cycle ranges are derived from its annotated segments, sequencing-primer landing
sites, primer direction and declared run lengths. Pages render those computed locations;
they do not maintain a second hand-written read layout. Whitelist references identify an
authoritative list but do not require committing third-party sequence lists to this repo.

The repository does not infer semantics from names such as “barcode”, from `N`/`X`, or
from the legacy `cbc`/`umi` presentation colours. Those signals are ambiguous: the same
colours are also used for non-identifier chemistry. Unknown encoding is stated as such.

**Why.** Identifier positions are useful beyond any one exchange format: they support
read-cycle explanation, processing configuration, protocol comparison and validation.
Keeping them in the construct makes every renderer consume the same molecular truth and
prevents prose tables from drifting when a construct changes.

**Review when.** Add an exporter only if a useful target format is chosen. It must consume
this model rather than becoming another source of protocol data.

## 2026-10-09 — CRISPR-MIP analysis is a downstream repository

**Decision.** This repository retains the published CRISPR-MIP protocol, reviewed source
note, chemistry model and schematic. Project-specific work around that protocol — GC-bias
analysis, dataset investigations, proposed experiments and live troubleshooting — belongs
in the sibling `chem_crisprmip` repository.

The dependency runs one way: `chem_crisprmip` may use models and primitives from `chem`,
but protocol pages and site builds in `chem` must not depend on the analysis repository.

**Why.** The chemistry catalogue should remain a reusable explanation of published library
construction. Active analysis has different data, tooling, privacy and revision needs and
should not enlarge or leak into the public protocol site.

**Review when.** Move material back only if it becomes a reusable chemistry primitive or a
necessary part of accurately documenting the published protocol.

## 2026-10-09 — SVG is the default schematic renderer

**Decision.** Protocol diagrams render as SVG by default. Sequence and annotation content
remains SVG text, not outlined paths, so readers can select and copy it. Each SVG keeps an
intrinsic width derived from its character grid and sits in a horizontal scroll port;
narrow screens must scroll rather than shrink bases below the established font size.

Reaction and construct diagrams use SVG. Ordering-oligo lists remain raw HTML text with
explicit 5′ and 3′ ends, because direct selection and copying into an order form matters
more there than geometric presentation.

A reaction diagram is not a character-grid construct placed inside an SVG element.
Reaction pages build a `Workflow`: one initial molecular state followed by ordered
transitions. Each transition consumes the workflow's current state, so its input cannot
drift from the preceding output. The shared renderer draws every state once and connects
them with labelled vector arrows. Plain construct panels remain appropriate for final
libraries, primer binding and other states that are not transformations.

The molecular drawing grammar is also semantic. `Scene` emits typed strand, binding-span
and process-arrow geometry; the SVG renderer turns those into directional strand boxes,
brackets and vector arrows. It never recognizes `^^^^` or `<--------` after the fact.
Concrete named `Segment` regions carry hover text generated by the model, including a
known length and a calculated Tm for sufficiently long unambiguous DNA. Symbolic or
ellipsis-length regions do not claim a length. Explanatory labels use the
proportional interface font; bases remain monospaced and selectable.

A reaction panel stays aligned to the normal text column when its intrinsic content fits.
Only an oversized panel may expand toward the viewport edge; it scrolls horizontally once
the viewport itself is narrower than the construct.

The original character-grid HTML renderer remains available through
`panel(..., renderer="legacy")` for comparison and possible future uses, but is not
duplicated invisibly in SVG pages.

**Why.** SVG permits clearer molecular geometry without forcing long libraries into
unreadably small responsive images. Raw oligo text preserves the simplest ordering
workflow, while conditional breakout gives complex reactions room without misaligning
every short diagram from the prose.

**Review when.** Revisit the visual grammar as richer topology is added. Keep reaction SVG
text selectable, retain fixed-scale horizontal scrolling and the legacy fallback, and keep
ordering oligos as copyable text.

## 2026-10-07 — One schematic per protocol

**Decision.** A public schematic represents one protocol, using that protocol's defining
paper as its primary authority. Do not publish aggregate “family” schematics. Related
papers may be linked when they clarify or extend the protocol, but their chemistry must
not silently enter the defining construct.

Commercial protocols without a defining paper are allowed. In that case, the current
authoritative vendor manual is the primary source.

**Why.** A family page obscures which source supports a base, oligo, or reaction and makes
closely related methods look interchangeable when they are not.

**Review when.** Revisit only if the site gains an explicit comparison view. A comparison
view would sit above protocol pages; it would not replace them.

## 2026-10-07 — Schematic pages are diagrams, not research reports

**Decision.** Public protocol pages contain the source link, key oligos, construct-changing
steps, final library structure, and read layout. Literature comparisons, evidence logs,
performance results, unresolved alternatives, implementation trivia, and extended
reasoning belong in research notes. A short caveat remains on the schematic when it is
necessary to prevent an inferred structure from looking published.

The public front page is a compact protocol catalogue. It does not expose evidence logs,
build statistics, research-note summaries, or explanatory project prose. A finished
schematic may offer one quiet **Research notes** link to its main note.

**Why.** The useful unit for a reader is the molecular construct and how it changes. The
research trail remains valuable, but it should be available on demand rather than compete
with the diagram.

**Review when.** Revisit after observing real navigation needs, not to accommodate more
metadata on the front page.

## 2026-10-07 — Enforce invariants in the functions that build the model

**Decision.** Properties the repository intends always to hold must be enforced by the
construction path, rather than restated as protocol tests.

Examples:

- derive complementary strands with `complement_segments` or `revcomp`;
- draw multi-strand states with `Scene`, which rejects non-pairing columns;
- assemble adapters from canonical `illumina` / `nextera` components;
- construct recurring end-prep, hairpin and dumbbell states through `endprep`, `hairpin`
  and `dumbbell`, so overhang compatibility and topology are enforced once;
- construct chromosome ends and phased terminal adapters through `telomere`, so every
  repeat phase is derived from one motif and every telomere drawing retains its native
  duplex-to-3′-overhang geometry;
- locate sequencing primers through `seqprimers` while building the page, so a missing
  site makes the build fail;
- construct uncertain segments through helpers that apply `inferred=True`, rather than
  relying on every caller to remember the evidence styling.

Tests remain appropriate for source transcription, behavior at an API boundary, and
scientific claims that cannot be made true by construction. They should not merely repeat
lengths, temperatures, step counts, or prose already stored in the model and notes.
Shared primitive behavior is checked once in the repository-level suite; protocol suites
do not rerun the same shared checks.

**Why.** An invariant checked only by a test can still be bypassed by another caller. An
invariant encoded in the API keeps the logic at the point where the state is created and
makes every consumer safer by default.

**Review when.** When a test looks like a restatement of model content, first ask whether
the model API can make the invalid state unrepresentable.

## 2026-10-08 — Reaction junctions and sequencing-primer sites are part of the drawing

**Decision.** Every covalent ligation shown in a construct is marked at the exact boundary
with `** ligation`. Callers name the two adjacent segments through `Scene.junction` or
`junction_row`; the renderer refuses a missing or non-adjacent boundary.

Every finished library also shows its declared sequencing primers annealed to the final
duplex, including extension direction. `seqprimers.diagram` derives those positions from
the same `seqprimers.locate` result that validates the library. A protocol for which a
primer sequence or site is genuinely unavailable must say so explicitly instead of
silently omitting the binding diagram.

Human-readable bracketed role tokens are orientation-independent display text. Reversing
a strand may reverse molecular sequence, but must never render prose backwards.

**Why.** Junctions and primer sites are chemically important relationships, not captions.
Anchoring them to named model segments keeps them visible and makes diagram drift a build
error.

**Review when.** Revisit the marker glyph or layout if readability demands it; retain the
named-boundary and computed-primer requirements.

## 2026-10-07 — Evidence and uncertainty remain visible but secondary

**Decision.** Research notes retain the detailed evidence vocabulary (verbatim, derived,
and unpublished). On a schematic, an inferred molecular region requires all of:

1. a concise page-level caveat;
2. an `INFERRED — …` caption on the affected step;
3. inferred styling on the affected bases.

Do not add uncertainty prose to unrelated steps. Exact bases unavailable from an
authoritative source may be represented as length-preserving placeholders instead of
copying a secondary sequence into the main construct.

**Why.** Hiding the research trail must not turn a model into an unsupported assertion.
The three signals make the boundary clear without turning the schematic back into a
research report.
