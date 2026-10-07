# Design decisions

This file records project-level choices that should remain visible and reviewable. It is
not a protocol reference and does not establish scientific facts. Dates record when a
decision was adopted; revise entries when the project changes direction rather than
letting implementation conventions become accidental policy.

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
