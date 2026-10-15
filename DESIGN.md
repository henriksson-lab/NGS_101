# Design decisions

This file records project-level choices that should remain visible and reviewable. It is
not a protocol reference and does not establish scientific facts. Dates record when a
decision was adopted; revise entries when the project changes direction rather than
letting implementation conventions become accidental policy.

## 2026-10-15 — Permit one-module protocol pages

**Decision.** A protocol module may declare a literal `PAGE_FILE` and expose
`protocol_graph()`, `render_page()` and `selftest()`. The central `tools/protocol.py`
command builds and tests it; site, index and catalogue discovery accept this alongside
existing local `build_page.py` and `selftest.py` scripts. A protocol with substantial
source checks may retain a separate self-test. Only the local protocol module chooses
its source-specific reaction graph and page emphasis.

**Why.** SMART-seq2, sci-ATAC-seq 2018 and FISSEQ produced byte-identical pages after
their local build and test wrappers were removed. Their separate checks were small
enough to keep as `selftest()` functions. The central command removes repetitive entry
points without relocating source choices into shared code. SMART-seq2 and sci-ATAC-seq
2018 subsequently moved their source-specific oligo lists, citations and read layout
into their own modules. They may still reuse predecessor molecular constructors; their
published pages no longer call a predecessor's page renderer.

**Review when.** Apply the layout to a protocol with a much larger page or source-test
suite before making it the default. Split a module when separate code makes its evidence
boundary easier to review; minimizing file count is not itself a chemistry invariant.

## 2026-10-15 — Put source discrepancies on every affected public protocol page

**Decision.** If two authoritative publications disagree, or a source appears
internally wrong about a sequence, length, reagent name or reaction, each affected
public protocol page states the exact source locations and conflicting claims, the
value chosen for the model, the evidence for that choice and its remaining
uncertainty. A detailed research note may hold the full comparison, but a link to
that note alone does not make the discrepancy visible on the schematic. A derived
protocol that inherits the disputed component needs its own page notice.

**Why.** A reader may use the schematic without opening maintainer notes. A quiet
correction would make the diagram look more certain than the source warrants. The
notice lets the reader distinguish a verified correction, an alternate published
variant and an unresolved conflict. Errors introduced by this repository are fixed
and described as model corrections, never attributed to the publication.

**Review when.** Recheck the notice when a correction, original supplement or vendor
revision becomes available; update every affected page and the local model together.

## 2026-10-15 — Attach page emphasis to graph nodes

**Decision.** A protocol may attach typed `Highlight`, `DetailView`, and
`SourceDiscrepancy` annotations to existing step or material IDs. A detail view uses the
target material's drawing state by default; an explicit view may use a state derived
from the same molecular constructor, such as FISSEQ's circular topology or its
sequencing-primer scene. Four simple page slots place annotations without making
reaction order or chemical truth depend on display order. The shared renderer escapes
annotation text and labels every rendered item with its graph target. Page assembly fails
if any declared annotation slot was omitted or repeated.

The graph remains the authority for inputs, outputs, reactions, terminal fates and
sequencing plans. Validation warnings come from the chemistry checker and cannot be
suppressed or promoted by an annotation. A source discrepancy records both located
claims, the chosen model, its supporting evidence and remaining uncertainty; the
public page renders those fields together.

**Why.** The sci-ATAC-seq 2018 Tn5 branch, SMART-seq2 TSO disagreement and FISSEQ
circle/rolony views need different emphasis, but each already has a stable graph target.
Keeping the emphasis beside that target removes protocol-name decisions from page
renderers while preserving specialist molecular drawings. A page can increasingly be
assembled from graph data and a short list of selected views. `OligoSection` and
`WorkflowSection` are presentation declarations: the former receives protocol-owned
segment constructors and the latter renders the existing graph. They add no chemical
claim or source authority.

The three pilot readout sections now use graph-linked view declarations as well:
`PrimerBindingSection` for an exact run, `ReadLayoutSection` for source-declared cycle
recipes beside exact primer sites, and `OffsetPrimerSection` for FISSEQ's five primers
on repeated rolony units. The last requires an explicit sequencing-plan limitation.
These views check that they target a declared plan, but only molecular states and the
independent chemistry checker can establish a reaction or binding claim.

**Review when.** Try a dense protocol such as Small-seq before adding more annotation
kinds or page slots. If detailed figures require a second hand-entered molecular state,
improve the graph's state model first. Do not treat annotations as evidence or a
replacement for independently reviewed source coverage.

## 2026-10-15 — Check narrow chemistry claims from existing declarations

**Decision.** Keep protocol declarations as ordinary Python objects. Analyze each
`ProtocolGraph` in a separate, small checker (`lib/chemistry_validation.py`); do not
introduce a protocol DSL or treat a publishable graph as a biochemical proof. A material
may supply `chemical_state` when its `state` is a rendering view that has lost molecular
identity, as with FISSEQ's circle rows. The chemical state must come from the same
molecular constructor as the drawing, not from a second hand-typed sequence.

The checker reports *scoped claims*, each established, conditional, unresolved, or
contradictory. Its initial kernel rechecks paired columns and exact PCR primer sites,
preservation of strand order on circle closure and non-sequence-changing operations,
rolling-circle repeat derivation, restriction-end fill and contact-junction sequence,
and a declared hairpin stem. An IUPAC or placeholder pairing is conditional because
existence of a compatible assignment is weaker than a known exact duplex. A checked
stem does not establish that it is the preferred thermodynamic fold. A checked primer
site does not establish the resulting PCR product. Unmodelled input-to-output
transitions are explicitly unresolved; the checker must not upgrade them because a
step has the right roles or a page draws a plausible product.

Source fidelity remains a separate obligation: this checker does not verify that a
paper stated a sequence, modification, reaction order, concentration, yield, or
temperature. The source-text oligo guard and human reading continue to cover that
boundary. Chemistry claims are currently an audit, not a publication gate; coverage
is too sparse for a green summary to imply that a protocol has been proved.

An explicitly unmodelled transition may use `StepKind.MAGIC` with a required reason
explaining why the gap is needed. Its input and output still belong to the provenance
graph and source-reviewed scope. The chemistry checker leaves its transition unresolved
and emits a warning; the audit and graph page render every such warning. Warnings do
not establish that the proposed product follows from its inputs.

The displayed `StepKind` is an operation label and a graph-shape contract. A transition
claims a molecular check only by declaring one or more `MolecularRelation` values on
the step, backed by molecular input/output states. The initial relation set covers
circle closure, restriction-end fill, rolling-circle copying, derived contact junctions,
hairpin stem pairing, and strand-order preservation. PCR geometry, drawn pairing,
restriction-cut coordinates, and rendered/chemical-state agreement are already typed
declarations and are checked directly. A missing state leaves its declared relation
unresolved and visible as a warning; a contradiction is also shown. A successful local
relation never establishes an entire named reaction.

Two-ended adapter reactions may also declare `END_COMBINATIONS`. A loaded-end pool
supplies the source-backed adapter sequences; the step supplies the actual downstream
PCR primers and assigns every derived unordered end pair to a retained or waste output.
The checker tests exact 3′ primer footprints and complete, non-overlapping assignment.
This establishes eligibility for the specified two-primer library, not Tn5 insertion
mechanism, product frequency, or absolute inability to amplify with other primers.

**Why.** The ten local-ownership pilots contain 78 steps across 17 `StepKind`s but
mostly lack exact molecular inputs and outputs for a reaction-level proof. Five
boundary protocols add gap-sealing capture (smMIP), circle copying (FISSEQ), RNA
ligation/dephosphorylation (Small-seq), a symbolic fold (MALBAC), and a restriction
fill/contact join (in situ Hi-C). An independent audit of all 15 currently yields
25 established, 29 conditional, 114 unresolved, and zero contradictory *local claims*.
The two FISSEQ sequence transitions and the in situ Hi-C fill/junction can be checked
without inventing bases; the other examples expose missing state or source evidence.

**Review when.** Add a primitive only after a real protocol provides explicit inputs
and outputs whose relation it can check. Make a contradiction test for each new rule.
Broaden molecular-state coverage before considering this audit for publication gating.
Keep the four statuses per claim; do not compress them into a single protocol pass.

## 2026-10-15 — Protocol directories assemble their own reaction graphs

**Decision.** Each documented protocol owns the assembly of its typed reaction graph in
its own `tools/` module. The local graph states the ordered operations, inputs, outputs,
selection and discard decisions, evidence boundary, and final sequencing material.
`build_page.py` may remain a small executable entry point; a `PAGE_FILE` module may instead
use the central command. Protocol-specific citations, oligo displays, read layouts and
evidence notices belong with the protocol even when its molecular constructors are
inherited from a published predecessor.

Reuse molecular constructors and bounded chemistry operations when they express the
same operation across protocols: a cDNA amplification front end, an adapter ligation,
or a family member's published molecular state are examples. Such helpers accept the
caller's source-specific states and choices; they do not choose a named protocol's
reaction order or evidence status. A protocol that inherits another published method
may import that method's molecular constructors explicitly. Its own graph still shows
where that branch enters and what happens after it.

A bounded operation may contain several graph steps when the same sequence is supported
for every caller. Its name, inputs, outputs and scope should make those steps clear at
the call site. Reused molecular states belong in a model module with a stable import
path, separate from the page renderer. Source-specific definitions must not be moved
into `lib/` merely to make that import convenient.

**Why.** The ten-protocol ownership pilot found that a local `protocol_graph()` which
only calls a family graph builder hides the reaction sequence outside the protocol
directory. It also lets a change to one family builder silently alter several protocol
schematics. Local assembly makes the source boundary and reaction order reviewable
where the protocol's oligos and final construct are defined.

**Review when.** Revisit this boundary if the same reaction sequence is independently
verified in several protocols and duplication causes demonstrable drift. Share the
smallest invariant operation, retaining each protocol's visible graph assembly and
source decisions. Graph completeness and source-text coverage remain separate checks:
a publishable graph does not prove that every oligo has been transcribed exactly.

## 2026-10-14 — Local source text guards exact graph oligo inputs

**Decision.** `tools/selftest_source_oligos.py` checks exact, non-inferred oligo and
adapter input sequences already declared in a protocol's `ProtocolGraph` against
locally extracted source text. It runs across the catalogue by default and skips a
protocol when the defining source text is absent. A protocol module can declare
`SOURCE_OLIGO_EXCEPTIONS = {"material_id": "reason"}` for a known case where
verbatim comparison is inappropriate; stale and reasonless entries fail. The
test is separate from the normal offline site build because downloaded source
documents are never committed.
`--details` lists the unrepresented graph input IDs so the aggregate coverage count
remains actionable without treating an unavailable sequence as a transcription error.

**Why.** A model constant and a hard-coded test expectation can share the same
transcription mistake. Matching the model against a separately extracted source
document catches that error. The check requires every fixed base to agree, while
allowing a model N position to match a concrete source base. It covers only whole
input oligos with exact sequence states; unmodelled inputs, placeholders, source
modifications, and reaction order still require source review. A hit in extracted
text is a pointer to evidence, not proof that the intended source or oligo role
was interpreted correctly.

## 2026-10-14 — Every protocol receives an oligo source review

**Decision.** `catalogue/protocol_oligo_audit.tsv` tracks an explicit source-level
review for every entry in `catalogue/ours.tsv`, including pages that make no visible
missing-sequence claim and entries still at notes status. `tools/audit_protocol_oligos.py`
syncs this census and records the sources read, oligo roles checked, open roles, and
conclusion. A `reviewed` row means the named source scope was examined; it does not
assert that every oligo is published or that the molecular model is complete.
`limited` records an audit stopped by unavailable or ambiguous primary evidence.
Unresolved oligos asserted on public pages continue to use the separate, digest-bound
`oligo_source_audit.tsv` below.

**Why.** A page without a visible caveat may still omit an oligo or rely on a shallow
reading of the paper. The source-level census gives those omissions a durable work
queue and preserves the precise boundary between published bases, inference, and
unavailable material. Raw-read consensus can support an explicitly labelled inference
only when the dataset accession, processing, and cross-checks are recorded; it cannot
establish chemical modifications absent from a synthesis source.

## 2026-10-13 — Missing oligo claims have a source-review ledger

**Decision.** A schematic that says a primer, adapter, linker, barcode or sequencing
entry sequence is unavailable creates an auditable source-review obligation. The
checked-in `catalogue/oligo_source_audit.tsv` records one row per affected protocol,
the exact limitation signals currently rendered, the defining and supplementary sources
actually checked, and whether the sequence was recovered, confirmed absent, mixed across
components, or found not to be an oligo gap. `tools/audit_missing_oligos.py` derives the
candidate set from built public pages and fingerprints their current signals; completed
reviews whose rendered evidence changes become stale instead of silently remaining valid.

The ledger is a work queue and evidence boundary, not an authority for bases. Recovered
sequences and modifications still enter protocol notes and molecular constructors only
after verification in the original paper, supplement, detailed protocol or authoritative
vendor document. A secondary schematic can identify what to search for but cannot turn an
unverified sequence green. Sources that remain inaccessible or omit the oligo are named
precisely; “not found” without a list of checked materials is not a completed review.

**Why.** Prose search misses protocols produced by shared factories, while isolated page
repairs give no guarantee that every asserted gap was revisited. Deriving candidates from
the rendered output covers both direct and factory-generated pages. A dated, digest-bound
review preserves the human source-reading step that construction functions cannot infer.

**Review when.** Replace text-signal discovery with typed source-gap objects once all
protocol factories can emit them directly. Keep the ledger and stale-review behavior: a
typed model should make discovery more exact, not erase the record of what was checked.

## 2026-10-12 — Protocol completeness is a typed material-provenance graph

**Decision.** A complete protocol is represented as a directed acyclic hypergraph. Material
nodes include the biological template, molecular intermediates, structural oligos and
adapters, catalysts, selection reagents, sequencing libraries and material-bearing side
products. Typed step nodes consume one or more materials and create one or more new
materials. Construction is topological: a step may consume only existing nodes and may
produce only new nodes, making cycles and duplicate producers unrepresentable.

Every material branch ends explicitly in sequencing, discard, QC, stored intermediate or
an explained unresolved disposition. At least one terminal is a sequencing library. Every
sequencing library must have directed, molecule-contributing provenance from a declared
primary template and an exact sequencing plan or an explicit source limitation. All steps
and materials belong to one weakly connected protocol graph; connection through a reagent
does not substitute for directed template provenance.

The graph is typed rather than merely connected. PCR, ligation, digestion, tagmentation,
reverse transcription, selection, splitting, pooling and other operations impose their own
input and geometry contracts. Catalysts and conditions participate in a step without being
treated as incorporated sequence. Relevant material-bearing side products are accounted
for; the graph does not enumerate salts, hydrolysis products or every failed molecule.
Population operations may group a documented class of non-library products in one waste
node rather than pretending their proportions are known.

Circle closure has structured geometry for the same reason adapter ligation does: a
closed product cannot reveal whether the source used direct single- or double-strand
intramolecular ligation, splint-assisted ssDNA closure, overlap assembly or gap sealing
of a copied probe. Every circularization step records one of those mechanisms or an
explicit source limitation; the word “circularize” is not itself evidence of an end
configuration.

Process validity is distinct from completeness relative to an authority. Each reviewed
paper/manual operation or deliberate scope item is recorded as a source item and mapped to
a graph material/step, or explicitly excluded with a reason. This cannot discover a source
step a curator never recorded, but it makes the reviewed scope auditable and prevents a
connected graph from being mistaken for a complete transcription by itself.

`ProtocolGraph.require_publishable()` is the eventual source of documented/publication
eligibility. Linear `Workflow` objects acquire a graph in lock-step as a migration aid, but
generic construct-changing transitions are not publishable. Rendering is derived from the
native graph once a protocol has migrated; parallel sequencing outputs receive separate
provenance paths. The graph renderer lists every non-sequencing terminal and draws its
molecular state when present. When no state is available, the page says so explicitly.
The page author need not redeclare those branches as annotations. When the model validates
distinct molecular off-product classes and can construct their states, represent each as
a separate graph material. The sci-ATAC-seq 2018 T5/T5 and T7/T7 products use the same
loaded-end constructors as its retained T5/T7 product.

**Why.** The earlier model strongly checked bases and individual reaction geometry, but a
page could still omit a branch, disconnect an intermediate from its inputs, or declare
itself documented manually. A typed graph combines local chemistry guarantees with global
material accounting and makes completeness an executable property.

**Review when.** Refine terminal categories or population semantics when a real protocol
cannot be represented without false detail. Do not weaken directed provenance or source
coverage merely to make an incomplete legacy page pass.

## 2026-10-10 — Catalogue entries require a protocol name

**Decision.** A method enters the protocol catalogue only when its authors or vendor give
the workflow an explicit name. A paper that happens to describe an interesting sequencing
procedure is not enough: do not invent a catalogue name from its assay, paper title or a
descriptive phrase. Acronyms, versioned product names and explicitly titled modifications
of an established method all qualify. Closely related named protocols remain separate.

**Why.** A stable, source-owned name makes the catalogue unit reviewable and searchable,
and prevents arbitrary boundaries around every one-off library preparation found in a
methods section.

**Review when.** If an important unnamed chemistry cannot be represented through a named
parent protocol or a reusable concept page, reconsider whether concept pages need a more
prominent role; do not silently lower the catalogue bar.

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

Every defining DOI has a row in the snapshot, including an explicitly unresolved row.
Adding or removing such a protocol therefore requires running the manual refresh command;
the offline index loader rejects catalogue/snapshot drift without making network access a
normal build dependency.

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

Read-content tables follow the same rule: derive R1, I1, I2 and R2 content from the final
construct, primer direction, feature annotations and declared run lengths. Index reads are
auxiliary demultiplexing reads; identifier-only reads are extracted rather than aligned;
payload reads enter sequence analysis; undeclared roles are shown as not acquired. A small
structured override may be added for exceptional processing, but cycle layouts must not be
duplicated as page prose. The opt-in pilot spans SMART-seq2, Drop-seq, 10x 3′ v3,
ATAC-seq, SMART-seq3, TIME-seq, 5′ direct-capture Perturb-seq, 10x single-cell ATAC,
sci-RNA-seq3, inDrop v1, scNOMe-seq and CITE-seq. Together these exercise ordinary
payload, identifier-only, mixed identifier/payload, index-carried cell identity,
multimodal final libraries, custom short and absent reads before catalogue-wide use.
Processing overrides use a controlled vocabulary and may state handling when clipping,
unresolved proprietary sequence or protocol-specific analysis makes it impossible to
infer safely from molecular features alone; they never duplicate cycle layouts.

Standard sequencing-primer identity is likewise derived from the completed construct.
Protocol modules declare which read roles were acquired; the shared model matches canonical
TruSeq/Nextera primers to the final sequence and derives strand and direction. A protocol may
declare a primer locally only when the source specifies a genuinely custom primer. Missing or
ambiguous standard matches fail during construction rather than inviting a local relabelling.
Existing protocol modules should migrate to this rule as their chemistry is reviewed.

A subsequent 50-protocol expansion made three further distinctions explicit. A primer
site present on a construct is not proof that the corresponding read was acquired, so
protocol models declare acquired roles separately when the run omits an otherwise valid
primer. A schematic placeholder with a source-known molecular length records
`length_bp`, preventing its drawing width from becoming a false cycle span. Finally,
technical-only index reads are not called demultiplexing reads unless a structured sample
index feature is actually present. These rules are enforced in the shared renderer and
consumed by protocol modules; pages do not restate their cycle maps in prose.

An insert or cDNA placeholder with no molecular length never expands to fill a declared
run length. Its row starts at the computed cycle, calls the extent molecule-dependent, and
warns that a short molecule can read through into the opposite library end. This is the
same geometry in Atrandi, 10x, Parse and other variable-insert libraries; protocol pages do
not independently guess how many cycles contain payload.

The read-content section is now a publication invariant rather than an opt-in feature.
The preferred tier is computed from exact primer geometry; the second tier is a structured
source-declared layout for proprietary or otherwise unavailable primer sequences. At the
site-build boundary, a legacy page that supplies neither receives one explicit unresolved
``Sequencing output`` row. That last tier is deliberately uninformative: it prevents silent
omission while never guessing read names, direction, cycle boundaries or analysis handling.
Source-dependent stand-in pages receive the same row. As legacy models gain authoritative
data, their structured declaration replaces the fallback automatically.

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

Linear placement is not used for an assembly whose meaning depends on a loop, fold or
branched hybrid. Such assemblies use a topology constructor: callers declare the intended
complementary spans and molecular connections, the constructor rejects incompatible
bases or impossible end order, and the renderer lays out that validated graph. The layout
is schematic rather than a thermodynamic structure prediction. This rule is independent
of polymer type and therefore applies to DNA, RNA and mixed hybrids; a future external
folding engine may supply coordinates, but not replace the declared-pairing constraints
that establish what the protocol intends.

The nucleotide-resolved fold is a planar graph: a candidate with crossing backbone,
pairing or junction edges is invalid. Layout hints describe biology that pairing alone
cannot choose (for example, keep a target RNA straight), while unpaired runs between
bound spans are flexible joints by default. Among satisfying layouts, the renderer favours
short routes, fewer bends and directions close to straight. Layout and tie-breaking are
deterministic; any future stochastic optimiser must take a fixed, checked-in seed.

Sticky-end adapter ligation is drawn as two states: the ligation-ready dA-tailed molecule,
then the sealed product. In the product only, subtle semantic background spans distinguish
the insert from ligated components. Background identity follows the reagent identity, so
two copies of the same adapter at opposite ends share a colour; it is not assigned by
left/right position or inferred from sequence punctuation.

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

### Preprint names that change at publication

When the same chemistry is renamed between its preprint and journal article, it keeps one
protocol page under the published name. The earlier name is an explicit catalogue alias,
so it remains searchable and its paper remains linked, but it must not create a duplicate
schematic. This applies only when the sources establish identity of the workflow; related
or versioned methods still receive separate pages. Spatial-CUT&Tag, introduced as
hsrChST-seq, is the first recorded instance.

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

An exact oligo list beside a schematic does not make a prose reaction state molecularly
complete. When an authoritative source gives the relevant oligo, a construct-changing
step must consume named `Segment`/`Scene` structures derived from it and show the annealed
region, exposed end and resulting junction. Plain `Row` prose is reserved for sample-level
operations that do not change a construct, or for a structured source-limitation state.
Hybridization-guided barcode adapters use `HybridizationLigationRound`, which reconstructs
the splint from the incoming and barcode arms and rejects incompatible sequence geometry.

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

### Adapter-end geometry

Every adapter-ligation step also shows the ligation-ready adapter end before sealing.  If
the authoritative source gives the overhang, the drawing includes its polarity, bases and
the substrate end it pairs with; blunt ends are stated and drawn as blunt.  Splinted
ligations show both ligating oligos annealed side-by-side on the splint, including the
unpaired adapter tails.  The sealed product remains a separate state with the covalent
junction marked at the exact segment boundary.

If the source identifies an adapter-ligation reaction but does not disclose the end
geometry, the page says that the adapter overhang is not reported.  It must not replace
the missing bases with a generic sticky end or infer them merely from the final library.
This availability is structured construction data, not an optional prose caveat, so the
shared page path can require and render one of the two states: known geometry or an
explicit source limitation.

**Why.** A ligation product alone hides the chemical compatibility that made ligation
possible.  Showing the reagent end distinguishes blunt, cohesive and splinted reactions;
an explicit unavailable state prevents a sparse drawing from looking like evidence that
the end was blunt.

### PCR primer-binding geometry

Every PCR transition shows the template before amplification, both primer-binding sites,
primer-extension direction and the resulting amplicon boundaries.  Protocol builders use
`Workflow.pcr` with structured `PCRGeometry`; they do not jump directly from an unrelated
input to an amplified product.  Exact sites are located from the primer and template
sequences.  If a source gives a primer pool and only names its target region, the drawing
may show that region but must explicitly say that exact coordinates are not reported.  If
a primer sequence or target is unavailable, that absence is a structured state rather
than a silent omission.

**Why.** PCR defines which part of a molecule survives and which adapter tails are added.
A finished amplicon by itself cannot distinguish those choices or reveal an incorrectly
oriented primer.

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
