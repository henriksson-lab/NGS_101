# Plan: sequence annotation tab

Status: implemented, 2026-10-21. The sections below retain the design rationale;
the implementation and verified index-export scope are recorded at the end.

## User experience

Add an **Annotate sequence** tab next to Index and About. A user pastes a sequence
or one FASTA record and receives an annotated sequence, an accessible match table,
and a bounded list of protocols sharing the matched components. Exact matching is
the initial scope, on both orientations. Protocol suggestions explain which shared
parts support them; they do not identify an experimental protocol from sequence alone.
The sequence stays in the browser, with no upload or persistence.

The front page searches protocol text and facets. The separate static annotation
page shares its navigation and keeps catalogue search intact.

## 1. Build a deduplicated component catalogue

Add `tools/build_sequence_index.py` and a small neutral extraction helper. Use the
existing graph discovery in `tools/selftest_source_oligos.py` as a reference, with
isolated protocol loading so local modules such as `chem.py` cannot leak between
protocols. Index only protocols whose public schematic built successfully; preserve
explicit diagnostics for skipped graphs and unsupported state types.

Initial high-value inputs:

- Full declared oligos and adapters from `MaterialKind.OLIGO` / `ADAPTER` inputs,
  their fixed named segments, and their representation as independent `Scene` strands.
- Forward/reverse primers from `PCRGeometry`, retaining exact versus unavailable
  geometry and any unverified use limitation.
- Sequencing primers from terminal `SequencingPlan`s. A compatible reference primer
  must remain labelled as a reference rather than a verified loaded reagent.
- Canonical components from `illumina.py` and `nextera.py`, with protocol associations
  only where the declarations actually contain or reference the component.

Start with these inputs rather than indiscriminately indexing every example insert
or intermediate. The catalogue inspection stage can establish which additional exact
`Segment` classes are useful. Do not infer component meaning from colour tags or names.
Use structured features and the role/location of the declaration where available.

A known fixed subsegment remains useful even when another part of its oligo is unknown.
Never join across a placeholder, unknown-length region, inferred segment or independent
strand. Placeholder letters that resemble DNA are not known bases. Do not generate
arbitrary k-mers as components. Full contiguous fixed oligos and their explicitly named
parts may both be indexed, with their containment relationship retained.

Default discovery threshold: 12 nt, configurable. Exclude low-complexity homopolymers,
arbitrary barcode/index examples and short incidental matches from ordinary results.
Inferred sequence instances can be retained as a separate, initially hidden class;
they must never silently acquire the evidence status of an exact source-qualified
occurrence of the same bases. This catalogue is an index of declarations, not a new
source-reading certification.

## 2. Deduplicate sequence identity, retain distinct occurrences

Use one canonical nucleotide record per sequence/reverse-complement pair, for example
`min(sequence, reverse_complement(sequence))`. Every occurrence retains its declared
orientation. A palindromic sequence is one identity, not two strand hits.

Separate records for:

- **Component:** stable ID, canonical bases, aliases, known parent/contained components.
- **Occurrence:** protocol, material/step/plan locator, strand and component orientation,
  local name, usage kind, source, inference/reference flags and available feature metadata.
- **Protocol:** ID, display title and public schematic URL.

Collapse repeated appearances of a component across steps within one protocol. Preserve
multiple meanings or source scopes as metadata instead of asserting equal chemistry.
Base equality does not imply equal modifications, molecular ends, experimental use or
identifier role. The same component can have both established and inferred occurrences.

Generate a compact JSON asset, a format version and a build report with component counts,
protocol coverage, excluded short/low-complexity parts and unsupported declarations.
No downloads or source-text parsing should occur during this build.

## 3. Match exact substrings in the browser

Create a pure matcher in `assets/sequence_annotation.js`, separate from the UI code,
and keep it worker-compatible. Begin with repeated `indexOf` on unique concrete strings
and their reverse complements; advance by one base to preserve overlapping matches.
Benchmark against the actual catalogue before choosing a more complex matcher.
If needed, add a seed lookup with full exact verification or an Aho–Corasick automaton.
Neither extension should introduce approximate matching or third-party dependencies.

Normalization rules:

- Accept plain sequence or one FASTA record; remove whitespace and uppercase.
- Reject multiple FASTA records and unexpected symbols with a useful message.
- Accept IUPAC ambiguity letters as barriers to concrete matches, not wildcards.
- Keep U distinct from T by default; any U→T conversion must be explicit and visible.
- Use zero-based half-open coordinates internally and display one-based inclusive
  nucleotide positions. Coordinates refer to the normalized sequence.
- Start with a 100,000-nt input limit. Cap stored/displayed hit counts explicitly;
  if a limit is reached, show truncation rather than imply a complete negative result.

Load the index on demand. Prefer a Web Worker if benchmarks show that realistic or
maximum-length inputs block the UI. Keep a simple main-thread path for small inputs
only if it meets the interaction budget. Benchmark common/repetitive inputs as well
as random sequences; exact matching alone does not guarantee bounded latency.

## 4. Render annotations and preserve overlap

Show the normalized submitted sequence with component spans, coordinates and direction.
Use a boundary sweep to render HTML spans; do not hand-compose a second sequence copy
for each hit. Clicking a region reveals all matching identities at that location.
Default to the longest informative matches in the overview, while retaining contained
and overlapping parts in the results table or expanded detail.

The accessible table lists position, orientation, names/aliases, matched length,
usage/evidence class and number of associated protocols. A matched component inside a
longer primer is a contained match, not independent extra evidence. Preserve both
records without showing a stack of redundant labels by default.

Results should say **exact sequence matches** and **protocols sharing matched parts**.
A substring hit does not establish annealing, biological origin, modification chemistry,
a preferred fold or protocol identity. Use structured identifier metadata where present;
use neutral “sequence component” labels where no role is declared.

## 5. Bound protocol suggestions and suppress ubiquitous components

Suggested configurable defaults:

- Show up to 8 protocol examples per component, with “8 of N protocols” and an
  expansion capped at 20 visible links.
- Show at most 20 protocols in the overall comparison list.
- Count each component once per protocol, regardless of repeated occurrences in the
  query or repeated graph states.
- Rank by informative components and matched-base coverage. Downweight a component by
  its protocol prevalence, such as an inverse-document-frequency weight.
- Collapse nested supports at the same query interval before scoring, so P5 inside a
  complete matched primer cannot inflate the score twice. Compute coverage by interval
  union rather than summing overlapping lengths.
- Use deterministic tie-breaking and show the supporting components beside each result.
- If only ubiquitous components such as P5 match, show their counts and a few examples
  under “common components”; do not present arbitrary protocols as strong suggestions.

Ranking is a navigation heuristic, not a confidence score or evidence of equal chemistry.
Define “common” from measured catalogue prevalence during the catalogue review, rather
than baking in a P5-specific exception.

## 6. Build and navigation integration

- Add the navigation link in `lib/page.py:site_header()` with the existing relative-path
  handling for root pages and protocol subdirectories.
- Generate `annotate.html` using shared page styling, and register it in `build_site.py`.
- Publish the component JSON, matcher/UI scripts and any worker explicitly. The current
  dependency walker follows HTML `href`/`src`; it will not discover an asset referenced
  only by JavaScript `fetch` or a worker constructor.
- Add the required asset types to the site allowlist if needed, and ensure these assets
  are included in the normal link/build checks. Use built schematic URLs, not source
  notes, caches or private provenance files.
- Keep generated runtime assets separate from authoritative protocol declarations.
  Avoid adding another hand-maintained primer table.

## 7. Meaningful verification

Use small synthetic fixtures, with no original-paper downloads needed:

- Extraction: repeated components, aliases, reverse orientations, palindromes,
  separate strands, unknown gaps, placeholders that resemble bases, inferred/reference
  occurrences and modification metadata.
- Matching: overlapping/nested/repeated hits, both orientations, normalization, invalid
  input, ambiguity barriers, U handling and explicit result truncation.
- Ranking: ubiquitous components do not dominate distinctive support; nested matches
  do not double-score; counts, ties and visible caps are deterministic.
- UI: keyboard-accessible annotations/table, escaped names/source strings, empty/no-hit
  states, index-load failures and long/repetitive input behavior.
- Integration: new page and assets are published, protocol links resolve, shared tests
  pass and the full static site/link build remains green.

Choose a lightweight way to run JavaScript matcher fixtures in CI using available
runtime tooling; do not add a package ecosystem solely for this feature. Confirm that
runtime before making the test harness part of the required build.

## Delivery order

1. **Catalogue first:** build the extractor and deduplicated index, review representative
   components and prevalence, publish a coverage report, and benchmark matching. This
   exposes declaration gaps before they become silent UI omissions.
2. **Usable tab:** normalization, exact matches, overview/table, bounded protocol links
   and normal static-site integration.
3. **Refinement only as needed:** worker/optimized matcher, overlap interaction,
   prevalence ranking and additional component classes based on measured behavior.

## Delivered implementation

- `lib/sequence_catalogue.py` extracts graph declarations; `tools/build_sequence_index.py`
  loads each protocol in a separate process. The normal site build supplies only public
  protocols whose schematics built, with their actual output URLs. JSON includes a
  per-protocol coverage report; missing graphs remain explicit.
- Canonical informative fixed-DNA constants are recognized inside accepted concrete runs.
  Their occurrence retains the parent locator, source and inference status. This asserts
  contained sequence identity, not use of a particular experimental primer.
- `assets/sequence_annotation.js` uses a cached Aho–Corasick matcher. Final measured
  100,000-base matching took about 0.4 seconds cold on random DNA and 0.7 seconds on
  catalogue-rich input (truncated at the hit cap), in Node on the development machine.
  Shorter typical queries use the cached matcher; a browser worker remains a possible
  refinement if browser measurements justify one. Limits are 100,000 input bases,
  10,000 retained hits, 20 suggested protocols and 1,000 index candidates, with visible
  truncation. Matches under 12 nt and low-complexity keys are excluded.
- `build_annotate.py`, the UI script and CSS provide the tab, paginated overview, overlap
  inspection, bounded protocol links, explicit inference opt-in and clipboard fallback.
  Inferred-only aliases and associations are excluded unless requested. Sequencing
  primers remain labelled references, rather than verified loaded reagents.
- Normal `build_site.py` publishes the page, JSON and browser assets. CI runs the Python
  extraction suite and dependency-free Node matcher, orientation and DOM integration tests.

## Illumina sample-sheet index export

`assets/illumina_samplesheet.js` contains verified instrument/workflow profiles, separate
from sequence extraction. Both selections are required. The UI lists resolved i7 in
`index` and i5 in `index2`, showing the transformation, individual copy buttons, a TSV
header/row and a CSV data-section snippet. A full data row requires i7; i5 alone remains
individually copyable. BaseSpace Run Planning provides setup inputs rather than a manual
CSV. Unknown combinations cannot default to another orientation.

Candidates require a structured `sample_index` feature and compatible matched template,
with at least 12 concrete flanking bases on both sides. Canonical contiguous P5/read-1
or P7/read-2 arms establish the Illumina-forward axis. Reverse query strands and raw
P7-primer orientation are normalized before applying the software profile. An ambiguous
index, unknown axis, insufficient anchoring or contradictory orientation remains
unresolved. A free barcode-length query is insufficient. UMIs and inline cell barcodes
are not exported as sequencer indexes. Selecting an i7 limits i5 choices to the same
protocol, without claiming that independently declared choices establish a valid pair.

Profiles cover original MiSeq, MiniSeq standard/Rapid, iSeq 100, NextSeq 500/550 and
1000/2000, NovaSeq 6000 v1.5 and X/X Plus, MiSeq i100 sequencing orders, and legacy HiSeq
manual bcl2fastq. Workflows distinguish manual bcl2fastq, standalone/DRAGEN BCL Convert,
BaseSpace Run Planning and BaseSpace BCL Convert v2. These specify index orientation,
not software compatibility for every instrument. Local Run Manager/IEM are not covered;
MiniSeq Rapid has a distinct LRM exception. Required `IndexOrientation,Forward` settings
are included in a `[Header]` section before `[BCLConvert_Data]` when applicable. Snippets
are not complete run configurations and must not replace other required settings.

Rules checked on 2026-10-21 against manufacturer documentation:

- [Illumina index orientation guidance](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000001800)
- [Connected Software i5 orientation table](https://help.connected.illumina.com/run-set-up/overview/index-orientation-guide/i5-index-orientation-table)
- [BCL Convert sample-sheet format](https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm)

The export transforms declared/observed nucleotide strings; it does not independently
verify the underlying protocol sources or establish index-pair compatibility.
