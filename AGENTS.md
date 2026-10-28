# AGENTS.md

Instructions for coding agents working in this repository.

Base-by-base documentation of NGS library chemistries (what the DNA construct looks like at
every protocol step), in the style of Teichlab's scg_lib_structs. `README.md` is the
short reader-facing introduction. `MAINTAINERS.md` holds the detailed build and handover
guide; this file is the condensed working guide for coding agents of any kind.
`CLAUDE.md` only imports it, so edit this file, not that one.

Project policy and its rationale live in `DESIGN.md`. Read it before changing site
presentation, protocol boundaries, evidence display, or model/test responsibilities.

## Commands

Pure Python 3 standard library — **no third-party dependencies, don't add any** (that is
why `lib/mdrender.py` is a home-grown Markdown renderer).

Approved browser-only exception (2026-10-21): the Anneal primers tab may load pinned
Pyodide 314.0.7 on demand to run the existing Python model in a Web Worker. Keep the
Python/site build offline and standard-library-only; add no Python packages, npm
packages or Rust toolchain. Runtime binaries are not checked in. See
[plans/annealing.md](plans/annealing.md).

```sh
# self-tests: shared primitives once, then only protocol suites affected by an edit
for t in */tools/selftest.py */*/tools/selftest.py; do python3 "$t"; done
python3 tools/protocol.py test fisseq__10.1126+science.1250212  # wrapper-free protocol
python3 small-seq__10.1038+nbt.3701/tools/selftest.py   # a single protocol's suite
python3 tools/selftest.py          # Markdown renderer, computed facts, docs build
python3 tools/selftest_sources.py  # source fetcher: DOI routing, download validation, fetchers (offline)
python3 tools/selftest_source_oligos.py [protocol-dir ...]  # local source-text transcription guard
python3 tools/check_chemistry.py --details  # scoped molecular claims across ten pilots + five boundary examples
python3 tools/lint_pairing.py [page.html ...]   # drawn duplex columns must base-pair

# build
for b in */tools/build_page.py; do python3 "$b"; done   # protocol diagram pages
python3 tools/protocol.py build fisseq__10.1126+science.1250212  # wrapper-free page
python3 build_docs.py              # every .md -> sibling .html + notes.html (only changed)
python3 build_docs.py --check      # exit 1 if any note page is stale; --force, --list
python3 build_index.py             # index.html (runs the self-tests to get check counts)
python3 build_site.py              # all of the above + _site/ (the website), link-checked
python3 -m http.server -d _site    # preview; CI deploys it (.github/workflows/pages.yml)
```

Citation counts are a checked-in OpenAlex snapshot, not a network build dependency. Run
`python3 tools/fetch_citations.py` whenever a protocol with a defining DOI is added or
removed (and whenever a deliberate refresh is requested; optional `OPENALEX_API_KEY`).
The index loader requires one snapshot row per defining DOI, including unresolved rows,
while normal builds remain offline.

There is no test framework: a selftest is a script using `checks.Check` (`check(label, got,
want)`), prints PASS/FAIL/SKIP and exits 1 on failure. To run one check, run its suite.

## Reading sources: the source tools

All articles, supplements and text conversions live in
`/husky/henriksson/NGS101_papers/sources/<slug>/`. This is the default in
`lib/storage.py`; `CHEM_DATA` explicitly overrides the storage root. Read and write
there directly. Do not recreate local copies or compatibility links in the checkout.
Older scattered references were moved under that root's `legacy/`, preserving
their former repository-relative paths. Use `storage.data_root()`,
`storage.sources_root()` and `storage.legacy_reference()` in Python readers.
For the shell examples below, set `CHEM_DATA` explicitly:

```sh
export CHEM_DATA=/husky/henriksson/NGS101_papers
```

Three standard-library CLIs in `tools/` for turning a paper into reference notes. Every
one has a full `--help`; this is when to use which. Their output is a **pointer into a
document, never a fact** — verify against the source text before anything goes into a
note or a protocol module, and mark it 🟢 only once you have read it there.

| Tool | Does | Typical call |
|---|---|---|
| `tools/get_sources.py` | Fetches everything legitimately reachable for a protocol in the catalogue (or a DOI): the upstream scg_lib_structs page (inspiration only, never a source — see below), PMC full text, the **PMC Cloud** bucket (JATS XML, PDF and supplements — also of author manuscripts whose PMC download page is gated), bioRxiv PDF + supplements, Nature-family supplements from Springer; links preprints and published versions both ways; then runs every special-case **fetcher** (below). Into `$CHEM_DATA/sources/<slug>/` with a `MANIFEST.tsv`. Every download is validated by type (a PDF must be a PDF, an XLSX a whole zip, an article page not a challenge page); a failed or invalid one becomes a `(manual)` row with its URL and reason. Skips what is on disk and valid. | `python3 tools/get_sources.py "sci-CAR-seq"` · `--slugs a,b` · `--doi 10.1101/829960` · `--todo` · `--all-ours` · `--revalidate` · `--list` · `--fetchers` |
| `tools/doctext.py` | PDF / DOCX / XLSX / HTML → a plain-text twin `FILE.txt` beside the original (XLSX as tab-separated rows under `## sheet <name>`). EPS figures are converted when Ghostscript is installed; unavailable conversion tools skip that file. `get_sources.py` runs it on everything it fetches. **Read and cite the `.txt`**: line numbers from the other tools point into it. Inspect original figures when layout or graphical marks matter. | `python3 tools/doctext.py DIR` · `FILE --stdout` |
| `tools/scrape_primers.py` | Lists every oligo-like sequence: the bases (PDF line-wrap rejoined; shorthand like `T30VN` and placeholders like `[8-bp barcode]` kept inside the oligo), the oligo **as written with its modifications** (`/5Phos/`, `rApp`, `-ddC`, `rGrG+G`, `*`, RNA `rN` or U), the parsed modifications, a nearby name, the `lib/` sequences it contains on either strand, other lines with the same oligo, and context. Tuned for **recall**: expect false hits. A family of table rows differing only in an index/barcode window is collapsed into one hit; files are ordered by their most oligo-like hit; `--max-hits` cuts the lowest scorers. `--find` treats N, U and IUPAC letters as base classes and sees through modifications and placeholders. | `python3 tools/scrape_primers.py $CHEM_DATA/sources/<slug>/` · `--tsv` · `--all` · `--known-only` · `--find SEQ` (either strand; `--strict` for letter-equal) |

**Special sources are fetchers, not hand downloads.** A source the generic routes miss
(a vendor user guide, a GitHub repository, protocols.io, a publisher API) gets either a
row in `catalogue/extra_sources.tsv` (slug, source, what, why — `source` is a URL or
`handler:arg`, e.g. `scg:data/BD/*`, `protocolsio:<doi>`, `doi:<doi>`, `pmccloud:<PMCID>`
or `pdf:<URL>` for an identified PDF whose URL lacks a PDF suffix)
or, for a new kind of source, a module in `tools/fetchers/` (how-to in
`tools/fetchers/__init__.py`; `--fetchers` lists them). For 10x Genomics and BD the
scg_lib_structs GitHub repo (`scg:` handler) mirrors the vendor guides, oligo tables and
barcode whitelists; it is fetched automatically for every protocol upstream has drawn.
Offline checks for all of this: `python3 tools/selftest_sources.py`.

A cached filename is not source identity. The fetcher refuses to reuse or overwrite
a manifest entry assigned to another URL; fix the filename mapping or declare a
distinct original, rather than bypassing that check. In particular, GEO accession
and record scope belong in filenames for its shared `acc.cgi` endpoint. Check the
manifest URL and original accession before counting a cached document as read.

The protocols.io fetcher follows attachment links in both the PDF and the public
metadata for the PDF's printed protocol integer ID. Legacy XLS order sheets may be
linked only in that metadata. A PDF with no attachment annotation is therefore not
evidence that the named workbook is inaccessible; do not substitute a newer protocol
revision without checking its identity.

What is still out of reach: Science / Cell / most publisher PDFs and supplements behind
Cloudflare or a paywall — get those by hand (the manifest has the URL) and drop them in
the same directory; `--revalidate` then checks them like any download. Do not work around a
challenge (no browser-fingerprint spoofing, headless challenge solving or proxy rotation).
`tools/fetchers/base.py` spaces requests per host across all fetcher processes (8 s for
publishers) and stops requesting a host for six hours once it serves a challenge or HTTP 429;
the cool-down markers are under `$CHEM_DATA/sources/.cache/ratelimit/`. The scraper
finds oligos, not reaction order: the methods still have to be **read**.

**Source revisions are part of identity.** Preserve a cited revision in retrieval
URLs and cache names. Check the returned original's printed revision and exported
identity before treating it as that source. A current export cannot stand in for
an earlier cited revision; matching selected components does not prove complete
revision equivalence. Rejected responses stay outside ordinary source reading.

**Published version first.** If a paper is published, use that article and its
published supplements; normally disregard its preprint and preprint supplements.
Do not merge versions or treat an ordinary revision change as a published error
or unresolved source conflict. A preprint may supply a specific necessary detail
omitted from the publication only with an explicit reason, verified applicability
and a provenance/uncertainty label; it never overrides a published value. Keep
historical audits as history, not current authority. For an excluded superseded
preprint, independent absence-audit readers check its identity and supersession;
they need not re-read its chemistry. See DESIGN.md.

**Capped breadth-first curation.** Work from the selected round ledger in
`catalogue/curation_rounds/`. Inspect defining originals and chemistry supplements
across the batch before following links. Investigation is limited to 10 minutes
and three originals per protocol; leave justified edits consistent and verified,
then move on. Record inspected sections and actual deferred targets. Do not
recursively chase references, build fetchers or expand into an architecture
migration. Later follow-up rounds stop at directly linked originals, depth one.
A completed capped pass is not a complete protocol review. The missing-information
gate below remains unchanged.

**Missing-information audit gate.** Do not call a protocol's missing sequence,
reagent, reaction detail, or read layout confirmed absent merely because a source
bundle or scraper returned no hit. Keep a checklist of the defining article,
supplements, detailed protocols, cited predecessor methods, vendor documents, and
other chemistry-relevant material linked by those sources. Record each source's
locator, availability, exact sections inspected, and whether linked information was
followed or why it was inaccessible or irrelevant. Two independent readers must
check the original sources and linked material, including any source judged
irrelevant, before a missing-information claim is accepted; an unavailable
relevant source leaves the claim blocked, not confirmed.
Record the two passes and their source checklist in `catalogue/source_reviews/`.
The existing oligo audit ledgers are historical work queues, not evidence of this
two-pass gate by themselves. Show the source checklist and linked-information
disposition on a public page when that page asserts confirmed missing information.

Each checklist row must identify an original document, repository or file. Keep a
desired component specification or an unassigned study-specific value as an open
question; do not invent an unavailable disclosure document for it. A genuinely cited
source with unresolved identity still blocks the associated absence claim. This
distinction does not waive source checks or establish that information is absent.

## Writing up a protocol (notes first)

1. `python3 tools/get_sources.py "<protocol>"` — fetch; check the manifest for `(manual)`.
2. **Always run** `python3 tools/scrape_primers.py $CHEM_DATA/sources/<slug>/`
   when reading an article for a protocol, including its extracted supplements. Review
   the hits and modifications before deciding that an oligo sequence is unavailable.
   Use `--find` for a suspected sequence and `--all` to inspect collapsed primer
   families. The scraper is an aid to reading, never page content or evidence by
   itself: follow its source-line pointers and verify each relevant hit in the
   original article or supplement before using that source fact in notes or a page.
   If the article points to a sequence table but its `.txt` contains only a caption
   or no bases, inspect the rendered original. Image-based oligo tables can produce
   zero scraper hits; an empty extraction is not an inspected table.
3. Read the methods in the `.txt` files. **`upstream_*.html.txt` is not a source.** The
   scg_lib_structs drawn page is inspiration — useful for orienting yourself in a
   chemistry, never evidence: do not cite it, do not take a sequence from it, and do not
   treat agreement with it as verification. Every base and every reaction comes from the
   paper, its supplements, or a vendor document. **A filename prefix is not provenance.**
   Publisher supplements and vendor guides *mirrored* in that repository are ordinary
   originals, but `scg_docs_*.md` is the drawn page's own markdown and `scg_adapters_*.fa`
   is scg's self-authored FASTA — upstream content under an `scg_` name — while a genuine
   vendor PDF can land as `upstream_*` when the page links it directly (the 10x CG00026
   guide is on disk under both prefixes with one md5). The `what` column of `MANIFEST.tsv`
   is the authority. `tools/source_oligo_check.py` excludes `upstream_*`, `scg_docs_*` and
   `scg_adapters_*`; a vendor original that landed under `upstream_*` is therefore invisible
   to that guard, and should be re-fetched through the `scg:` handler so it lands as `scg_*`.
4. Write `<slug>/01_<name>.md` (`<slug>` is the `slug` column of
   `catalogue/scg_lib_structs.tsv`): citation, sources read, what is new, oligos
   **verbatim with modifications**, how they interlock (compute with `lib/`, do not
   eyeball), step by step, final library, read layout, open questions. Example:
   `astar-seq__10.1101+829960/01_astar-seq.md`.
   Curate method-development or component-adaptation links alongside the notes in
   `catalogue/protocol_lineage.json` (see MAINTAINERS.md, "Protocol relationships").
   Multiple parents are allowed. Record inspected evidence, scope and limits, and
   distinguish documented ancestry from design-based inference. Shared primers
   suggest links to investigate; computed similarity never establishes ancestry
   or authorizes inheriting chemistry. Leave uncurated origins unspecified.
5. Add the directory to `catalogue/ours.tsv` with status `notes`, `section` `published`
   (or `wip` for our own unpublished work; the catalogue self-test pins the wip set, so
   update its `WIP` too) and `modality` `DNA` / `RNA` / `multi`. Re-run
   `python3 catalogue/tools/fetch_scg_lib_structs.py` (cheap: reuses resolved papers),
   `python3 build_docs.py`, and the self-tests. The protocol is listed by name on the public
   front page, but its notes remain maintainer-facing. A `tools/` module, self-test and
   finished schematic page (status `draft` → `documented`) come later; only then does its
   title link to a public page — see MAINTAINERS.md "Adding a protocol".

## Architecture

- **`lib/`** — shared, protocol-agnostic model. `chemdraw.py` (Segment / Construct / Scene,
  complement/revcomp, Tm), canonical sequences in `illumina.py` and `nextera.py`, chemistry
  helpers (`rt.py`, `padlock.py`, `crispr.py`, `plasmid.py`), `seqprimers.py` (locates
  Read 1/2 + index primers on a library by computation), `mdfacts.py` + `mdrender.py`
  (Markdown with computed facts), `page.py` (shared HTML/CSS), `checks.py` (harness +
  `run_common`, retained for legacy suites; new protocol suites do not repeat it).
  Directory placement does not make a helper neutral. The completed ownership
  migrations are recorded in `catalogue/chemistry_decoupling.tsv`; do not reintroduce
  retired source-selecting recipes. A readout label check is not reaction validation;
  distinguish the chemical product from the base reported by sequencing.
- **`protocol_families/`** — remaining shared Chromium partition helper. Legacy
  family recipes were migrated to their protocol owners and retired (see
  `catalogue/chemistry_decoupling.md`). A family name or matching endpoint
  does not establish shared chemistry. Shared operations must receive the molecular
  inputs, reagents, geometry and outputs explicitly; they must not supply biological
  defaults or source conclusions. Protocol modules assemble their own graphs.
  Sharing an operation does not establish equal substrates, ends, conditions or
  products. Keep those choices source-qualified locally; derive comparisons from
  declarations without using apparent similarity to fill gaps in either protocol.
  A pure function with explicit parameters can still encode a biological assumption.
  Keep expected products, off-products and physical versus computational fates local;
  shared checks establish only their stated molecular relation. Unknown values are
  not evidence that two protocols agree.
- **`<protocol>__<doi with / → +>/`** — one directory per protocol, named for its defining
  paper. `catalogue/ours.tsv` is the single record mapping protocol → directory, and the
  catalogue self-test rebuilds directory names from it (so renames must go through it).
  Each contains `0*.md` notes, `ref/` (with `MANIFEST.md`), and `tools/`. Most have
  `<chem>.py`, `build_page.py`, and `selftest.py`; a wrapper-free protocol instead has
  one module declaring `PAGE_FILE`, `protocol_graph()`, `render_page()`, and `selftest()`.
  `tools/protocol.py` builds or tests it, and `build_site.py` discovers it automatically.
  Keep that protocol's citations, displayed oligos, read layout and evidence notices in
  its module. Import a predecessor molecular constructor only when local evidence
  establishes inheritance of that specific component and version; record that scope
  beside the import/use. Similar names or matching sequences alone do not establish it.
  Status `notes` directories hold reference notes only. Scripts bootstrap
  `sys.path` with `HERE` and `HERE.parents[1] / "lib"` (see MAINTAINERS.md "Adding a protocol").
  Construct-building pages use status `construct` and an explicit
  `ProtocolGraph(..., required_output=Disposition.CONSTRUCT)` with a constructed-DNA
  terminal. Do not invent a sequencing output for them. The graph renderer and site
  build check their construct endpoint separately from NGS readiness.
- **`catalogue/`** — scraped worklist of protocols from scg_lib_structs
  (`tools/fetch_scg_lib_structs.py`, network, cached). `properties.tsv` drives the public
  chemistry finder: every `ours.tsv` directory must choose a profile and record overrides;
  `tools/properties.py` rejects missing facets and vocabulary drift during the normal build.
- **Website** — `build_site.py` builds every page but publishes only the front page,
  finished schematic pages, their main research note (linked quietly from the schematic),
  and files those pages deliberately link to. Other reference notes, evidence logs,
  catalogue internals and work in progress stay out of public navigation
  and search; never publish `_data/`, `pdf/`, caches or `ref/` data files.
  Broken relative links fail the build. A diagram page whose build needs a missing source
  is omitted (a stand-in page explains why). `index.html` is a searchable list of the
  published protocol names from `catalogue/ours.tsv`; names link only to successfully built
  schematics, never to working notes.
- **`ref/concepts/`** — recurring chemistry (Tn5, template switching, small-RNA ligation,
  RT, padlock) written once; most new protocols are a new front end on one of these.

## Rules that are easy to break

- **Never hand-type sequence columns.** Any multi-strand panel is a `chemdraw.Scene`:
  give strands 5'→3' and place them with `sc.anneal(..., pair=(segA, segB))`; the Scene
  computes columns and raises if bases don't pair. Anchor marks/arrows to segments.
  `Row(indent=...)` and `markers([(col, ...)])` are legacy — don't add new uses.
- **A reaction is a `chemdraw.Workflow`, not several unrelated panels.** Start from one
  `MolecularState` and call `react()` for each operation; it automatically uses the prior
  output as the next input. Render it once with `workflow_panel()`. Reserve `panel()` for
  a molecular state such as a final library or primer-binding view.
- **Adapter ends are part of the ligation-ready state.** Use
  `AdapterLigationGeometry` with `Workflow.ligate_adapter()`: draw source-known cohesive,
  blunt, splinted or single-strand end geometry (including polarity and bases when known),
  or construct an explicit `unavailable(...)` state when the authoritative source does
  not disclose it. The sealed product is a separate state with its exact junction marked.
- **SVG geometry comes from model semantics, not punctuation.** `Scene` marks strands,
  binding spans and process arrows as typed row visuals; shared rendering supplies arrow
  boxes, proportional-font comments and segment hover name/length/Tm. Do not draw `^^^^` or
  `<--------` in new SVG-specific code or hand-write tooltip values in page builders.
- **Import canonical sequences** from `illumina` / `nextera` rather than retyping them.
- **Facts in Markdown come from the model:** `{{= module.EXPR }}` and
  ```` {{oligo: module.fn()}} ```` are expanded by `build_docs.py`; every lib and protocol
  module is in scope, and an unresolvable fact is a build error. Don't hard-code numbers
  in prose that the model can compute.
- **Sequencing primers are declared** by reference (`seqprimers`), and `sp.section()`
  derives where they bind; a declared primer without a site fails the build unless the
  protocol says why.
- **Identifier meaning is structured data.** Attach a `chemdraw.feature(...)` to every
  barcode, sample index and UMI segment. Do not infer meaning from `cbc`/`umi` colours,
  placeholder letters or prose names. Preserve feature identity through transforms and
  derive read-cycle spans with `seqprimers.feature_spans()` from declared run lengths.
- **Put invariants in construction functions, not tests.** Derive complements and adapter
  assemblies, use `Scene` to reject invalid pairing, and let sequencing-page generation
  fail when primers do not land. Tests are for source transcription or behavior that the
  model cannot enforce; do not restate constants or prose as checks.
- **Drawing convention:** the bottom strand is the plain complement (not revcomp),
  written left to right as `3'-…-5'`. Placeholders (barcodes, linkers, inserts) are never
  complemented as bases — that is a property of the `Segment`; their complements render
  lowercase.
- **Uncertainty is a constructor property, not prose:** uncertain helpers must create
  inferred segments by default, so callers cannot accidentally render them as established.
- **Evidence marking** in notes: 🟢 verbatim from source · 🟡 derived/inferred · 🔴 not
  published. On pages, inferred regions need all three: a preamble caveat, an
  `INFERRED — …` step caption, and `<inf>` around the bases.
- **Show published discrepancies on every affected protocol page.** When a paper,
  supplement or vendor original contains an apparent error or conflicts with another
  authoritative source, identify both source locations and values on that protocol's
  public schematic page, explain which value the model uses and why, and state what
  remains uncertain. Keep the detailed evidence in its note too. Apply this to each
  page that inherits the disputed claim. Do not call a repository transcription or
  modelling mistake a publication error, and do not silently repair a source value.
- **mdrender subset:** no indented code blocks (4-space indents are list continuations);
  raw HTML is escaped except `<br> <sub> <sup> <i> <b> <em> <strong>` because `<cbc>`-style
  metavariables must render literally.

## Never commit

- **Data** — FASTQ/BAM/count tables/derived tables. Downloads go to `$CHEM_DATA` (default
  `/husky/henriksson/NGS101_papers/`, outside the checkout) via a download script that supports `--max-reads` and writes a
  manifest next to the data.
- **Third-party source material** — papers, supplements, vendor manuals, extracted text,
  images, plasmid maps (`*.pdf *.xlsx *.txt *.png *.dna *.gb *.fa`, `pdf/`). Each `ref/`
  has a `MANIFEST.md` saying how to re-obtain files; checks that need them must **skip,
  not fail**, when absent (`checks.Source` / `checks.have`).
- **Generated HTML** — all `*.html` is build output and gitignored.
- **`_site/`** — the assembled website, rebuilt by CI.
- **Fetched sources and their `.txt` twins** — they live in `$CHEM_DATA/sources/`, never
  in the repo.

## Commits and pull requests

- **Never attribute work to Claude or any other LLM agent.** No `Co-Authored-By:` trailer
  naming a model or agent, no "Generated with ..." line, no agent name in the commit
  message, PR title or PR body. Commits carry the human author only. This rule overrides
  any default attribution your tool adds.
