# NGS protocol chemistry

Base-by-base documentation of sequencing library chemistries — what the DNA construct
actually looks like at every step of a protocol, drawn in the idiom of the Teichmann lab's
[scg_lib_structs](https://teichlab.github.io/scg_lib_structs/) pages.

Project-level presentation and modelling choices are recorded in [DESIGN.md](DESIGN.md)
so they can be reviewed independently of the implementation.

Each protocol gets a set of reference notes and one generated HTML page. The sequences and
diagrams are **generated from a segment table, not hand-written**. Intended invariants are
enforced by the functions that construct and render the model; tests cover only properties
that cannot be made true by construction.

## Why it is built this way

Hand-aligning a dozen near-identical duplex ladders is exactly the kind of work that goes
wrong silently: a reverse complement written backwards, a primer arrow two columns off, an
annotation row that drifts the moment a segment length changes. Several real errors were
caught this way during the first protocol — a primer line double-counting an overlap
between two adjacent sequences, an adapter strand written 5'→3' where it should have been
3'→5', a junction base on the wrong strand.

So: define the construct once, derive every diagram from it, and make invalid structures
unrepresentable in the construction API.

**No hand-typed columns.** Any panel with more than one strand is a `chemdraw.Scene`: each
strand is given 5'→3' as ordered, and placed by naming which of its segments pairs with
which segment of a strand already drawn (`sc.anneal(..., pair=("oligo-dT", "polyA"))`).
The Scene computes every column, draws the antiparallel strand reversed, and **raises** if
any drawn column does not base-pair. Marks and arrows are anchored to segments, not columns.
`Row(indent=...)` and `markers([(col, ...)])` are legacy; don't add new ones.
`tools/lint_pairing.py` (run by `tools/selftest.py`) checks every generated page's labelled
duplexes as a backstop. All seven protocol pages are converted; nothing in `lib/` or any
`tools/build_page.py` positions sequence by a typed column any more.

**Facts in prose come from the model too.** A Markdown note writes `{{= crisprmip.PROBE_LEN }}`
or ```` {{oligo: crisprmip.probe_construct()}} ```` and `build_docs.py` substitutes the
computed value (`lib/mdfacts.py`). Every lib and protocol module is in scope. A fact that
stops resolving is a build error, so prose cannot quietly disagree with the chemistry it
describes -- which is how the notes came to state a 6-nt index that was really 8, and a
closed circle of 238 nt that was really 237.

**Sequencing primers are declared, not described.** Read 1 / Index 1 / Index 2 / Read 2 live
once in `lib/illumina.py` and `lib/nextera.py`; a protocol lists which it uses by reference
(`lib/seqprimers.py`) and `sp.section()` computes, from that page's own final library, which
strand each primer anneals to, what it covers and the first bases each read reports. A
declared primer with no site -- 3' end exact over >= 18 nt -- fails the build unless the
protocol declares why.

**One directory per protocol, named for its defining paper.**

```
<protocol>__<doi, with every "/" replaced by "+">
smart-seq2__10.1038+nmeth.2639
```

Protocol names collide and get reused; a DOI does not, and every preprint has one from the
day it is posted. `catalogue/ours.tsv` is the single record of which protocol lives in
which directory, and the catalogue self-test rebuilds every directory name from it, so the
name on disk and the paper it claims cannot drift apart. `catalogue/README.md` explains the
scheme; `catalogue/scg_lib_structs.tsv` is the worklist of what is out there.

## Layout

```
.
├── index.html               landing page: searchable protocol list, work in progress (generated)
├── build_index.py           builds it from catalogue/ours.tsv + the notes; runs the suites
├── build_site.py            builds every page and assembles the website in _site/
├── .github/workflows/pages.yml   builds and deploys the website on push to main
├── lib/                     shared, protocol-agnostic
│   ├── chemdraw.py          Segment / Construct, duplex rendering, Tm
│   ├── illumina.py          canonical Illumina + NEBNext sequences (verbatim, sourced)
│   ├── nextera.py           Tn5 / Nextera: ME, s5, s7, index + read primers, 9-bp gap
│   ├── rt.py                reverse transcription, oligo-dT, TSO / template switching
│   ├── plasmid.py           GenBank parsing, restriction digest, primer binding, amplicons
│   ├── crispr.py            U6/scaffold landmarks, guide oligos, cloning simulation
│   ├── padlock.py           padlock / MIP capture, gap-fill, circularisation
│   ├── endprep.py           end repair, opposed dA tails, compatible TA ligation
│   ├── hairpin.py           self-pairing oligos and site-specific loop opening
│   ├── dumbbell.py          covalently closed hairpin-ended sequencing templates
│   ├── restriction.py       typed restriction cuts, cohesive-end fill and junctions
│   ├── telomere.py          chromosome 3′ overhangs, repeat phases and affinity capture
│   ├── seqprimers.py        Read 1/2 + index primers, located on a library by computation
│   ├── mdfacts.py           {{= ...}} computed facts inside the Markdown notes
│   ├── page.py              shared stylesheet and page scaffolding
│   └── checks.py            self-test harness + checks true for every protocol
├── ref/                     shared reference material
│   ├── concepts/            chemistry that recurs across protocols, written once
│   │   ├── tn5-tagmentation.md
│   │   ├── template-switching.md
│   │   ├── small-rna-ligation.md
│   │   ├── reverse-transcription.md
│   │   └── padlock-circularization.md
│   ├── scg_lib_structs_style.md   the house style, reverse-engineered
│   ├── SPLiT-seq.html             archived style exemplar
│   └── page_format.css            its stylesheet
├── catalogue/               what exists out there, and what we have tackled
│   ├── scg_lib_structs.tsv  the scraped worklist (protocol x paper), built by tools/
│   ├── ours.tsv             our own coverage -- and the source of the names below
│   └── README.md            the naming scheme and the worklist, with live counts
│
├── <protocol>__<doi>/       one directory per protocol, named for its defining paper
│   ├── 0*.md                reference notes, one per source
│   ├── *.html               the generated page
│   ├── ref/                 protocol-specific reference material, where there is any
│   └── tools/
│       ├── <chem>.py        segment definitions for this chemistry
│       ├── build_page.py    emits the HTML
│       └── selftest.py      shared checks + this protocol's checks
│                            (the seven that exist are listed in catalogue/ours.tsv)
```

## Using it

```sh
# self-tests -- run first, and after any edit. Every suite is found, so nothing here
# needs updating when a protocol is added or a directory is renamed.
for t in */tools/selftest.py */*/tools/selftest.py; do python3 "$t"; done
python3 tools/selftest.py                 # Markdown renderer, computed facts, docs build
python3 tools/lint_pairing.py             # every page: drawn duplex columns must base-pair

# the diagram pages, one build script each (they carry base-by-base figures)
for b in */tools/build_page.py; do python3 "$b"; done

# reading sources for a protocol (details: AGENTS.md, "Reading sources")
python3 tools/get_sources.py "sci-CAR-seq"    # paper, preprint, supplements -> $CHEM_DATA/sources/<slug>/
python3 tools/doctext.py DIR                  # PDF/DOCX/XLSX/HTML -> FILE.txt (get_sources does this)
# every oligo-looking sequence, as written (modifications kept), with its nearby name, the
# lib/ sequences it contains, and its context. Heuristic, tuned for recall: a pointer into
# the document for a person or an LLM to read, never a fact
python3 tools/scrape_primers.py DIR             # or files; --tsv, --known-only, --max-hits
python3 tools/scrape_primers.py DIR --find AGATGTGTATAAGAGACAG  # where does an oligo occur

# the catalogue of what exists (network; cached, so re-running is cheap)
python3 catalogue/tools/fetch_scg_lib_structs.py

# every .md note -> a sibling .html, plus notes.html; then the front page
python3 build_docs.py                     # regenerates all note pages
python3 build_index.py                    # regenerates index.html

# or all of the above plus the deployable website in _site/ (see "Website")
python3 build_site.py
```

Run the self-tests first and after any edit. Each one runs the shared checks plus its own,
so a change to a shared sequence is caught by every protocol that depends on it.

### Keeping the HTML current

Every `.md` in the repo renders to a sibling `.html`, and `notes.html` indexes them all.
One command does the lot, and it only rewrites files whose content actually changed:

```sh
python3 build_docs.py            # render anything that changed
python3 build_docs.py --check    # exit 1 and name the stale pages; nothing written
python3 build_docs.py --force    # re-render everything
python3 build_docs.py --list     # list the notes it would render
python3 build_docs.py crispr-mip__10.1101+2024.03.28.587082/01_crispr-mip.md
```

`--check` is the one to put in CI or a pre-commit hook: it fails if a `.md` was edited
without regenerating its page, so the HTML cannot silently drift from the source.

The renderer is `lib/mdrender.py` — a small dependency-free Markdown subset, because none
of python-markdown, mistune, commonmark or markdown-it is installed and the repo takes no
third-party dependencies. It covers what the notes use: ATX headings with anchors, GFM
tables with alignment, fenced code, blockquotes, nested lists, rules, and inline
code/bold/italic/links/autolinks. Two deliberate departures, both forced by the content:

- **No indented code blocks.** Every code block here is fenced; all 84 four-space indents
  in the notes are list continuations, so an indented line is continuation text.
- **Raw HTML is escaped except for `<br> <sub> <sup> <i> <b> <em> <strong>`.** The notes
  use angle brackets as *metavariables* far more often than as markup — `<cbc>`,
  `<seqid>`, `<align>`, `<inf>` — and those must render literally rather than vanish into
  the DOM.

Links to `.md` files are rewritten to `.html`, so the generated site navigates the same
way the sources do. `tools/selftest.py` pins the supported subset by example and then
renders every real note in the repo, asserting that tags balance, that each source table
and fence becomes exactly one element, and that no markdown is left unconsumed.

No dependencies beyond the Python standard library.

## Website

The generated pages are published as a static site on GitHub Pages. One command builds it:

```sh
python3 build_site.py              # every page, then _site/, then a link check
python3 -m http.server -d _site    # preview at http://localhost:8000/
```

`build_site.py` runs every `*/tools/build_page.py`, then `build_docs.py` and
`build_index.py`, and copies into `_site/` (gitignored) only the public front page, finished
diagram pages, a `.nojekyll`, and the few files of ours those pages deliberately link to.
The main research note for a finished schematic is available through a quiet link on that
schematic page. Other reference notes and work in progress remain in the repository and are
not public site navigation or search content. It never copies third-party material --
nothing from `_data/`, `pdf/`, the catalogue's download cache, archived exemplars such as
`ref/SPLiT-seq.html`, any `ref/`
data file, or a file type `.gitignore` treats as source material. It then checks that every
relative link and anchor in `_site/` resolves, and fails if one does not.

**The front page** (`index.html`) is built from data, not a hand-written list:

- **Published protocols** -- every row of `catalogue/ours.tsv` with `section` = `published`,
  in a compact list led by chemistry filters rather than a name search. The filter data are
  `catalogue/properties.tsv`: sixteen closed facets covering indexing, assay, platform,
  partitioning, input, fragmentation, adapter installation, amplification, topology,
  strand handling, identifiers, reads, selection, conversion and availability. Profiles
  supply shared facts and each row records only its differences. `properties.load()` expands
  and validates the records; a missing protocol, empty facet or unknown value stops the
  build. Options within a facet are OR, while facets combine with AND. An unobtrusive text
  search remains for titles, aliases, papers and descriptions. Filter state is kept in the
  URL, so a result set can be linked.

  A title opens its schematic when one exists; entries without one remain listed by name
  rather than exposing their working notes. The index is JSON embedded in the page and
  filtered by vanilla JavaScript: no library or network, so it works from `file://` too,
  and without JavaScript the full list is still shown. Reference-note bodies are never
  embedded.

Each entry's title is its directory's `protocol` in ours.tsv; its blurb is the first
paragraph of its first note that describes the method (preferring the `## 1. What it is`
section), kept in the search data rather than printed on every entry. A new protocol appears
by name from `ours.tsv`; its title becomes a link when `tools/build_page.py` produces a
finished schematic.

When adding a protocol, add its row to both `catalogue/ours.tsv` and
`catalogue/properties.tsv`. Choose the closest named profile, then override every facet
where the protocol differs. In particular, keep **library amplification** separate from
**how an index is introduced**, and keep **adapter architecture** (TruSeq, Nextera, custom)
separate from **adapter installation chemistry** (ligation, Tn5, PCR, template switching).

**In CI** (`.github/workflows/pages.yml`, on every push to `main` and on demand) the runner
has only what is committed. The self-tests skip what needs a third-party file, and the
diagram pages whose build scripts need one (plasmid maps, for CRISPR-MIP, CRISPR-UMI and the
pooled-screen page) fail with their "missing source" message; the site is built without them
-- the index does not link them, and a short stand-in page at their address explains why,
so a link from another page still lands. Those pages are therefore only complete in a local
build that has the sources.

**Turning it on** (once the repository has a GitHub remote): Settings → Pages → Build and
deployment → Source: **GitHub Actions**. The next push to `main` (or a manual run of the
"Pages" workflow) builds and deploys the site.

## Adding a protocol

1. `mkdir -p <name>/tools <name>/ref`
2. Write the reference notes first, one `.md` per source document, marking evidence as you
   go (see *Evidence marking* below). The notes are where the thinking happens; the page is
   a rendering of conclusions already reached.
3. Write `<name>/tools/<name>.py` defining the construct as a list of `Segment`s. Import
   canonical sequences from `illumina` rather than retyping them — that is the single
   biggest source of avoidable error, and it means a correction propagates everywhere.
4. Put chemistry invariants in the model and rendering functions used by the page. Add a
   self-test only for source transcription or behavior that cannot be enforced by those
   functions; do not restate model constants as checks.
5. Write `<name>/tools/build_page.py` to emit the page from the construct.

Each script starts with the same three-line bootstrap:

```python
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
```

## Conventions

**Drawing.** Duplex ladders draw the bottom strand as the plain *complement* of the top
(not the reverse complement), written left to right and labelled `3'-…-5'`, so that paired
bases share a column. `revcomp()` exists separately, for working out primer sequences.

**Placeholders.** Barcodes, linkers and inserts are stand-ins (`AAAAAAAA`, `LLLL`,
`XXX...XXX`), never complemented as if they were bases — this is a property of the
`Segment`, never of the characters, so `AAAAAAAA` (barcode A) can never become eight
thymines. Their complements render lowercase, which keeps them visually distinct from real
bases and from a dA/dT ligation junction sitting next to them.

**Evidence marking.** Notes mark every claim 🟢 verbatim from the source · 🟡 derived or
inferred · 🔴 not published anywhere. On the page, inference is marked three ways, all
required: a standing caveat in the preamble, an `INFERRED — …` caption on each affected
step, and an `<inf>` element (dotted underline) around the bases themselves. Colour says
what a region *is*; `<inf>` says how well it is *known*. The two are orthogonal and nest
freely.

**Open questions are switches, not prose.** Where a length or structure is genuinely
uncertain, make it a module constant with both options documented
(e.g. `D_ARM_INCLUDES_CCGATCT`). Flip it, rebuild, and every diagram re-aligns — and a
check asserts what the flip should and should not change.

## Protocols

| Directory | Chemistry | Status |
|---|---|---|
| `atrandi-wgs__10.1101+2025.06.20.660799/` | DNA. Semi-permeable capsules, PTA amplification, 4-round split-pool barcoding, NEBNext FS library prep | notes + page built; barcode cassette architecture still modelled, not documented |
| `smart-seq__10.1038+nbt.2282/` | RNA. Original SMART-seq | schematic + shared construct checks |
| `smart-seq2__10.1038+nmeth.2639/` | RNA. SMART-seq2 | schematic + shared construct checks |
| `smart-seq3__10.1038+s41587-020-0497-0/` | RNA. SMART-seq3 | schematic + shared construct checks |
| `smart-seq3xpress__10.1038+s41587-022-01311-4/` | RNA. SMART-seq3xpress | schematic + shared construct checks |
| `flash-seq__10.1038+s41587-022-01312-3/` | RNA. FLASH-seq | schematic + shared construct checks |
| `small-seq__10.1038+nbt.3701/` | RNA. Small-seq: single-cell small RNA / miRNA. TruSeq **Small RNA** adapters, sequential ligation, UMI in the ligated 5' adapter, 5.8S rRNA masking oligo | notes + constructs + page + 167 checks; the SRX index primer is not published, so its arm is bracketed from two published library sizes |
| `lenticrispr-v1-screening__10.1126+science.1247005/` | DNA. lentiCRISPR v1 screening; v2 and the one-/two-PCR readouts have separate pages | notes, plasmid-map provenance, and a focused schematic |
| `crispr-umi-schmierer__10.15252+msb.20177834/` | DNA. A lineage UMI cloned into the guide library &mdash; Schmierer and Michlits | notes + page + checks against the real parent map |
| `astar-seq__10.1101+829960/` | DNA + RNA. C1 chip; Tn5 first, then RT; cDNA biotinylated by PCR and pulled away from the ATAC fragments | reference notes only (status `notes`): no model or page yet |
| `crispr-mip__10.1101+2024.03.28.587082/` | DNA. Padlock/MIP capture with a UMI, replacing the screen readout PCR | notes + full probe from Table S2 + end-to-end checks to the 269 bp library; no page yet |

### A few things recur almost everywhere

Worth reading before adding a protocol, because most new chemistries are a new front end
bolted to one of these:

- **[Tn5 tagmentation](ref/concepts/tn5-tagmentation.md)** — fragmentation and adapter
  addition in one step. The back end of SMART-seq, SPLiT-seq, ATAC-seq and most
  tagmentation-based single-cell methods. Always the same 19-bp mosaic end, the same s5/s7
  entry points, the same 9-bp gap and 72 °C fill-in, and the same rule that only the s5/s7
  heteroduplex amplifies.
- **[Template switching and the TSO](ref/concepts/template-switching.md)** — how a defined
  handle reaches the 5' end of a cDNA without ligation, and the place a UMI can be attached
  before amplification.
- **[Sequential ligation, and masking oligos](ref/concepts/small-rna-ligation.md)** — how a
  22-nt RNA with no cap, no poly(A) and no handle gets two defined ends with no polymerase
  involved: a pre-adenylated 3' adapter, an exonuclease step that destroys the leftovers so
  adapter dimers do not take over, and a 5' adapter that can only ligate to a 5'-phosphate
  (which is what excludes capped mRNA). Plus rRNA depletion by occupying an end rather than
  by pulldown.
- **[Reverse transcription, RNA vs DNA](ref/concepts/reverse-transcription.md)** — the fork
  that determines nearly everything downstream, including why UMIs belong to RNA protocols
  and not to single-cell WGS.
- **[Padlock probes and circularisation](ref/concepts/padlock-circularization.md)** —
  capture by two arms rather than one primer, and the trick of letting an exonuclease throw
  away every failure mode because only the successes are circular.

## Data policy — never check in data

**No sequencing data, count table, or downloaded archive is ever committed to this
repository.** Not a FASTQ, not a subsample, not a "small" test slice. The repo holds *code
and documentation only*; anything that came off an instrument or out of an archive lives
outside it and is reproduced on demand.

What that means in practice:

- Every external dataset gets a **download script** that fetches it from its accession, and
  that script is the artefact we keep. If a result cannot be regenerated from a script in
  this repo plus a public accession, it is not finished.
- Downloads go to a directory named by `$CHEM_DATA` (default `./_data/`), which is
  **gitignored**. The scripts create it; nothing in it is tracked.
- Download scripts must support **subsetting** (`--max-reads`), because the real files here
  run to tens of gigabytes and almost every question can be answered on a slice.
- Scripts record what they fetched — URL, accession, byte count, md5 where the archive
  provides one — into a manifest *next to the data*, not in the repo.
- Derived tables (count matrices, fitted parameters) are also data: they go to `$CHEM_DATA`
  too. Only the script that produces them is tracked.

CRISPR-MIP analysis download scripts live in the companion `chem_crisprmip` repository. A
new analysis keeps its download tooling beside that analysis rather than in this chemistry
catalogue.

### Nor any third-party source document

Papers, vendor manuals, supplementary tables, and anything **derived from them by
extraction** — page-text dumps, figures, gel photographs, spreadsheets, SnapGene exports —
are copyrighted and not ours to redistribute, so none of it is committed either
(`*.pdf`, `*.xls(x)`, `*.doc(x)`, `*.txt`, image formats, `*.dna`). Generated HTML is
ignored for a different reason: it is output, rebuilt by the build scripts, and committing
it would bury a one-line sequence change under thousands of lines of re-aligned diagram.

This costs nothing in checkability, which is the point of the repo: **what we computed from
a source lives in the model and is pinned in the self-test**, so the derived fact is
verifiable and citable without the document. Where a value had to be transcribed, the
transcription itself is quoted once in a comment and the self-test parses that comment and
compares it to what the model derives — see `tools/atrandi.py` and its self-test.

Two conventions make this workable:

- **Every `ref/` directory has a `MANIFEST.md`** naming each file it expects, what it is,
  the exact source to re-obtain it (DOI, supplementary file name, vendor document number
  *and revision*, Addgene ID), and which checks depend on it.
- **Checks that read a third-party file skip instead of failing** when it is absent, naming
  the file and its origin (`checks.Source` / `checks.have`). A fresh clone therefore runs
  green with a loud, specific list of what it could not check — a missing source document
  is an incomplete clone, not a wrong chemistry. Build scripts, which genuinely cannot
  produce a page without their data, exit with that same message instead.

## CRISPR-MIP analysis

This repository keeps the published CRISPR-MIP chemistry and protocol schematic. GC-bias
work, dataset investigations, proposed experiments and live protocol troubleshooting live
in the companion `chem_crisprmip` repository, which may depend on this one. The dependency
never runs in the other direction: this site builds without the analysis checkout.

## A note on what the checks are for

They are regression tests, not discovery. The useful ones assert relationships that *must*
hold if the model is right — for instance, the Atrandi Read-2 layout is derived from the
segment table and then asserted against the offsets hard-coded in Bascet
(`0, 12, 24, 36`, trim `45`), which were arrived at independently from sequencing data.
When a model reproduces a number nobody fed it, that is worth locking down.
