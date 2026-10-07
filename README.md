# NGS protocol chemistry

Base-by-base documentation of sequencing library chemistries — what the DNA construct
actually looks like at every step of a protocol, drawn in the idiom of the Teichmann lab's
[scg_lib_structs](https://teichlab.github.io/scg_lib_structs/) pages.

Each protocol gets a set of reference notes and one generated HTML page. The sequences and
diagrams are **generated from a segment table, not hand-written**, and a self-test pins
every verified fact so a later edit cannot quietly break it.

## Why it is built this way

Hand-aligning a dozen near-identical duplex ladders is exactly the kind of work that goes
wrong silently: a reverse complement written backwards, a primer arrow two columns off, an
annotation row that drifts the moment a segment length changes. Several real errors were
caught this way during the first protocol — a primer line double-counting an overlap
between two adjacent sequences, an adapter strand written 5'→3' where it should have been
3'→5', a junction base on the wrong strand.

So: define the construct once, derive every diagram from it, and assert the things you
checked by hand.

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
smart-seq-family__10.1038+nbt.2282
```

Protocol names collide and get reused; a DOI does not, and every preprint has one from the
day it is posted. `catalogue/ours.tsv` is the single record of which protocol lives in
which directory, and the catalogue self-test rebuilds every directory name from it, so the
name on disk and the paper it claims cannot drift apart. `catalogue/README.md` explains the
scheme; `catalogue/scg_lib_structs.tsv` is the worklist of what is out there.

## Layout

```
.
├── index.html               landing page linking the protocol pages (generated)
├── build_index.py           builds it; check counts come from running the suites
├── lib/                     shared, protocol-agnostic
│   ├── chemdraw.py          Segment / Construct, duplex rendering, Tm
│   ├── illumina.py          canonical Illumina + NEBNext sequences (verbatim, sourced)
│   ├── nextera.py           Tn5 / Nextera: ME, s5, s7, index + read primers, 9-bp gap
│   ├── rt.py                reverse transcription, oligo-dT, TSO / template switching
│   ├── plasmid.py           GenBank parsing, restriction digest, primer binding, amplicons
│   ├── crispr.py            U6/scaffold landmarks, guide oligos, cloning simulation
│   ├── padlock.py           padlock / MIP capture, gap-fill, circularisation
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
├── gcbias/                  ALL the GC-bias work, self-contained
│   ├── datasets/            which public data exists per paper, and if it is reusable
│   ├── download/            download scripts (data is NEVER committed)
│   ├── lib/umimodel.py      lineage-UMI molecule counting: collisions, dropout, Chao1
│   ├── tools/               extract, diagnostics, model fit, simulation, controls
│   └── R/                   plots and the GAM non-linearity test
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
└── to_debug/                a live problem: ira1.md + check.py reproduce the analysis
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

# the catalogue of what exists (network; cached, so re-running is cheap)
python3 catalogue/tools/fetch_scg_lib_structs.py

# every .md note -> a sibling .html, plus notes.html; then the front page
python3 build_docs.py                     # regenerates all note pages
python3 build_index.py                    # regenerates index.html
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
python3 build_docs.py gcbias/README.md   # just these files
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

## Adding a protocol

1. `mkdir -p <name>/tools <name>/ref`
2. Write the reference notes first, one `.md` per source document, marking evidence as you
   go (see *Evidence marking* below). The notes are where the thinking happens; the page is
   a rendering of conclusions already reached.
3. Write `<name>/tools/<name>.py` defining the construct as a list of `Segment`s. Import
   canonical sequences from `illumina` rather than retyping them — that is the single
   biggest source of avoidable error, and it means a correction propagates everywhere.
4. Write `<name>/tools/selftest.py`: call `run_common(check)` from `lib.checks`, then add
   the checks specific to the chemistry. **Encode every identity you verified by hand.**
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
| `smart-seq-family__10.1038+nbt.2282/` | RNA. Template switching, full-length cDNA, Nextera tagmentation | notes + constructs + checks; no page yet |
| `small-seq__10.1038+nbt.3701/` | RNA. Small-seq: single-cell small RNA / miRNA. TruSeq **Small RNA** adapters, sequential ligation, UMI in the ligated 5' adapter, 5.8S rRNA masking oligo | notes + constructs + page + 167 checks; the SRX index primer is not published, so its arm is bracketed from two published library sizes |
| `lenticrispr-gecko-screen__10.1126+science.1247005/` | DNA. Pooled CRISPR screening: lentiCRISPR/GeCKO, Broad GPP, and the single-cell screens | notes, plasmid maps, page, and primer/amplicon checks against real maps |
| `crispr-umi-schmierer__10.15252+msb.20177834/` | DNA. A lineage UMI cloned into the guide library &mdash; Schmierer and Michlits | notes + page + checks against the real parent map |
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

Download scripts live in [`gcbias/download/`](gcbias/download/) for the GC-bias work; a new
analysis adds its own `download/` beside its tools.

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

## GC-bias analysis

[`gcbias/`](gcbias/README.md) is a separate track from the protocol pages: not chemistry but
an analysis, and **self-contained** — its own `datasets/` catalogue, its own
`lib/umimodel.py`, its own download scripts, tools and R plots. Nothing in it is needed to
build a protocol page, and nothing in a protocol page depends on it.

It was started to answer a reviewer question about PCR GC bias in CRISPR-MIP. Two things in
it are worth knowing even from outside:

- **The obvious statistic is wrong.** Reads per distinct UMI invents a GC trend where none
  exists — by up to +1.4 log2 per unit GC in simulation — because distinct-UMI counts
  saturate and saturation tracks abundance. `lib/umimodel.py` models the collisions and the
  unseen molecules so that "expected reads" is a prediction rather than a restatement.
- **Both published CRISPR-UMI datasets fail, for opposite reasons** — one UMI is 82 %
  saturated, the other was never deposited. `datasets/` records which public data exists per
  paper, whether it is reusable, and the criteria a dataset has to meet, so the search does
  not have to be redone.

On the Schmierer plasmid input the real effect is small, strongly non-linear, and
concentrated above ~65 % GC; on genomic DNA it is roughly 1.7× larger.

## A note on what the checks are for

They are regression tests, not discovery. The useful ones assert relationships that *must*
hold if the model is right — for instance, the Atrandi Read-2 layout is derived from the
segment table and then asserted against the offsets hard-coded in Bascet
(`0, 12, 24, 36`, trim `45`), which were arrived at independently from sequencing data.
When a model reproduces a number nobody fed it, that is worth locking down.
