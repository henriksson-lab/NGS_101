# AGENTS.md

Instructions for coding agents working in this repository.

Base-by-base documentation of NGS library chemistries (what the DNA construct looks like at
every protocol step), in the style of Teichlab's scg_lib_structs. `README.md` is the
authoritative, detailed description; this file is the condensed working guide for coding
agents of any kind (Claude Code, Codex, Cursor, Gemini, ...). `CLAUDE.md` only imports it,
so edit this file, not that one.

## Commands

Pure Python 3 standard library — **no third-party dependencies, don't add any** (that is
why `lib/mdrender.py` is a home-grown Markdown renderer). R is used only in `gcbias/R/`.

```sh
# self-tests: run first and after any edit (each suite = shared checks + its own)
for t in */tools/selftest.py */*/tools/selftest.py; do python3 "$t"; done
python3 small-seq__10.1038+nbt.3701/tools/selftest.py   # a single protocol's suite
python3 tools/selftest.py          # Markdown renderer, computed facts, docs build
python3 tools/lint_pairing.py [page.html ...]   # drawn duplex columns must base-pair

# build
for b in */tools/build_page.py; do python3 "$b"; done   # protocol diagram pages
python3 build_docs.py              # every .md -> sibling .html + notes.html (only changed)
python3 build_docs.py --check      # exit 1 if any note page is stale; --force, --list
python3 build_index.py             # index.html (runs the self-tests to get check counts)
python3 build_site.py              # all of the above + _site/ (the website), link-checked
python3 -m http.server -d _site    # preview; CI deploys it (.github/workflows/pages.yml)
```

There is no test framework: a selftest is a script using `checks.Check` (`check(label, got,
want)`), prints PASS/FAIL/SKIP and exits 1 on failure. To run one check, run its suite.

## Reading sources: the source tools

Three standard-library CLIs in `tools/` for turning a paper into reference notes. Every
one has a full `--help`; this is when to use which. Their output is a **pointer into a
document, never a fact** — verify against the source text before anything goes into a
note or a protocol module, and mark it 🟢 only once you have read it there.

| Tool | Does | Typical call |
|---|---|---|
| `tools/get_sources.py` | Fetches everything legitimately reachable for a protocol in the catalogue (or a DOI): the upstream scg_lib_structs page, PMC full text, bioRxiv PDF + supplements, Nature-family supplements from Springer; matches a journal paper to its preprint. Into `$CHEM_DATA/sources/<slug>/` with a `MANIFEST.tsv`; files it could not get are listed there as `(manual)` with their URL. Skips what is on disk. | `python3 tools/get_sources.py "sci-CAR-seq"` · `--doi 10.1101/829960` · `--todo` · `--list` |
| `tools/doctext.py` | PDF / DOCX / XLSX / HTML → a plain-text twin `FILE.txt` beside the original (XLSX as tab-separated rows under `## sheet <name>`). `get_sources.py` runs it on everything it fetches. **Read and cite the `.txt`**: line numbers from the other tools point into it. | `python3 tools/doctext.py DIR` · `FILE --stdout` |
| `tools/scrape_primers.py` | Lists every oligo-like sequence: the bases (PDF line-wrap rejoined), the oligo **as written with its modifications** (`/5Phos/`, `rApp`, `-ddC`, `rGrG+G`, `*`, RNA `rN`), the parsed modifications, a nearby name, the `lib/` sequences it contains on either strand, other lines with the same oligo, and ~400 characters of context. Tuned for **recall**: expect false hits. Files are ordered by their most oligo-like hit and each hit has a `score`; `--max-hits` (default 200/file) cuts the lowest scorers (motif tables, guide libraries). | `python3 tools/scrape_primers.py $CHEM_DATA/sources/<slug>/` · `--tsv` · `--known-only` · `--find SEQ` (where does an oligo occur, either strand) |

What they cannot do: PMC supplements of non-open-access papers sit behind a JavaScript
download gate, and Science / Cell / publisher PDFs behind Cloudflare or a paywall — get
those by hand (the manifest has the URL) and drop them in the same directory. The scraper
finds oligos, not reaction order: the methods still have to be **read**.

## Writing up a protocol (notes first)

1. `python3 tools/get_sources.py "<protocol>"` — fetch; check the manifest for `(manual)`.
2. `python3 tools/scrape_primers.py $CHEM_DATA/sources/<slug>/` — find the oligo table(s).
3. Read the methods in the `.txt` files; for a protocol upstream has drawn, the
   `upstream_*.html.txt` page is a complete second source to check against.
4. Write `<slug>/01_<name>.md` (`<slug>` is the `slug` column of
   `catalogue/scg_lib_structs.tsv`): citation, sources read, what is new, oligos
   **verbatim with modifications**, how they interlock (compute with `lib/`, do not
   eyeball), step by step, final library, read layout, open questions. Example:
   `astar-seq__10.1101+829960/01_astar-seq.md`.
5. Add the directory to `catalogue/ours.tsv` with status `notes`, `section` `published`
   (or `wip` for our own unpublished work; the catalogue self-test pins the wip set, so
   update its `WIP` too) and `modality` `DNA` / `RNA` / `multi`. Re-run
   `python3 catalogue/tools/fetch_scg_lib_structs.py` (cheap: reuses resolved papers),
   `python3 build_docs.py`, and the self-tests. The protocol is listed by name on the public
   front page, but its notes remain maintainer-facing. A `tools/` module, self-test and
   finished schematic page (status `draft` → `documented`) come later; only then does its
   title link to a public page — see README "Adding a protocol".

## Architecture

- **`lib/`** — shared, protocol-agnostic model. `chemdraw.py` (Segment / Construct / Scene,
  complement/revcomp, Tm), canonical sequences in `illumina.py` and `nextera.py`, chemistry
  helpers (`rt.py`, `padlock.py`, `crispr.py`, `plasmid.py`), `seqprimers.py` (locates
  Read 1/2 + index primers on a library by computation), `mdfacts.py` + `mdrender.py`
  (Markdown with computed facts), `page.py` (shared HTML/CSS), `checks.py` (harness +
  `run_common`, checks true for every protocol).
- **`<protocol>__<doi with / → +>/`** — one directory per protocol, named for its defining
  paper. `catalogue/ours.tsv` is the single record mapping protocol → directory, and the
  catalogue self-test rebuilds directory names from it (so renames must go through it).
  Each contains `0*.md` notes, `ref/` (with `MANIFEST.md`), and `tools/` with
  `<chem>.py` (segment definitions), `build_page.py`, `selftest.py` — except status
  `notes` directories, which hold reference notes only. Scripts bootstrap
  `sys.path` with `HERE` and `HERE.parents[1] / "lib"` (see README "Adding a protocol").
- **`catalogue/`** — scraped worklist of protocols from scg_lib_structs (`tools/fetch_scg_lib_structs.py`, network, cached).
- **`gcbias/`** — separate, self-contained analysis (PCR GC bias via lineage UMIs) with
  its own lib, download scripts, tools and R plots; protocol pages never depend on it.
  Workflow is in `gcbias/README.md`.
- **Website** — `build_site.py` builds every page but publishes only the front page,
  finished schematic pages, and files those pages deliberately link to. Reference notes,
  evidence logs, catalogue internals and work in progress stay out of public navigation
  and search; never publish `_data/`, `pdf/`, caches, `ref/` data files or `to_debug/`.
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
- **Import canonical sequences** from `illumina` / `nextera` rather than retyping them.
- **Facts in Markdown come from the model:** `{{= module.EXPR }}` and
  ```` {{oligo: module.fn()}} ```` are expanded by `build_docs.py`; every lib and protocol
  module is in scope, and an unresolvable fact is a build error. Don't hard-code numbers
  in prose that the model can compute.
- **Sequencing primers are declared** by reference (`seqprimers`), and `sp.section()`
  derives where they bind; a declared primer without a site fails the build unless the
  protocol says why.
- **Encode every hand-verified identity as a check** in the protocol's selftest.
- **Drawing convention:** the bottom strand is the plain complement (not revcomp),
  written left to right as `3'-…-5'`. Placeholders (barcodes, linkers, inserts) are never
  complemented as bases — that is a property of the `Segment`; their complements render
  lowercase.
- **Uncertainty is a switch, not prose:** make an uncertain length/structure a module
  constant (e.g. `D_ARM_INCLUDES_CCGATCT`) and add a check of what flipping it changes.
- **Evidence marking** in notes: 🟢 verbatim from source · 🟡 derived/inferred · 🔴 not
  published. On pages, inferred regions need all three: a preamble caveat, an
  `INFERRED — …` step caption, and `<inf>` around the bases.
- **mdrender subset:** no indented code blocks (4-space indents are list continuations);
  raw HTML is escaped except `<br> <sub> <sup> <i> <b> <em> <strong>` because `<cbc>`-style
  metavariables must render literally.

## Never commit

- **Data** — FASTQ/BAM/count tables/derived tables. Downloads go to `$CHEM_DATA` (default
  `./_data/`, gitignored) via a download script that supports `--max-reads` and writes a
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
