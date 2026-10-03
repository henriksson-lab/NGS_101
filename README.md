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

## Layout

```
.
├── lib/                     shared, protocol-agnostic
│   ├── chemdraw.py          Segment / Construct, duplex rendering, Tm
│   ├── illumina.py          canonical Illumina + NEBNext sequences (verbatim, sourced)
│   └── checks.py            self-test harness + checks true for every protocol
├── ref/                     shared reference material
│   ├── scg_lib_structs_style.md   the house style, reverse-engineered
│   ├── SPLiT-seq.html             archived style exemplar
│   └── page_format.css            its stylesheet
├── atrandiWGS/              protocol: Atrandi semi-permeable capsules + PTA
│   ├── 0*.md                reference notes, one per source
│   ├── atrandi_wgs.html     the generated page
│   ├── ref/                 protocol-specific reference (barcode whitelist)
│   └── tools/
│       ├── atrandi.py       segment definitions for this chemistry
│       ├── build_page.py    emits the HTML
│       └── selftest.py      shared checks + this protocol's checks
└── crisprmip/               protocol: to come
```

## Using it

```sh
python3 atrandiWGS/tools/selftest.py     # 88 checks; run this first and after any edit
python3 atrandiWGS/tools/build_page.py   # regenerates atrandi_wgs.html
```

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
| `atrandiWGS/` | Atrandi semi-permeable capsules, PTA amplification, 4-round split-pool barcoding, NEBNext FS library prep | page built; barcode cassette architecture still modelled, not documented |
| `crisprmip/` | — | planned |

### atrandiWGS

Documents the protocol in [Gourlé *et al.*, *bioRxiv*
2025.06.20.660799](https://doi.org/10.1101/2025.06.20.660799) — Atrandi's Single-Microbe
DNA Barcoding Kit with PTA substituted for the kit's MDA (and the debranching step
therefore dropped).

| File | Source |
|---|---|
| `01_barcoding_kit.md` | Atrandi barcoding kit guide, DGPM02323198001 V3 |
| `02_library_prep.md` | Atrandi library-prep guide, DGPM02323206001 V3 |
| `03_our_protocol.md` | the preprint, plus read structure from [Bascet](https://github.com/henriksson-lab/bascet) |
| `04_nebnext_illumina.md` | NEB and Illumina primary documents |
| `05_barcode_cassette_model.md` | the one genuinely inferred part of the construct |
| `06_pta.md` | the PTA literature |

The vendor PDFs the notes cite are gitignored; each note carries its source's document
number and revision so it can be re-obtained.

## A note on what the checks are for

They are regression tests, not discovery. The useful ones assert relationships that *must*
hold if the model is right — for instance, the Atrandi Read-2 layout is derived from the
segment table and then asserted against the offsets hard-coded in Bascet
(`0, 12, 24, 36`, trim `45`), which were arrived at independently from sequencing data.
When a model reproduces a number nobody fed it, that is worth locking down.
