#!/usr/bin/env python3
"""Generate index.html -- the landing page linking the protocol pages together.

Check counts are read by actually running each protocol's self-test, so the index cannot
drift from reality: if a suite is failing, the page says so.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))

import catalogue as cat  # noqa: E402
from page import caveat, head, info, legend  # noqa: E402

OUT = ROOT / "index.html"

PROTOCOLS = [
    dict(dir="atrandi-wgs__10.1101+2025.06.20.660799", name="Atrandi SPC + PTA",
         page="atrandi-wgs__10.1101+2025.06.20.660799/atrandi_wgs.html", kind="DNA",
         blurb="Single-microbe whole-genome sequencing. Cells are held in semi-permeable "
               "capsules through lysis and amplification, then given a cell barcode by "
               "four rounds of split-pool ligation.",
         chem=["semi-permeable capsules", "PTA amplification", "split-pool barcoding",
               "NEBNext FS library prep"],
         note="The barcode cassette architecture is modelled, not documented &mdash; "
              "Atrandi publish the barcode sequences but not the oligo design."),
    dict(dir="florian-pta-rnaseq", name="florian-pta-rnaseq",
         page="florian-pta-rnaseq/florian-PTA-rnaseq.html", kind="RNA",
         blurb="An attempt at single-cell RNA-seq read out through Atrandi's SPC/PTA "
               "workflow: a Smart-seq3xpress front end for the 5' UMI, and a back end "
               "where all four rounds of split-pool barcoding are replaced by one "
               "pre-annealed duplex.",
         chem=["template switching", "oligo-dT priming", "PTA amplification",
               "dA-tail adapter ligation"],
         note="Built from the ordered oligos only &mdash; the PTA chemistry and the "
              "indexing primers are not specified yet, so the final library is "
              "predicted, and three design questions are still open."),
    dict(dir="smart-seq-family__10.1038+nbt.2282", name="SMART-seq family",
         page="smart-seq-family__10.1038+nbt.2282/smartseq.html", kind="RNA",
         blurb="Full-length single-cell mRNA. SMART-seq and SMART-seq2 put the same handle "
               "on both ends of the cDNA; SMART-seq3 deliberately breaks that symmetry to "
               "get UMI counting as well as coverage.",
         chem=["template switching", "oligo-dT priming", "Tn5 tagmentation"],
         note=""),
    dict(dir="small-seq__10.1038+nbt.3701", name="Small-seq", page="small-seq__10.1038+nbt.3701/smallseq.html", kind="RNA",
         blurb="Single-cell small-RNA / miRNA sequencing. Two adapters are ligated on in "
               "sequence -- no polymerase touches either end -- and the UMI rides inside "
               "the ligated 5' adapter, where it is the first thing sequenced.",
         chem=["TruSeq Small RNA adapters", "sequential ligation", "UMI in the 5' adapter",
               "rRNA masking oligo"],
         note="The SRX index primers are in a supplementary table we do not have, so "
              "everything right of the 3' adapter is modelled &mdash; its length "
              "bracketed by two independent published library sizes."),
    dict(dir="lenticrispr-gecko-screen__10.1126+science.1247005", name="Pooled CRISPR screening",
         page="lenticrispr-gecko-screen__10.1126+science.1247005/crisprscreen.html", kind="DNA",
         blurb="lentiCRISPR / GeCKO, Broad GPP and the single-cell screens. Mostly a study "
               "of why published readout primers are not interchangeable between vectors.",
         chem=["BsmBI guide cloning", "sgRNA scaffolds", "staggered readout PCR"],
         note="Cloning and readout both drawn from real Addgene maps &mdash; including "
              "why KERMIT, BEAKER and the Joung pair are not interchangeable."),
    dict(dir="crispr-umi-schmierer__10.15252+msb.20177834", name="CRISPR-UMI (Schmierer)",
         page="crispr-umi-schmierer__10.15252+msb.20177834/crisprumi.html", kind="DNA",
         blurb="Two 2017 methods that put a lineage label in the guide library itself, so "
               "counting labels counts clones rather than reads &mdash; Schmierer (6-bp RSL, "
               "AU-flip vector) and Michlits (10-nt barcode, PacI enrichment).",
         chem=["UMI in the library", "AU-flip and F+E scaffolds", "BbsI Golden Gate",
               "PacI pre-PCR enrichment", "lineage dropout analysis"],
         note="Three methods share the name &ldquo;CRISPR-UMI&rdquo;. These two change the "
              "vector; CRISPR-MIP does not."),
    dict(dir="crispr-mip__10.1101+2024.03.28.587082", name="CRISPR-MIP",
         page="crispr-mip__10.1101+2024.03.28.587082/crisprmip.html", kind="DNA",
         blurb="A padlock probe carrying a UMI replaces the readout PCR of a pooled "
               "screen, so the result counts molecules instead of reads. The UMI is on the "
               "probe, not in the library.",
         chem=["padlock capture", "circularisation", "exonuclease selection", "UMI counting"],
         note=""),
]

CONCEPTS = [
    ("Tn5 tagmentation", "ref/concepts/tn5-tagmentation.md",
     "Fragmentation and adapter addition in one step. The back end of SMART-seq, SPLiT-seq, "
     "ATAC-seq and most tagmentation methods &mdash; always the same 19-bp mosaic end, the "
     "same 9-bp gap and 72&nbsp;&deg;C fill-in, and the same rule that only the s5/s7 "
     "heteroduplex amplifies."),
    ("Template switching and the TSO", "ref/concepts/template-switching.md",
     "How a defined handle reaches the 5' end of a cDNA without ligation, by exploiting the "
     "untemplated C's reverse transcriptase adds when it runs off the template &mdash; and "
     "the place a UMI can be attached before amplification."),
    ("Sequential ligation, and masking oligos", "ref/concepts/small-rna-ligation.md",
     "How a 22-nt RNA with no cap, no poly(A) and no handle gets two defined ends with no "
     "polymerase involved &mdash; a pre-adenylated 3' adapter, an exonuclease step that "
     "destroys the leftovers, and a 5' adapter that can only ligate to a 5'-phosphate. "
     "Plus rRNA depletion by occupying an end rather than by pulldown."),
    ("Reverse transcription, RNA vs DNA", "ref/concepts/reverse-transcription.md",
     "The fork that determines nearly everything downstream, including why UMIs belong to "
     "RNA protocols and not to single-cell WGS."),
    ("Padlock probes and circularisation", "ref/concepts/padlock-circularization.md",
     "Capture by two arms rather than one primer, and the trick of letting an exonuclease "
     "throw away every failure mode because only the successes are circular."),
]


def checks(d: str) -> tuple[int, bool]:
    """Run a protocol's self-test and report (count, passed)."""
    st = ROOT / d / "tools" / "selftest.py"
    if not st.exists():
        return 0, False
    r = subprocess.run([sys.executable, str(st)], capture_output=True, text=True)
    m = re.search(r"all (\d+) checks passed", r.stdout)
    if m:
        return int(m.group(1)), True
    m = re.search(r"(\d+) of (\d+) checks FAILED", r.stdout)
    return (int(m.group(2)) if m else 0), False


EXTRA_CSS = """<style>
.cards { display:grid; gap:14px; grid-template-columns:repeat(auto-fit,minmax(330px,1fr));
         margin:1.2em 0; align-items:start; }
.card { border:1px solid var(--rule); border-radius:6px; padding:16px 18px;
        background:var(--surface); display:flex; flex-direction:column; gap:8px; }
.card h3 { margin:0; font-size:1.08rem; }
.card h3::before { content:none; }
.card p { margin:0; font-size:.95rem; }
.kind { font-size:.72rem; letter-spacing:.08em; text-transform:uppercase;
        color:var(--ink-muted); font-weight:600; }
.tags { display:flex; flex-wrap:wrap; gap:5px; }
.tag { font-size:.76rem; padding:2px 8px; border-radius:999px;
       background:var(--surface-2); color:var(--ink-muted); }
.meta { font-size:.82rem; color:var(--ink-muted); font-variant-numeric:tabular-nums;
        margin-top:auto; padding-top:6px; border-top:1px solid var(--rule); }
.go { font-weight:600; }
.nopage { color:var(--ink-muted); font-style:italic; }
.concept { margin:1.1em 0; }
.concept h3 { margin:0 0 .2em; font-size:1rem; }
.concept h3::before { content:none; }
.concept p { margin:0; max-width:76ch; font-size:.95rem; }
</style>"""


def card(p: dict) -> str:
    n, ok = checks(p["dir"])
    tags = "".join(f'<span class="tag">{c}</span>' for c in p["chem"])
    link = (f'<a class="go" href="{p["page"]}">Open the chemistry page &rarr;</a>'
            if p["page"] else '<span class="nopage">No diagram page yet</span>')
    note = f'<p style="font-size:.88rem;color:var(--ink-muted)">{p["note"]}</p>' if p["note"] else ""
    status = f"{n} checks passing" if ok else (f"&#9888; {n} checks, suite FAILING" if n else "no suite")
    return f"""<div class="card">
<span class="kind">{p['kind']}</span>
<h3>{p['name']}</h3>
<p>{p['blurb']}</p>
<div class="tags">{tags}</div>
{note}
{link}
<div class="meta"><code>{p['dir']}/</code> &middot; {status}</div>
</div>"""


# Suites that are not a protocol page but still gate the build.
EXTRA_SUITES = ["catalogue", "gcbias"]


def main() -> None:
    total = (sum(checks(p["dir"])[0] for p in PROTOCOLS)
             + sum(checks(d)[0] for d in EXTRA_SUITES))
    n_notes = len([q for q in ROOT.rglob("*.md")
                   if not {".git", "_data", "pdf", "__pycache__"}
                   & set(q.relative_to(ROOT).parts)])
    cards = "\n".join(card(p) for p in PROTOCOLS)
    concepts = "\n".join(
        f'<div class="concept"><h3>{t}</h3><p>{d} <span style="font-size:.85rem">'
        f'<code>{f}</code></span></p></div>' for t, f, d in CONCEPTS)
    body = f"""<div class="wrap">
<h1>NGS protocol chemistry</h1>

{info("""What the DNA actually looks like at every step of a sequencing protocol, drawn
base-by-base in the idiom of the Teichmann lab's
<a href="https://teichlab.github.io/scg_lib_structs/">scg_lib_structs</a> pages.""")}

{legend(f"""<b>Everything here is generated, not transcribed.</b> Each protocol's construct
is defined once as a table of segments; the diagrams, oligo lists and final structures are
derived from it, and a self-test asserts the facts that were checked by hand &mdash;
<b>{total} checks</b> across the protocols. Where a real plasmid exists, primer sites
and amplicon sizes are computed against the actual map rather than quoted from a paper.""")}

<h2>Protocols</h2>
<div class="cards">
{cards}
</div>

<h2>The catalogue &mdash; what is out there</h2>
{info(f"""Documenting a protocol starts with knowing it exists. <b>{cat.n_protocols()}
protocols</b> and <b>{cat.n_papers()} papers</b> are catalogued in
<code>catalogue/scg_lib_structs.tsv</code>, scraped from
<a href="https://github.com/Teichlab/scg_lib_structs">Teichlab/scg_lib_structs</a> and
resolved against Crossref and NCBI &mdash; {cat.n_documented()} with a drawn page upstream
to check ourselves against, {cat.n_todo()} named but never drawn. One row per
(protocol, paper), because a method can have several papers and a paper can define several
methods. <a href="catalogue/README.html"><b>read the worklist &rarr;</b></a>""")}

<h2>Chemistry that recurs</h2>
{info("""Most new protocols are a new front end bolted onto an old back end. These four
pieces turn up again and again, so they are written once and shared rather than repeated
per protocol.""")}
{concepts}

<h2>Written notes</h2>
{info(f"""Alongside the diagram pages, the repository carries {n_notes} Markdown notes
&mdash; protocol readings, dataset inventories, the enzyme/buffer comparison, and working
debug logs. All of them are rendered to HTML by <code>build_docs.py</code>:
<a href="notes.html"><b>browse all notes &rarr;</b></a>""")}

<h2>How it is built</h2>
{info("""<code>lib/</code> holds the protocol-agnostic machinery: a renderer for
character-aligned duplex diagrams, the canonical Illumina and NEBNext sequences, Tn5 and
reverse-transcription building blocks, a GenBank parser with restriction digestion and
circular-aware PCR prediction, and the padlock capture model. Each protocol directory holds
its own notes, its segment definitions and its self-test, which runs the shared checks plus
its own. The Markdown notes go through <code>lib/mdrender.py</code>, a small
dependency-free renderer, driven by <code>build_docs.py</code>; run
<code>python3 build_docs.py --check</code> to find out whether any page has drifted from
its source.""")}

{caveat("""<b>On trusting this.</b> The checks are regression tests, not proof. The useful
ones assert relationships that <i>must</i> hold if a model is right &mdash; for instance the
Atrandi read layout is derived from the segment table and then asserted against offsets
hard-coded in a demultiplexer that was written from sequencing data. When a model reproduces
a number nobody fed it, that is worth locking down. Anything inferred rather than documented
is marked as such on the page it appears on.""")}
</div>"""
    OUT.write_text(head("NGS Protocol Chemistry") + "\n" + EXTRA_CSS + "\n" + body,
                   encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)  {total} checks total")


if __name__ == "__main__":
    main()
